#!/usr/bin/env python3
"""The leading-order chains of lo.py in TeNPy (charge-conserving DMRG).
Run with ~/.venvs/scp-tenpy/bin/python.

Charges per local state: 2S^z (, 2I^z for W) (, N_hole for mobile sites).
The sector is fixed by the initial product state, so no penalties.
usage: tp_chain.py M mode L_bulk target t_over_J0 chi [sweeps [p11]]
  mode: nohole | hole | mob | mobK (K >= 2 mobile holes on L_bulk + K sites);
  target: comma list 2S^z[,2I^z]; p11: the (1,1) link pricing of lo.py
  (zero | core | upper | lower; default zero)
Prints E/J0, the S^z (and I^z) profile, and for mob the hole density and
the hole-centred correlation C(d)."""
import sys, time, warnings
import numpy as np
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here)
import graphm as GM, lo
import tenpy
import tenpy.linalg.np_conserved as npc
from tenpy.networks.site import Site
from tenpy.networks.mps import MPS
from tenpy.models.lattice import Lattice
from tenpy.models.model import CouplingMPOModel, CouplingModel, MPOModel
from tenpy.algorithms import dmrg
warnings.filterwarnings("ignore")
tenpy.tools.misc.setup_logging(to_stdout=None)

def make_site(L, mobile):
    M = L.M
    names = ['2Sz'] + (['2Iz'] if M == 4 else []) + (['Nh'] if mobile else [])
    chinfo = npc.ChargeInfo([1]*len(names), names)
    q = []
    for k in range(L.d):
        row = list(L.lab[k])
        if mobile: row.append(int(L.N[k] == 1) if L.kind == 'mob' else 0)
        q.append(row)
    leg = npc.LegCharge.from_qflat(chinfo, np.array(q, dtype=int))
    site = Site(leg, [str(k) for k in range(L.d)], sort_charge=False)
    site.add_op('Sz', L.Sz)
    if M == 4: site.add_op('Iz', L.Iz)
    if L.kind == 'mob': site.add_op('nh', L.nh)
    return site

class Chain(CouplingModel, MPOModel):
    def __init__(self, locs, P, g, t, J0, mobile):
        sites = [make_site(L, mobile) for L in locs]
        # one shared charge info is required: rebuild legs on a common chinfo
        lat = Lattice([1], sites, bc='open', bc_MPS='finite')
        CouplingModel.__init__(self, lat)
        n = len(locs)
        for i in range(n - 1):
            terms = lo.link_terms(locs[i], locs[i + 1], P, g) + lo.hop_terms(locs[i], locs[i + 1], t)
            for k, (c, A, B) in enumerate(terms):
                na, nb = f"A{i}_{k}", f"B{i}_{k}"
                sites[i].add_op(na, A) if na not in sites[i].opnames else None
                sites[i + 1].add_op(nb, B) if nb not in sites[i + 1].opnames else None
                self.add_coupling_term(c/J0, i, i + 1, na, nb)
        MPOModel.__init__(self, lat, self.calc_H_MPO())

def run(M, mode, Lb, target, tJ, chi, sweeps=20):
    P = GM.P_W if M == 4 else GM.P_2
    g = 0.1*GM.DELTA; J0 = 2*g*g/GM.DELTA
    Lh = Lb//2
    nholes = 0
    if mode == 'nohole': kinds = ['end'] + ['full']*Lb + ['end']
    elif mode == 'hole': kinds = ['end'] + ['full']*Lh + ['hole'] + ['full']*Lh + ['end']
    elif mode.startswith('pair'):
        dd = int(mode[4:]); kinds = ['end'] + ['full']*(Lb + 2) + ['end']
        c0 = len(kinds)//2 - dd//2; kinds[c0] = 'hole'; kinds[c0 + dd] = 'hole'
    elif mode.startswith('mob') and mode != 'mob':
        nholes = int(mode[3:]); kinds = ['end'] + ['mob']*(Lb + nholes) + ['end']
    else:
        kinds = ['end'] + ['mob']*(Lb + 1) + ['end']; nholes = 1
    locs = [lo.Local(M, k) for k in kinds]
    mobile = mode.startswith('mob')
    model = Chain(locs, P, g, tJ*J0, J0, mobile)
    sites = model.lat.mps_sites()
    # product state in the target sector (greedy, as lo_dmrg.py)
    L = len(locs); sgn = 1; init = []
    for Lc in locs:
        w = np.array([(1 if Lc.kind in ('end', 'hole') else 2)*sgn]*Lc.lab.shape[1]); sgn = -sgn
        cost = np.abs(Lc.lab - w[None, :]).sum(1) + (0 if Lc.kind != 'mob' else 10*(Lc.N != 2))
        init.append(int(np.argmin(cost)))
    if mobile:
        spots = [L//2] if nholes == 1 else [(k + 1)*L//(nholes + 1) for k in range(nholes)]
        for mid in spots: init[mid] = [k for k in range(locs[mid].d) if locs[mid].N[k] == 1][0]
    def tot(init): return np.sum([locs[i].lab[init[i]] for i in range(L)], 0)
    for _ in range(4*L):
        diff = np.array(target) - tot(init)
        if not np.any(diff): break
        best = None
        for i in range(L):
            Lc = locs[i]
            for k in range(Lc.d):
                if mobile and Lc.kind == 'mob' and Lc.N[k] != Lc.N[init[i]]: continue
                nd = np.abs(diff - (Lc.lab[k] - Lc.lab[init[i]])).sum()
                if best is None or nd < best[0]: best = (nd, i, k)
        init[best[1]] = best[2]
    assert not np.any(np.array(target) - tot(init)), "sector not reached"
    psi = MPS.from_product_state(sites, [str(k) for k in init], bc='finite')
    opts = {'max_N_sites_per_ring': None, 'max_trunc_err': 1.0, 'mixer': True, 'max_E_err': 1e-10, 'max_S_err': 1e-7, 'min_sweeps': 6, 'max_sweeps': sweeps,
            'trunc_params': {'chi_max': chi, 'svd_min': 1e-12}}
    info = dmrg.run(psi, model, opts)
    psi.canonical_form()
    E = info['E']
    out = {'E': E, 'Sz': psi.expectation_value('Sz'), 'trunc': float(np.max(info['sweep_statistics']['max_trunc_err'][-1:]))}
    if M == 4: out['Iz'] = psi.expectation_value('Iz')
    if mobile:
        nh = np.zeros(L); C = np.zeros(L)
        idx = [i for i in range(L) if locs[i].kind == 'mob']
        nh[idx] = psi.expectation_value('nh', sites=idx)
        for i in idx:
            row = psi.correlation_function('nh', 'Sz', sites1=[i], sites2=list(range(L)))[0]
            for j in range(L): C[abs(j - i)] += row[j]
        out['nh'] = nh; out['C'] = C
        if nholes >= 2:
            Pr = np.zeros(L)
            for i in idx:
                js = [j for j in idx if j > i]
                if not js: continue
                row = psi.correlation_function('nh', 'nh', sites1=[i], sites2=js)[0]
                for j, v in zip(js, row): Pr[j - i] += v
            out['Pr'] = Pr
    return out, locs

if __name__ == '__main__':
    t0 = time.time()
    M, mode, Lb = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
    target = tuple(int(x) for x in sys.argv[4].split(','))
    tJ, chi = float(sys.argv[5]), int(sys.argv[6])
    sw = int(sys.argv[7]) if len(sys.argv) > 7 else 20
    if len(sys.argv) > 8: lo.P11 = sys.argv[8]
    out, locs = run(M, mode, Lb, target, tJ, chi, sw)
    line = f"M={M} mode={mode} Lb={Lb} target={target} t/J0={tJ} chi={chi} p11={lo.P11} tenpy E/J0={out['E']:.10f} sumSz={out['Sz'].sum():.5f} trunc={out['trunc']:.1e}"
    if M == 4: line += f" sumIz={out['Iz'].sum():.5f}"
    if mode.startswith('mob'): line += f" holes={out['nh'].sum():.5f}"
    print(line + f" [{time.time()-t0:.0f}s]")
    print("Sz " + " ".join(f"{x:+.5f}" for x in out['Sz']))
    if M == 4: print("Iz " + " ".join(f"{x:+.5f}" for x in out['Iz']))
    if mode.startswith('mob'):
        print("nh " + " ".join(f"{x:+.5f}" for x in out['nh']))
        print("C " + " ".join(f"{x:+.5f}" for x in out['C']))
    if 'Pr' in out:
        print("Pr " + " ".join(f"{x:+.5f}" for x in out['Pr']))
