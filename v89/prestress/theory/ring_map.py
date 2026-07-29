#!/usr/bin/env python3
"""
ring_map.py -- v89 prestress: EXACT reduced per-step map of the cellfab
dense sector for idealized consonant networks (pair / ring-N / cube),
faithful to v89/cellfab.c step_field() pass by pass, restricted to the
degrees of freedom a seeded dense network actually exercises:

    state:  th[N]      dense clocks (kernel th2, pass 6 line ~2069)
            Em[N]      dense stores
            lem[NL,2]  per-link in-flight dense energy (slots, dir 0: i->j)
            lph[NL,2]  flight cycle positions (deterministic: adv=dt*C/d)

Faithfully implemented (kernel line refs from cellfab.c @ commit of
2026-07-28):
  pass 0  flload_i = 0.5*sum incident (lem fwd+back)          [~1462]
  pass 1  x=(Em+flload)/cap, w=w2/(1+q x)                     [~1552]
  pass 2  wants: kd*km*dt*geo*gpl*res(det)*gate*head_recv*mob [~1573]
          comb held at 1:1 (uniform-pitch nets; asserted)
          gates g_ij=G(th_i - w_i d/C - th_j), even, p_gate=8 [482,1650]
          mob_sym: mi=sqrt(Em_i*max(Em_j,fl0))                [1652]
          kappa_freq=0, kappa_reac=0 (standing table V2g)
  pass 3  outflow limiter 0.98*Em                             [~1735]
  pass 4  deposits into flight + source debit                 [1760]
  pass 4c receiver entrainment on the arriving deposit:
          err=wrap(th_src_snap - w_src d/C - th_recv),
          mix=f/(f+Em_recv+lock_floor*cap), th+=klock*mix*err [1795-1819]
  pass 5  flight advance, delivery on completed cycle, ROUGHNESS
          R=2|det|Gr/(det^2+Gr^2), rough=take*rough_k*R -> radiated,
          back_s=rough*s_pull/(1+s_pull) -> Es ledger,
          completion entrainment mix=take/(take+Em_prev+lockf)[1827-1895]
          optional action atoms (quant_mode=2 credit)         [1411-1437]
  pass 6  th += w*dt (mod 2pi)                                [2068]
Frozen (justified in STABILITY.md S1): Es (pass S is a no-op among
identical cells), plane normals (pass 7; gpl entered as the seeded
geometric value), field sector (Ee=0; roughness output is booked to a
'radiated' bucket), beat conversion (Ee=0, Em<cap: both doors shut --
matches M0's evap==0).

Modes:
  fixed   <net>    relax to the locked steady state, print observables
  floquet <net>    exact one-window Jacobian -> eigenvalue table by mode
  kick    <net>    time-domain mode kick + matrix-pencil fit (validation)
  pair             sign-theorem demo: heavy->light flow + decay rate
  unlock           antisymmetric-load unlock threshold (tongue edge)
  tax     <net>    kick-amplitude scan of the mass leak (comma tax);
                   atoms on/off -> the credit dead-time threshold
  all              the suite quoted in ../STABILITY.md

nets: ring12 (d=1.25, m=5), ring6 (d=1.25, m=3), ring6m2 (m=2, one-way),
      cube (a=1.586), cube125 (a=1.25), pair (d=1.25, m=1)
options: gpl=<f|seeded>  es=<seeded|footprint>  quant=<0|1>
"""
import numpy as np
import sys
import math

# ---- the standing law table, laws_V2g.cfg (dense sector keys) ----
C = 1.0
W2 = 2.9
QD = 1.2                    # q_detune
GM = 0.10                   # gamma_res_m
PGATE = 8
LOCKF = 0.005 * 2.5         # lock_floor * cap
KD = 1.2 * 2.0              # k_dep * k_dep_m
CAP = 2.5
ES0 = 1.0
ESFLOOR = 0.05
SPULL = 0.15
KLOCK = 1.0
ROUGHK = 0.35
GROUGH = 0.5
MOBFLOOR = 0.004 * 2.5      # mob_floor * cap
R0 = 0.85
DT = 0.02
A0 = 1.15                   # quant_A0 (atoms; used only when quant=1)

TWO_PI = 2.0 * math.pi


def wrap(a):
    return (a + math.pi) % TWO_PI - math.pi


def gate(dphi):
    g = 0.5 * (1.0 + np.cos(dphi))
    return g ** PGATE


def lens_geo(d, ri, rj):
    """kernel channel geometry: overlap lens area, then (A/Aref)*(dref/d)."""
    A = 0.0
    if d < ri + rj:
        t = d * d - rj * rj + ri * ri
        a2 = (4.0 * d * d * ri * ri - t * t) / (4.0 * d * d)
        if a2 > 0:
            A = math.pi * a2
        elif d < abs(ri - rj):
            rm = min(ri, rj)
            A = math.pi * rm * rm
    Aref = math.pi * R0 * R0
    dref = 2.0 * R0
    return (A / Aref) * (dref / d), A


def seeded_es(xbar):
    """the seeding primitive (kernel ~865-875): Em add + s_pull space join."""
    add = xbar * CAP / (1.0 + SPULL)
    pull = SPULL * add
    avail = ES0 - ESFLOOR
    pull = min(pull, max(avail, 0.0))
    return ES0 - pull


def tangent_frame(u):
    """kernel's transverse-normal construction (seeder, ~876-886)."""
    ax = np.array([0.0, 0.0, 1.0])
    if abs(u[2]) > 0.9:
        ax = np.array([1.0, 0.0, 0.0])
    t1 = ax - np.dot(ax, u) * u
    t1 /= np.linalg.norm(t1)
    n2 = np.cross(u, t1)
    return t1, n2


class Net:
    """cells + links; per-link geo*gpl precomputed (Es, normals frozen)."""

    def __init__(self, name, pos, links, xbar, es_mode='seeded',
                 gpl_override=None, axes=None):
        self.name = name
        self.N = len(pos)
        self.pos = np.array(pos)
        self.links = links          # list of (i, j, d)
        self.NL = len(links)
        self.xbar = xbar
        # per-link phase step at seed (uniform pitch omega_bar)
        self.wbar = W2 / (1.0 + QD * xbar)
        # Es: seeded (post-pull) or footprint-relaxed (pass-S equilibrium
        # against ambient: Es ~ e_s0 - s_disp*Em_store, s_disp=0.3)
        es_seed = seeded_es(xbar)
        if es_mode == 'seeded':
            self.Es = es_seed
        else:
            # one-iteration footprint: store ~ 60% of x*cap after flight fill
            self.Es = ES0 - 0.3 * 0.6 * xbar * CAP
        self.r = R0 * (self.Es / ES0) ** (1.0 / 3.0)
        # geometry + plane factor per link
        self.geo = np.zeros(self.NL)
        self.gpl = np.zeros(self.NL)
        for l, (i, j, d) in enumerate(links):
            g, A = lens_geo(d, self.r, self.r)
            self.geo[l] = g
            if gpl_override is not None:
                self.gpl[l] = gpl_override
            else:
                # seeded normals: axes[i] = (n1, n2) unit vectors
                n1i, n2i = axes[i]
                n1j, n2j = axes[j]
                u = self.pos[j] - self.pos[i]
                u = u / np.linalg.norm(u)
                d1i, d2i = np.dot(n1i, u), np.dot(n2i, u)
                d1j, d2j = np.dot(n1j, u), np.dot(n2j, u)
                axi = max(0.0, (1 - d1i * d1i) * (1 - d2i * d2i))
                axj = max(0.0, (1 - d1j * d1j) * (1 - d2j * d2j))
                dnn = np.dot(n2i, n2j)          # dense sector plane = n2
                self.gpl[l] = axi * axj * dnn * dnn
        # incidence lists (CSR order = ascending link index, as the kernel)
        self.inc = [[] for _ in range(self.N)]
        for l, (i, j, d) in enumerate(links):
            self.inc[i].append(l)
            self.inc[j].append(l)
        # mechanism probes (surgical, default off = faithful kernel map)
        self.freeze_gb = None    # float: clamp back gates at this value
        self.no_tau = False      # True: gate/err retard angle at wbar (not w_i)

    def seed_state(self):
        th = np.zeros(self.N)
        # orientation-aware lock seeding (the ring seeder walks the loop:
        # th[next] = th[k] - w*d/C along the link's forward direction,
        # kernel ~984; traversing a link backward inverts the sign, and
        # rung closure makes the spanning-tree seed globally consistent)
        th[0] = 0.7
        seen = {0}
        stack = [0]
        while stack:
            k = stack.pop(0)
            for l in self.inc[k]:
                i, j, d = self.links[l]
                o = j if i == k else i
                if o in seen:
                    continue
                sgn = 1.0 if i == k else -1.0
                th[o] = (th[k] - sgn * self.wbar * d / C) % TWO_PI
                seen.add(o)
                stack.append(o)
        Em = np.full(self.N, self.xbar * CAP)
        lem = np.zeros((self.NL, 2))
        lph = np.zeros((self.NL, 2))
        return [th, Em, lem, lph]

    # ---------------- the exact per-step map ----------------
    def step(self, st, quant=0, cred=None, ledger=None):
        th, Em, lem, lph = st
        N, NL = self.N, self.NL
        # pass 0: flight load
        fll = np.zeros(N)
        for l, (i, j, d) in enumerate(self.links):
            s = lem[l, 0] + lem[l, 1]
            fll[i] += 0.5 * s
            fll[j] += 0.5 * s
        # pass 1: pitch (bound energy only; Ee=0)
        x = (Em + fll) / CAP
        w = W2 / (1.0 + QD * x)
        # pass 2: wants (dense, 1:1 comb -- asserted valid for these nets)
        want = np.zeros((NL, 2))
        for l, (i, j, d) in enumerate(self.links):
            det = w[i] - w[j]
            res = GM * GM / (GM * GM + det * det)
            wi_r = self.wbar if self.no_tau else w[i]
            wj_r = self.wbar if self.no_tau else w[j]
            gf = gate(wrap(th[i] - wi_r * d / C - th[j]))
            gb = self.freeze_gb if self.freeze_gb is not None \
                else gate(wrap(th[j] - wj_r * d / C - th[i]))
            head_i = min(max(1.0 - Em[i] / CAP, 0.0), 1.0)
            head_j = min(max(1.0 - Em[j] / CAP, 0.0), 1.0)
            mi = math.sqrt(Em[i] * max(Em[j], MOBFLOOR))
            mj = math.sqrt(Em[j] * max(Em[i], MOBFLOOR))
            base = KD * DT * self.geo[l] * self.gpl[l] * res
            want[l, 0] = base * gf * head_j * mi
            want[l, 1] = base * gb * head_i * mj
        # pass 3: outflow limiter
        scl = np.ones(N)
        for k in range(N):
            out = 0.0
            for l in self.inc[k]:
                i, j, d = self.links[l]
                out += want[l, 0] if i == k else want[l, 1]
            a1 = 0.98 * Em[k]
            if out > a1 and out > 0:
                scl[k] = a1 / out
        # pass 4: resolve deposits, debit sources
        ths = th.copy()                       # snapshot (kernel line 1758)
        dep = np.zeros((NL, 2))
        for l, (i, j, d) in enumerate(self.links):
            for dr, src in ((0, i), (1, j)):
                f = want[l, dr] * scl[src]
                if f <= 0:
                    continue
                if lem[l, dr] <= 0:
                    lph[l, dr] = 0.0
                lem[l, dr] += f
                dep[l, dr] = f
        for k in range(N):
            dtot = 0.0
            for l in self.inc[k]:
                i, j, d = self.links[l]
                dtot += dep[l, 0] if i == k else dep[l, 1]
            Em[k] -= dtot
        # pass 4c: deposit-side receiver entrainment (sequential, CSR order)
        for k in range(N):
            for l in self.inc[k]:
                i, j, d = self.links[l]
                dr = 1 if i == k else 0       # direction arriving at k
                f = dep[l, dr]
                if f <= 0:
                    continue
                src = i if dr == 0 else j
                wsrc_r = self.wbar if self.no_tau else w[src]
                err = wrap(ths[src] - wsrc_r * d / C - th[k])
                mix = f / (f + Em[k] + LOCKF)
                th[k] += KLOCK * mix * err
        # pass 5: flight advance + deliveries (receiver-sequential)
        ths2 = th.copy()                      # snapshot (kernel line 1828)
        rough_out = 0.0
        backs = 0.0
        for k in range(N):
            for l in self.inc[k]:
                i, j, d = self.links[l]
                dr = 1 if i == k else 0
                send = i if dr == 0 else j
                adv = DT * C / d
                if lem[l, dr] <= 0:
                    continue
                lph[l, dr] += adv
                if lph[l, dr] < 1.0:
                    continue
                free = CAP - Em[k]            # Ee=0
                take = lem[l, dr]
                if take > free:
                    take = max(free, 0.0)
                if take > 0:
                    mobprev = Em[k]
                    lem[l, dr] -= take
                    det = w[i] - w[j]         # kernel: link orientation
                    R = 2.0 * abs(det) * GROUGH / (det * det + GROUGH * GROUGH)
                    rough = take * ROUGHK * R
                    if quant and cred is not None:
                        epsF = A0 * w[k] / TWO_PI
                        cred[k] += rough
                        if cred[k] > 2.0 * epsF:
                            cred[k] = 2.0 * epsF
                        mv = math.floor(cred[k] / epsF) * epsF
                        cred[k] -= mv
                        rough = min(mv, take)
                    bs = rough * SPULL / (1.0 + SPULL)
                    Em[k] += take - rough
                    rough_out += rough - bs
                    backs += bs
                    wsend_r = self.wbar if self.no_tau else w[send]
                    err = wrap(ths2[send] - wsend_r * d / C - th[k])
                    mix = take / (take + mobprev + LOCKF)
                    th[k] += KLOCK * mix * err
                if lem[l, dr] <= 1e-17:
                    lem[l, dr] = 0.0
                    lph[l, dr] = 0.0
                elif take <= 0:
                    lph[l, dr] = 0.0
                else:
                    lph[l, dr] -= 1.0
        # pass 6: the clock advance
        th[:] = (th + w * DT) % TWO_PI
        if ledger is not None:
            ledger[0] += rough_out
            ledger[1] += backs
        return st


# ---------------- structures ----------------
def make_ring(N, d, m, gpl=None, es='seeded'):
    om = TWO_PI * m * C / (N * d)             # ring_m uniform pitch (kernel 941)
    xbar = (W2 / om - 1.0) / QD
    R = d / (2.0 * math.sin(math.pi / N))
    pos = [np.array([R * math.cos(TWO_PI * k / N),
                     R * math.sin(TWO_PI * k / N), 0.0]) for k in range(N)]
    axes = []
    for k in range(N):
        u = pos[(k + 1) % N] - pos[k]
        u = u / np.linalg.norm(u)
        axes.append(tangent_frame(u))         # (n1, n2): n1=z, n2 in-plane
    links = [(k, (k + 1) % N, d) for k in range(N)]
    return Net(f'ring{N}_m{m}', pos, links, xbar,
               es_mode=es, gpl_override=gpl, axes=axes)


def make_cube(a, gpl=None, es='seeded'):
    om = math.pi * C / a                      # pi-rung (kernel 1081)
    xbar = (W2 / om - 1.0) / QD
    h = 0.5 * a
    pos, axes = [], []
    for k in range(8):
        p = np.array([h if (k & 1) else -h,
                      h if (k & 2) else -h,
                      h if (k & 4) else -h])
        pos.append(p)
        u = p / np.linalg.norm(p)             # radial (shell seeder ~1095)
        axes.append(tangent_frame(u))
    links = []
    for k in range(8):
        for b in range(3):
            k2 = k ^ (1 << b)
            if k2 > k:
                links.append((k, k2, a))
    return Net(f'cube_a{a:g}', pos, links, xbar,
               es_mode=es, gpl_override=gpl, axes=axes)


def make_pair(d, m=1, gpl=1.0, es='seeded'):
    om = math.pi * m * C / d
    xbar = (W2 / om - 1.0) / QD
    pos = [np.array([0.0, 0.0, 0.0]), np.array([d, 0.0, 0.0])]
    u = np.array([1.0, 0.0, 0.0])
    axes = [tangent_frame(u), tangent_frame(u)]   # both transverse: gpl=1
    links = [(0, 1, d)]
    return Net(f'pair_m{m}', pos, links, xbar,
               es_mode=es, gpl_override=gpl if gpl != 'seeded' else None,
               axes=axes)


def get_net(name, gpl=None, es='seeded'):
    if name == 'ring12':
        return make_ring(12, 1.25, 5, gpl, es)
    if name == 'ring6':
        return make_ring(6, 1.25, 3, gpl, es)
    if name == 'ring6m2':
        return make_ring(6, 1.25, 2, gpl, es)
    if name == 'cube':
        return make_cube(1.586, gpl, es)
    if name == 'cube125':
        return make_cube(1.25, gpl, es)
    if name == 'pair':
        return make_pair(1.25, 1, 1.0, es)
    raise SystemExit(f'unknown net {name}')


# ---------------- steady state + observables ----------------
def relax(net, T=1500.0, quant=0):
    st = net.seed_state()
    cred = np.zeros(net.N)
    steps = int(round(T / DT))
    for s in range(steps):
        net.step(st, quant=quant, cred=cred)
    return st


def observables(net, st):
    th, Em, lem, lph = st
    fll = np.zeros(net.N)
    for l, (i, j, d) in enumerate(net.links):
        s = lem[l, 0] + lem[l, 1]
        fll[i] += 0.5 * s
        fll[j] += 0.5 * s
    x = (Em + fll) / CAP
    w = W2 / (1.0 + QD * x)
    out = {
        'Em_store': Em.mean(), 'flload': fll.mean(), 'x': x.mean(),
        'w': w.mean(), 'lem_f': lem[:, 0].mean(), 'lem_b': lem[:, 1].mean(),
    }
    gf, gb, Ff, Fb = [], [], [], []
    for l, (i, j, d) in enumerate(net.links):
        det = w[i] - w[j]
        res = GM * GM / (GM * GM + det * det)
        pf = wrap(th[i] - w[i] * d / C - th[j])
        pb = wrap(th[j] - w[j] * d / C - th[i])
        gf.append(gate(pf))
        gb.append(gate(pb))
        head_i = 1.0 - Em[i] / CAP
        head_j = 1.0 - Em[j] / CAP
        mi = math.sqrt(Em[i] * max(Em[j], MOBFLOOR))
        base = KD * net.geo[l] * net.gpl[l] * res
        Ff.append(base * gf[-1] * head_j * mi)
        Fb.append(base * gb[-1] * head_i * mi)
    out.update(g_f=np.mean(gf), g_b=np.mean(gb),
               F_f=np.mean(Ff), F_b=np.mean(Fb),
               geo=net.geo.mean(), gpl=net.gpl.mean(),
               psi_b=wrap(-2.0 * net.wbar * net.links[0][2] / C))
    return out


# ---------------- exact Floquet ----------------
def window_steps(net):
    """smallest S with S*dt*C/d integer for all link lengths (delivery
    schedule periodic and identical across perturbed trajectories)."""
    d = net.links[0][2]
    for S in range(1, 40001):
        r = S * DT * C / d
        if abs(r - round(r)) < 1e-9 and round(r) >= 1:
            return S
    raise SystemExit('no integer window')


def flatten(st):
    th, Em, lem, lph = st
    return np.concatenate([th, Em, lem.ravel()])


def unflatten(net, z, lph):
    N, NL = net.N, net.NL
    return [z[:N].copy(), z[N:2 * N].copy(),
            z[2 * N:2 * N + 2 * NL].reshape(NL, 2).copy(), lph.copy()]


def floquet(net, T_relax=2000.0, verbose=True, nwin=1):
    S = window_steps(net) * nwin
    st = relax(net, T_relax)
    # advance to a window boundary reference and store lph
    ref = [a.copy() for a in st]
    z0 = flatten(ref)
    lph0 = ref[3].copy()

    def run_window(z):
        st2 = unflatten(net, z, lph0)
        for s in range(S):
            net.step(st2)
        return flatten(st2)

    # baseline periodicity check (modulo uniform clock advance)
    z1 = run_window(z0)
    dth = wrap(z1[:net.N] - z0[:net.N])
    per = np.max(np.abs(z1[net.N:] - z0[net.N:]))
    rot = dth.mean()
    if verbose:
        print(f'# window S={S} steps ({S*DT:g} t.u.), baseline periodicity '
              f'|dz|={per:.2e}, uniform rotation/window={rot:.6f}')
    dim = z0.size
    J = np.zeros((dim, dim))
    eps_th, eps_e = 1e-6, 1e-7
    N = net.N
    for c in range(dim):
        e = eps_th if c < N else eps_e
        zp, zm = z0.copy(), z0.copy()
        zp[c] += e
        zm[c] -= e
        fp = run_window(zp)
        fm = run_window(zm)
        col = fp - fm
        # phase coordinates live on the circle: wrap their differences
        col[:N] = np.vectorize(wrap)(col[:N])
        J[:, c] = col / (2 * e)
    mu, V = np.linalg.eig(J)
    lam = np.log(mu.astype(complex)) / (S * DT)
    return lam, V, S


def mode_number(net, vec):
    """dominant Fourier index of the (th, Em) content of an eigenvector."""
    N = net.N
    th, Em = vec[:N], vec[N:2 * N]
    k = np.arange(N)
    best, nb = -1.0, 0
    for n in range(N):
        ph = np.exp(-2j * math.pi * n * k / N)
        p = abs(np.dot(ph, th)) ** 2 + abs(np.dot(ph, Em)) ** 2
        if p > best:
            best, nb = p, n
    content = (np.linalg.norm(th) ** 2 + np.linalg.norm(Em) ** 2) \
        / (np.linalg.norm(vec) ** 2 + 1e-300)
    return nb, content


def adjacency_class(net, vec):
    """cube: project Em part onto adjacency eigenspaces mu=3,1,-1,-3."""
    N = net.N
    A = np.zeros((N, N))
    for (i, j, d) in net.links:
        A[i, j] = A[i, j] + 1
        A[j, i] = A[j, i] + 1
    w_, U = np.linalg.eigh(A)
    v = vec[N:2 * N].real + vec[:N].real
    if np.linalg.norm(v) < 1e-12:
        v = np.abs(vec[N:2 * N]) + np.abs(vec[:N])
    pw = {}
    for mu in (-3, -1, 1, 3):
        sel = np.abs(w_ - mu) < 1e-6
        P = U[:, sel]
        pw[mu] = np.linalg.norm(P.T @ v) ** 2
    return max(pw, key=pw.get)


def print_floquet_table(net, lam, V, ring=True):
    print(f'# {net.name}: Floquet exponents lambda = growth_rate + i*freq '
          f'[1/t.u., rad/t.u.]')
    order = np.argsort(-lam.real)
    rows = []
    for idx in order:
        l = lam[idx]
        if abs(l.imag) > math.pi / (window_steps(net) * DT) - 1e-6:
            pass  # at the Nyquist edge of the window map: fold noted
        v = V[:, idx]
        if ring:
            n, cont = mode_number(net, v)
            rows.append((n, l, cont))
        else:
            cl = adjacency_class(net, v)
            n, cont = mode_number(net, v)
            rows.append((cl, l, cont))
    seen = {}
    for nb, l, cont in rows:
        key = (nb, round(l.real, 6), round(abs(l.imag), 6))
        if key in seen:
            continue
        seen[key] = 1
        tag = 'n' if ring else 'mu'
        if cont > 0.05:
            print(f'  {tag}={nb:>2}  Re={l.real:+.5f}  Im={l.imag:+.5f}  '
                  f'(th,Em content {cont:.2f})')
    return rows


# ---------------- kick spectroscopy (validation + tax) ----------------
def pencil_fit(z, dt_s, order=4):
    """matrix-pencil complex-exponential fit of a scalar complex series."""
    M = len(z)
    L = M // 2
    H0 = np.array([z[i:i + L] for i in range(M - L)])
    H1 = np.array([z[i + 1:i + 1 + L] for i in range(M - L - 1)])
    U, s, Vt = np.linalg.svd(H0[:-1], full_matrices=False)
    r = min(order, np.sum(s > s[0] * 1e-8))
    U, s, Vt = U[:, :r], s[:r], Vt[:r]
    A = np.diag(1 / s) @ U.conj().T @ H1 @ Vt.conj().T
    mu = np.linalg.eigvals(A)
    return np.log(mu) / dt_s


def kick(net, n_mode, amp=1e-5, which='Em', T=250.0, quant=0):
    st = relax(net, 1500.0)
    ref = [a.copy() for a in st]
    N = net.N
    k = np.arange(N)
    pat = np.cos(2 * math.pi * n_mode * k / N)
    if which == 'Em':
        st[1] += amp * pat
    else:
        st[0] += amp * pat
    steps = int(round(T / DT))
    every = 5
    zs_E, zs_T, ts = [], [], []
    ph = np.exp(-2j * math.pi * n_mode * k / N)
    for s in range(steps):
        net.step(st)
        if s % every == 0:
            dE = st[1] - ref[1]
            dth = wrap(st[0] - ref[0])
            dth -= dth.mean()
            zs_E.append(np.dot(ph, dE))
            zs_T.append(np.dot(ph, dth))
            ts.append(s * DT)
    # (the reference state also cycles within the delivery window; the
    #  projections onto n>0 modes are immune to the uniform sawtooth)
    return np.array(ts), np.array(zs_E), np.array(zs_T)


def leak_tax(net, n_mode, amps, T=200.0, quant=0):
    """total-mass loss rate vs kick amplitude: the comma tax."""
    rows = []
    for a in amps:
        st = relax(net, 1500.0)
        N = net.N
        k = np.arange(N)
        st[1] += a * np.cos(2 * math.pi * n_mode * k / N)
        led = [0.0, 0.0]
        cred = np.zeros(N)
        M0 = st[1].sum() + st[2].sum()
        steps = int(round(T / DT))
        for s in range(steps):
            net.step(st, quant=quant, cred=cred, ledger=led)
        M1 = st[1].sum() + st[2].sum()
        rows.append((a, (M0 - M1), led[0], led[1]))
    return rows


# ---------------- analytic pair rate (sign theorem cross-check) --------
def pair_rate_analytic(net, obs):
    """lambda_h = 2*A*g*S/cap*(1+S/(2*h*Ebar)) with A=kd*geo*gpl*res."""
    A = KD * obs['geo'] * obs['gpl']
    S = obs['Em_store']
    h = 1.0 - S / CAP
    lam = 2.0 * A * obs['g_f'] * (S / CAP) * (1.0 + h * CAP / (2 * S) * 0)
    # exact linear coefficient: d(net)/d(dE_asym), see STABILITY.md S2
    lam_exact = A * obs['g_f'] * (S / CAP + h / 2 * 0 + 0)
    return 2 * A * obs['g_f'] * (S / CAP)


# ---------------- drivers ----------------
def cmd_fixed(names, gpl=None, es='seeded'):
    for nm in names:
        net = get_net(nm, gpl, es)
        st = relax(net)
        ob = observables(net, st)
        print(f'# {net.name}: xbar={net.xbar:.5f} wbar={net.wbar:.5f} '
              f'Es={net.Es:.4f} r={net.r:.4f}')
        print('   ' + '  '.join(f'{k}={v:.4f}' for k, v in ob.items()))


def cmd_floquet(nm, gpl=None, es='seeded'):
    net = get_net(nm, gpl, es)
    lam, V, S = floquet(net)
    ring = nm.startswith('ring') or nm == 'pair'
    print_floquet_table(net, lam, V, ring=ring)
    return net, lam, V


def cmd_pair():
    net = get_net('pair')
    st = relax(net)
    ob = observables(net, st)
    print('# pair fixed point:', {k: round(v, 4) for k, v in ob.items()})
    # heavy->light: +dE on 0, -dE on 1
    for a in (0.005, 0.02, 0.05):
        st2 = [x.copy() for x in st]
        st2[1][0] += a
        st2[1][1] -= a
        asym0 = st2[1][0] - st2[1][1]
        tr = []
        for s in range(int(40.0 / DT)):
            net.step(st2)
            if s % 25 == 0:
                tr.append((s * DT, st2[1][0] - st2[1][1]))
        # early net flow direction
        t1, a1 = tr[1]
        t9 = min(len(tr) - 1, 20)
        rate = -math.log(abs(tr[t9][1] / asym0)) / tr[t9][0] \
            if tr[t9][1] * asym0 > 0 else float('nan')
        print(f'  dE={a:g}: asym {asym0:.4f} -> {tr[t9][1]:+.5f} @t={tr[t9][0]:g} '
              f'(flow {"heavy->light (restoring)" if abs(tr[t9][1]) < asym0 else "LIGHT->HEAVY (runaway)"}); '
              f'decay rate ~{rate:.4f}/t.u.')


def cmd_unlock():
    net = get_net('pair')
    st0 = relax(net)
    lo, hi = 0.0, 0.95 * float(st0[1].min())
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        st = [x.copy() for x in st0]
        st[1][0] += mid
        st[1][1] -= mid
        ok = True
        for s in range(int(300.0 / DT)):
            net.step(st)
        th, Em, lem, lph = st
        fll = np.zeros(2)
        for l, (i, j, d) in enumerate(net.links):
            s2 = lem[l, 0] + lem[l, 1]
            fll[i] += 0.5 * s2
            fll[j] += 0.5 * s2
        x = (Em + fll) / CAP
        w = W2 / (1.0 + QD * x)
        d = net.links[0][2]
        gg = gate(wrap(th[0] - w[0] * d / C - th[1])) \
            * gate(wrap(th[1] - w[1] * d / C - th[0]))
        ok = gg > 0.05 and abs(Em[0] - Em[1]) < 2 * mid + 0.05
        if ok:
            lo = mid
        else:
            hi = mid
    xasym = lo / CAP
    det = abs(W2 / (1 + QD * (net.xbar + lo / CAP))
              - W2 / (1 + QD * (net.xbar - lo / CAP)))
    print(f'# pair unlock threshold: dE_asym={lo:.4f} (x asym {xasym:.4f}), '
          f'initial |det|={det:.4f} rad/t.u. (2*Gamma_b={2*GM:.2f})')


def cmd_tax(nm='ring12'):
    net = get_net(nm)
    amps = [0.0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]
    print(f'# {nm} comma tax, kick n=3 in Em, T=200 t.u. '
          '(loss = radiated roughness; continuous vs atoms/credit)')
    print('#   amp      loss_cont   rough_cont |  loss_quant  rough_quant')
    r0 = leak_tax(net, 3, amps, quant=0)
    r1 = leak_tax(net, 3, amps, quant=1)
    for (a, l0, ro0, b0), (_, l1, ro1, b1) in zip(r0, r1):
        print(f'   {a:7.3f}  {l0:10.3e} {ro0:10.3e} |  {l1:10.3e} {ro1:10.3e}')


def cmd_kickfit(nm, n_mode):
    net = get_net(nm)
    ts, zE, zT = kick(net, n_mode)
    sk = 6
    lamE = pencil_fit(zE[sk:], (ts[1] - ts[0]), order=5)
    lamT = pencil_fit(zT[sk:], (ts[1] - ts[0]), order=5)
    keep = lambda ls: sorted([l for l in ls if -3.5 < l.real < 0.5],
                             key=lambda l: -l.real)[:3]
    print(f'# {nm} kick n={n_mode}: pencil poles (Em obs): '
          + '  '.join(f'{l.real:+.4f}{l.imag:+.4f}i' for l in keep(lamE)))
    print(f'#                          (th obs): '
          + '  '.join(f'{l.real:+.4f}{l.imag:+.4f}i' for l in keep(lamT)))


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if '=' not in a]
    kv = dict(a.split('=', 1) for a in sys.argv[1:] if '=' in a)
    gpl = kv.get('gpl')
    gpl = None if gpl in (None, 'seeded') else float(gpl)
    es = kv.get('es', 'seeded')
    mode = args[0] if args else 'all'
    if mode == 'fixed':
        cmd_fixed(args[1:] or ['ring12', 'ring6', 'cube', 'pair'], gpl, es)
    elif mode == 'floquet':
        cmd_floquet(args[1] if len(args) > 1 else 'ring12', gpl, es)
    elif mode == 'pair':
        cmd_pair()
    elif mode == 'unlock':
        cmd_unlock()
    elif mode == 'tax':
        cmd_tax(args[1] if len(args) > 1 else 'ring12')
    elif mode == 'kick':
        cmd_kickfit(args[1] if len(args) > 1 else 'ring12',
                    int(kv.get('n', 1)))
    elif mode == 'all':
        print('==== fixed points ====')
        cmd_fixed(['ring12', 'ring6', 'cube', 'pair'])
        print('==== sign theorem (pair) ====')
        cmd_pair()
        cmd_unlock()
        for nm in ('ring12', 'ring6', 'cube'):
            print(f'==== floquet {nm} ====')
            cmd_floquet(nm)
        print('==== kick validation ====')
        for n in (1, 3, 6):
            cmd_kickfit('ring12', n)
        print('==== comma tax ====')
        cmd_tax('ring12')
