# V24-180B Results: UUD/UDD in the 180-degree Anti-Phase State

## Parameters

mu=-20, kappa=20, m_U=1.0, A=0.8, sigma=3.0, Nx=4000, xmax=100, tfinal=10000

## Mass Ordering Summary

| m_D  | E_UUD   | E_UDD   | E_UDD - E_UUD | fc_UUD | fc_UDD | Stable? |
|------|---------|---------|---------------|--------|--------|---------|
| 1.00 | 1.2638  | 1.2638  |  0.0000       | 1.000  | 1.000  | Both    |
| 0.95 | 1.2421  | 1.1959  | -0.0461       | 1.000  | 1.000  | Both    |
| 0.90 | 0.4753  | 0.4351  | -0.0402       | 0.195  | 0.017  | Neither |
| 0.85 | 0.3520  | 0.2946  | -0.0574       | 0.236  | 0.198  | Neither |
| 0.80 | 0.2807  | 0.2323  | -0.0484       | 0.212  | 0.274  | Neither |
| 0.70 | 0.2907  | 0.1329  | -0.1578       | 0.174  | 0.196  | Neither |

## Key Findings

### 1. Mass ordering: UUD > UDD at ALL m_D values (proton heavier)

The delta E_UDD - E_UUD is NEGATIVE everywhere m_D < m_U. This means
E_UDD < E_UUD, i.e., the UDD configuration (neutron analog) is LIGHTER
than UUD (proton analog). This is the OPPOSITE of the physical
neutron > proton mass ordering.

At m_D = 1.0, the two are degenerate by symmetry (as expected).

### 2. Stability window is narrow

Only m_D = 1.00 and m_D = 0.95 produce stable oscillons (fc > 0.99).
For m_D <= 0.90, both configurations disperse (fc < 0.25). The 180-degree
state is sensitive to mass asymmetry -- even a 10% mass difference
destabilizes the oscillon.

### 3. The ordering does NOT flip in the 180-degree regime

The key question was whether the UDD > UUD ordering (seen at 120 degrees)
persists in the 180-degree regime where binding is strong. Answer: NO.
The ordering is REVERSED. In the 180-degree state:
- E_UDD < E_UUD (UDD lighter, opposite to nature)
- This holds for both stable (m_D=0.95) and dispersed (m_D<0.90) cases

### 4. Physical interpretation

The reversal makes sense: UDD has TWO down fields with smaller mass m_D,
so its total rest-mass content is inherently lower. The binding energy
(Ep) is similar for both configurations at given m_D, so the mass
ordering is dominated by the free-field mass contribution, not the
interaction.

At m_D = 0.95 (only stable asymmetric point):
- omega_UUD = 0.852 vs omega_UDD = 0.828 (both below mass gap)
- UUD binding: E drops from ~3.37 to 1.24 (63% radiated)
- UDD binding: E drops from ~3.17 to 1.20 (62% radiated)
- The 3.7% energy difference tracks the 5% mass difference

### 5. Comparison with 120-degree results

The V24-UDD result at 120 degrees showed correct neutron > proton ordering.
That result relied on differential coupling strength (the interaction
term breaks the degeneracy differently when P ~ 0). In the 180-degree
state, P is large and the dominant effect is simply the rest-mass
asymmetry, which favors UDD being lighter.

## Conclusion

The UDD > UUD mass ordering does NOT persist in the 180-degree anti-phase
regime. The 180-degree state produces UUD > UDD (proton heavier), opposite
to nature. The binding energy is too symmetric between UUD and UDD to
overcome the bare mass asymmetry from having more light fields in UDD.

## Output Files

- `data/uud180_mD{val}_ts.tsv` — UUD time series per m_D
- `data/udd180_mD{val}_ts.tsv` — UDD time series per m_D
- `data/mass_ordering.tsv` — summary table
- `data/*_spectrum.tsv` — DFT spectra per configuration
