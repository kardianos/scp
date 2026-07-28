Resumed session. Your v88/GROK_V88_DESIGN.md + fabric_harmonic.c are good work and
the theory section is the right structure. But the instrument does not yet show
what it needs to show, and I want you to attack that with a different method.
DO NOT modify sfa/sim/scp_sim.c or sfa/format/sfa.h. Work under v88/.

=== PUSHBACK 1: your own smoke test FAILS your own P1, and is getting worse ===
Built and ran fabric_harmonic smoke (d=3, L=22, N=10648, M=3):

  t=3.2   n_lumps=469  mean_Nc=2.11  sd/mean=1.488  peakI=0.315
  t=6.4   n_lumps=629  mean_Nc=3.21  sd/mean=2.859  peakI=0.207

P1 requires CV < 0.15 inside a peak. This is 1.49 -> 2.86, i.e. 10-20x outside,
and RISING with time. Every "species" is Nc ~ 1-2 cells, i.e. single cells, not
clusters. theta ~ -0.10 to -0.17, so the densification coupling has barely
engaged and the "tighter cell" mechanism is not doing anything yet.

=== PUSHBACK 2: your lock vectors may be counting your truncation ===
The recurring locks are (2,-1,0) and (-1,2,-1). In a 3-mode tower those are the
LOWEST-ORDER channels available. That is exactly your own tension #2 ("if results
depend sensitively on m, the cutoff is physical"). I am running m=2 and m=4
builds now. If the species set changes with m, the cutoff is doing the work and
the integers are an artifact of truncation, not of the fabric.

=== PUSHBACK 3: single cells are not particles ===
Nc ~ 1 means you are labelling individual cells by their internal mode mix. That
is not an emergent particle -- it is a per-cell state vector. A particle must be
a MULTI-CELL cluster whose interior is locked and whose exterior is detuned.
Nothing in the smoke run shows interior/exterior structure at all.

=== WHAT I WANT NEXT: SYMBOLIC REFINEMENT WITH FITTED BACK-TRACING ===
Stop postulating the Hamiltonian and simulating forward. Do the inverse.

1. RUN the instrument (or a reduced version) and MEASURE the effective dynamics:
   the actual secular drift of I_alpha and phi_alpha, the actual theta response,
   the actual inter-cell transfer as a function of detuning.
2. FIT symbolic expressions to those measured quantities -- amplitude equations,
   averaged/secular forms, effective potentials. The symbolic terms should be
   FITTED to the numerics, not assumed.
3. BACK-TRACE: from the fitted effective terms, work backwards to the MINIMAL
   symbolic structure that must be present in the microscopic model to produce
   them. State which terms are forced, which are redundant, and which of your
   current postulates are NOT needed at all.
4. Use that to REDESIGN the model: keep only what back-tracing shows is
   necessary, and add whatever the fit says is missing (in particular whatever
   is needed to get multi-cell clusters and interior/exterior detuning, which is
   absent now).

Use sympy or maxima for the symbolic work; both are available. Show the fitted
coefficients with their residuals -- a fit you cannot quantify is a guess.

=== SPECIFIC QUESTIONS THE BACK-TRACING SHOULD ANSWER ===
 (a) What is the effective interaction between two ADJACENT locked cells, as a
     fitted symbolic form? Is it attractive? If not, clusters can never form and
     the model is dead as written -- say so.
 (b) What actually sets Nc? Fit it, do not assert it.
 (c) Does the theta (density) sector do any work in the current model, or is it
     decorative? theta ~ -0.1 suggests decorative. If GAM/thmax need to be orders
     larger, derive the scale from the fit rather than guessing.
 (d) Is the commensurate-vs-detuned tower choice (your tension #3) actually
     controlling species, or is it irrelevant? Test symbolically.
 (e) Which of your [P] postulates survive back-tracing, and which are dead
     weight? I want the list shortened, not lengthened.

Deliverable: v88/GROK_V88_SYMBOLIC.md with the fitted forms, residuals, the
back-traced minimal structure, and a revised model spec. Be blunt about what is
not working -- I have run your instrument and it is not working yet.
