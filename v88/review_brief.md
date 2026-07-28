You are reviewing a physics design from FIRST PRINCIPLES. Read
v88/DESIGN_CLEAN.md. That document is the entire design. Do not look for or
rely on project history, earlier versions, or how it was arrived at — it is
deliberately written without them, and prior context would bias the review.

Work in /home/d/code/scp/v88. You may write new files there. Do NOT modify
sfa/sim/scp_sim.c or sfa/format/sfa.h.

=== THE TWO METHODOLOGICAL CONSTRAINTS, WHICH ARE ABSOLUTE ===

1. FIRST PRINCIPLES ONLY.
   Judge the design on whether it is coherent and whether its claims follow.
   Do not appeal to what is conventional, to what other models do, or to
   authority. If a claim does not follow, say why. If it does follow, say from
   what.

2. TRUE EXPERIMENT ONLY — NO SIMPLIFICATION, NO REDUCED DIMENSIONALITY.
   Any numerical work must be the real thing:
     * FULL 3D. No 1D chains, no 2D reductions, no "for tractability" cuts.
       A reduced dimension removes the transverse plane and therefore removes
       handedness, which is load-bearing here (§2, §3). A 1D or 2D test does
       not probe this design at all and its results are inadmissible.
     * FULL NONLINEAR. No linearised proxies. A linear model can show
       propagation and gating but cannot bind; linear results must not be
       offered as evidence about binding, species, or mass.
     * NO SUBSTITUTE MECHANISMS. Do not replace a design element with a
       simpler one that "behaves similarly". If an element is hard to
       implement, implement it or state plainly that you did not.
   If you cannot run the true experiment within your budget, say so and report
   what you did instead, labelled as inadmissible. A simplified result reported
   as if it tested the design is the failure mode I most want to avoid.

=== WHAT I WANT FROM YOU ===

A. COHERENCE. Is the design internally consistent? Specifically:
   - Does "energy crosses only in complete cycles" actually follow from
     link-mediated transfer, or does it need something further?
   - Does c emerge as claimed, or is something smuggled in?
   - Is "mass is a cell count" consistent with objects that have a sharp
     boundary but also an exterior distortion (§5.2)? Where does the mass of
     the distortion live?
   - Is the 3-from-7 claim (§3) load-bearing, decorative, or unsupported?
     The cross product existing only in 3 and 7 dimensions is a fact; whether
     it implies this fabric has 7 degrees of freedom and 3 emergent dimensions
     is what I want examined.

B. DOES IT PRODUCE PARTICLES? §7 lists what the design must produce. For each:
   is it plausible, is it derivable, and what would demonstrate it? Be
   concrete about which are consequences of the design and which are hopes.

C. THE HARD PROBLEM: BINDING. §7 requires an attractive channel that vanishes
   with separation. Handedness (§2) gives two distinct channels. Does an
   attractive one exist? Derive it if you can. If the design as stated cannot
   bind, say so plainly — that is the single most valuable thing you can tell
   me, and it is better said early than after a large simulation.

D. IMPLEMENTATION. If and only if you can do it in full 3D and fully
   nonlinearly, build it and report real measurements. Otherwise give a
   specification precise enough to implement, and say explicitly that you did
   not run it.

E. WHAT WOULD FALSIFY THE DESIGN. Give the cheapest true experiment that could
   kill it. Not a proxy — a real one.

=== HOW TO ANSWER ===

Separate clearly what you DERIVE, what you POSTULATE, and what you GUESS.
Quantify anything you assert numerically, with residuals where you fit.
If part of the design is incoherent or unworkable, say so directly with the
reason — that is more useful to me than agreement. Do not pad.

Write to v88/REVIEW_FIRST_PRINCIPLES.md and list any code you create.
