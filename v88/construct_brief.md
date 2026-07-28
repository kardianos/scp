Read v88/PRINCIPLE.md. That is the foundation. Then read v88/DESIGN_CLEAN.md for
the partial structure already built on it (noting its header: it is subordinate,
and anything in it that assumes a background is superseded).

=== YOUR ROLE ===
You are a PROPONENT, not a reviewer. Previous seats were asked to attack this;
you are not. Your job is to CONSTRUCT this universe more fully — to make it
work. Take the principle as given and build outward from it: find the
structures that make it coherent, derive what follows, supply what is missing.

Being a proponent does not mean being uncritical. It means that when you hit a
gap you try to FILL it rather than report it. If you find something that cannot
be made to work, say so — but only after genuinely trying to construct a way
through, and say what you tried.

You may write new files under v88/. Do NOT modify sfa/sim/scp_sim.c or
sfa/format/sfa.h.

=== THE ONE CONSTRAINT THAT IS ABSOLUTE ===
NO BACKGROUND. Space is a MODE of energy, not a stage. There is nothing for
anything to move through, sit in, or be located against.

This is the constraint every previous attempt violated, including mine and every
instrument in v88. The failure is subtle and it is always the same shape:
someone writes "the cells are space" and then writes a model with a permanent
index set — coordinates that persist while contents evolve. `psi[i]` on a fixed
lattice is a background. `x(i) = a*i + xi(i)` is a background (a*i is a
permanent reference lattice and xi is a wiggle away from it). Fixed
connectivity is a background. Immortal cell labels are a background.

If your construction has anything that persists and is merely re-valued, you
have reintroduced the stage. Check for this explicitly and say how you avoided
it.

Related: NO IMPORTED FIELDS OR SPECIES. Nothing may be added to carry an effect.
If an effect requires a new field, that is a failure of the construction, not a
licence to add one.

=== WHAT TO CONSTRUCT ===
PRINCIPLE.md section 7 lists what is open. In priority order:

1. THE REGULAR METHOD. This is the heart. There are no equations of motion
   because nothing has a trajectory; what replaces them is a rule for which
   successions of structure are permitted. Construct it. What is the state?
   What is a permitted transition? How does conservation of energy across modes
   constrain it? A structure-to-structure rule, not a differential equation on
   an index set.

2. DIRECTION AND RATE OF CONVERSION. Conservation says energy changes mode; it
   does not say which way a transformation runs or when. Construct the
   condition. The natural candidate from the design is that conversion proceeds
   where chiral cycles can close and stalls where they cannot — develop that, or
   find something better.

3. THE MODE STRUCTURE. Space, dense pattern (mass), field pattern (EMF) are
   established. What distinguishes the patterning of mass from that of EMF?
   Is the list complete? What forces the distinction?

4. CURVATURE, QUANTITATIVELY. Section 4.3 derives that space curves because
   matter is converted space and energy is never destroyed. Make it
   quantitative. What is the relation between converted energy and the
   distortion of the remaining space? Does anything recognisable follow?

5. WHAT A PARTICLE IS, under this principle. Not a trapped lump (that is the
   background picture). A self-sustaining pattern of conversion. Construct it.
   What makes such a pattern persist? What makes it discrete? What sets mass?

=== USEFUL MEASURED FACTS ===
* Link-mediated transfer gives an emergent speed: c ~ 1.05*g*a, constant over a
  4x range in coupling. Transfer through a mediating harmonic is rate-limited by
  cycle completion, and a detuned link suppresses delivery 227x. [M]
* On a discrete fabric, per-component charge conserved to 13 significant figures
  while momentum conserved only to ~1e-4. Read as retrodiction: charge is a
  property of structure and does not reference position; momentum references a
  background and there is none. [M]
* Chiral link actuation in the two transverse planes gives a signed handedness
  chi = 2 Im(u1* u2), well-defined only where a transverse plane exists. [M]

=== HOW TO WORK ===
Construct first, then check. Distinguish what you DERIVE from the principle,
what you POSTULATE to make it work, and what you GUESS. A proponent may
postulate freely — but must label it, so the structure can be audited later.

Mathematics is welcome and encouraged. If a construction suggests a computation
that does not need a background, do it. If you cannot see how to compute
something without reintroducing a stage, say that explicitly — it is important
information.

Note: Quantum Ring Theory (Guglinski) has been mentioned as conceptually
adjacent. Treat as a pointer only; do not import its framework.

Write to v88/CONSTRUCTION.md. List any code you create.
