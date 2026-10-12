# v103 FRAME CARD — read this before any design, and after every compaction

A SessionStart hook (`.claude/settings.json`) prints this card into
context at startup, resume, clear and compaction. It exists because the
user said, verbatim:

> "Be sure that each compaction in this experiment keeps these
> "something something" v "something nothing" retains, because your
> training fails us here, so we have to activly do work to ensure key
> concepts are retained."

Standard physics training pulls back toward the ether frame. This card
pulls the other way.

## THE FRAME (R14): something with gaps, NOT something in something

The user, verbatim (2026-10-11):

> "The trick is in this frame the medium isn't somewhere, it isn't en
> ether. In fact, it is the opposite of the ether. I think we can
> construct around this. ... You keep coming back to the wrong frame on
> the ether. Ether had particles rolling around in an ether "something
> in something" we propose "something with gaps" which is our key out."

| ether frame (WRONG) | gap frame (OURS) |
|---|---|
| a container (space), a substance in it (ether), particles in that | only the medium. Nothing holds it |
| a particle is something in something | a particle is a gap: medium that is not there, "something / nothing" |
| the medium is somewhere | the medium is what "where" is made of. Distance, volume and neighbourhood are counted through the medium |
| a void is an empty slot in a space | a gap is missing where. There is no slot to be empty |
| drawn medium goes somewhere else | there is no elsewhere. The medium closes over its gaps, and there is less space between gaps |
| a body pushes through the medium (drag, ether wind) | a gap moves because the medium re-forms. Nothing pushes through anything |

## Never do these

1. **Never ask an ether-frame question or raise one as an objection:**
   - where does drawn or consumed medium go?
   - what about drag?
   - how does a body move relative to the medium?
   - what is the medium's velocity through space?
2. **Never build a container into a model:**
   - a fixed lattice of slots or cells whose occupancy is "the
     medium";
   - a density per coordinate cell treated as physics;
   - coordinates treated as physical positions.

   Computational indices may only be bookkeeping of quanta. They are
   never places.
3. **Never describe a gap as an object in the medium.** It is an
   absence of medium.

## Frame check: every PREREG must answer these before any run

1. **What are the places?** They must be the medium's own quanta or
   relations.
2. **Can a place exist with no medium in it?** The answer must be no.
3. **How is distance counted?** It must be counted through the medium:
   bonds, or quanta along a path.
4. **What is a gap, and how is it read?** It must be an absence that
   the medium's own state shows. No label (R13).
5. **Which earlier results came from a container?** State their caveat.
   c1, a1 and b1 used fixed lattices of slots or cells.

## The other standing rules (user)

- **R13 no label:** "Any encoding by label is an instant fail
  intrinsically. If the setup of quanta cannot survive without label,
  then it is dead."
- **Light is not matter:** "The underlying field will compose light
  and matter, but light is not matter." Score conversion, never
  occupation.
- **No edge:** "the simulation needs to simulate no edge". Systems are
  closed or periodic, and the closure is the medium's own.
- **Every decision is made from all eight combinations:** none, a, b,
  c, ab, ac, bc, abc.
- **Code is C or Go.** Python only as a thin driver.
- **Commits carry no Co-Authored-By footer.**
- **Do not consult any version before v89.**
- **Do not map onto measured constants.**
- **Do not seed a profile.**
- **Do not scan a constant to place an outcome.**

## Where the work stands

See `README.md` for the status table and the step DECIDE files.

- **Decided (1-D or container, caveated):**
  - c1 E: Coulomb law, signed by sublattice;
  - b1 C;
  - c2 U by the bars, read as like-with-like;
  - a1 P(2).
  - a1, b1 and c2 ran on 1-D rings, and c1, a1 and b1 on containers.
    Spin, exchange, sharpness, binding and 1/r need three dimensions.
- **The user's reading:** two things are going on. One is a signed
  charge. The other is an unsigned gap–gap interaction through 1/r².
- **Now: g1, the 3-D geometry gate, `gate/g1/`. REGISTERED, NOT RUN.**
  The user said "Do not yet run". Read `gate/g1/STATUS.md` for the run
  state, and do not start it without the user's go.
  - The medium is bare quanta and contacts only.
  - A gap is a missing quantum that the medium closes over, read from
    the medium's own shape.
