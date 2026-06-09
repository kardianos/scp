import Lake
open Lake DSL

package «scp_v61» where

/-
v61 Lean modules (nonlinear/curved gravity + EW-vev home).

BUILD (same method as v60): Mathlib v4.29.0 is already compiled in the v59 tree;
check a module against that environment without a fresh Mathlib build:

  cd ../../v59/furey_construction/lean
  lake env lean ../../../v61/lean/CurvedBackreaction.lean

(empty output = clean). Axiom-check by appending `#print axioms <thm>`.
-/

@[default_target]
lean_lib CurvedBackreaction where
  -- GEN1: Schwarzschild vacuum from rho_grav; weak-field = v60 Newtonian; deflection 4GM/b

lean_lib Backreaction where
  -- GEN2: m'(r)=4πr²ρ from G^t_t=-8πρ; M=∫ρ_grav; Newtonian Poisson; Schwarzschild matching

lean_lib EwVevHome where
  -- GEN3: R1 home = End(L) Frobenius Higgs; v=784a²; 783 Goldstones; democracy not selected

lean_lib GravitationalWaves where
  -- GEN4: 2 TT modes = h+,hx (helicity ±2 = LIGO); quadrupole 32/5; GWs at c — LIGO closed

lean_lib Perihelion where
  -- GEN5: Schwarzschild orbit GR term (coeff 3) -> precession 6πGM/(c²a(1-e²)); Mercury 43″/cy

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.29.0"
