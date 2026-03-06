import itertools
import numpy as np
import os

print("HFKT v13: Algebraic Exploration of 3-Fold Knot States")
print("=====================================================")

os.makedirs("v13/results", exist_ok=True)

# In Cl+(3,0,1), the bulk quaternion (rotor) has 4 components: 1 scalar, 3 bivectors.
# A topological knot involves a mapping from S^3 (physical space) to S^3 (rotor space).
# The fundamental traits of a fractional knot (quark) in this algebra:
# 1. Winding fraction (Charge / B-number contribution)
# 2. Orientation (Iso-spin / Flavor projection)
# 3. Phase/Color (Permutation symmetry constraint for multi-knot stability)

# Let's define the fundamental fractional states (quarks).
# We need 2 distinct flavors (Up, Down) to make protons (uud) and neutrons (udd).

class FractionalKnot:
    def __init__(self, name, b_frac, isospin, phase_idx):
        self.name = name
        self.b_frac = b_frac       # Contribution to total Baryon number B=1
        self.isospin = isospin     # +1/2 for Up, -1/2 for Down
        self.phase_idx = phase_idx # 0, 1, 2 (analogous to Red, Green, Blue)

    def __repr__(self):
        return f"{self.name}(B={self.b_frac:.2f}, I={self.isospin}, P={self.phase_idx})"

# Define the base states
def get_base_states():
    states = []
    for phase in [0, 1, 2]:
        # Up quark: B=1/3, Isospin = +1/2
        states.append(FractionalKnot('U', 1/3, 0.5, phase))
        # Down quark: B=1/3, Isospin = -1/2
        states.append(FractionalKnot('D', 1/3, -0.5, phase))
    return states

# Combinatorics of 3-fold states
def explore_3fold_states():
    base_states = get_base_states()
    
    # We want to find stable combinations of 3 knots.
    # Rules for a stable macroscopic knot (Baryon):
    # 1. Total B = 1 (sum of b_frac = 1)
    # 2. Must be totally asymmetric in phase to cancel internal boundaries (Color singlet: P=0,1,2 must all be present)
    
    valid_combinations = []
    
    out_tsv = "v13/results/algebraic_combinations_raw.tsv"
    print(f"\nStreaming intermediate algebraic combinations to {out_tsv}...")
    
    with open(out_tsv, "w") as f:
        f.write("Knot1_Phase\tKnot2_Phase\tKnot3_Phase\tTotal_B\tTotal_Iso\tUnique_Phases\tStable\n")
        
        for combo in itertools.combinations(base_states, 3):
            b_total = sum(k.b_frac for k in combo)
            phases = set(k.phase_idx for k in combo)
            iso_total = sum(k.isospin for k in combo)
            
            # Stability criteria
            is_stable = (abs(b_total - 1.0) < 1e-6) and (len(phases) == 3)
            
            k1, k2, k3 = [f"{k.name}{k.phase_idx}" for k in combo]
            f.write(f"{k1}\t{k2}\t{k3}\t{b_total:.3f}\t{iso_total:.1f}\t{len(phases)}\t{is_stable}\n")
            
            if is_stable:
                flavors = sorted([k.name for k in combo])
                flavor_str = "".join(flavors)
                # Keep original tracking
                valid_combinations.append((combo, flavor_str, iso_total))
            
    print("\nValid Stable 3-Fold Knots (B=1, Color Singlet):")
    print("-" * 50)
    
    # Group by flavor string
    flavor_groups = {}
    for combo, f_str, iso in valid_combinations:
        if f_str not in flavor_groups:
            flavor_groups[f_str] = {"iso": iso, "count": 0, "examples": []}
        flavor_groups[f_str]["count"] += 1
        # Store just one example
        if len(flavor_groups[f_str]["examples"]) == 0:
             flavor_groups[f_str]["examples"].append(combo)
            
    for f_str, data in flavor_groups.items():
        print(f"Flavor Configuration: {f_str}")
        print(f"  Total Isospin: {data['iso']}")
        print(f"  Number of unique phase permutations: {data['count']}")
        print(f"  Example: {data['examples'][0]}")
        print("")

    return flavor_groups

# Vibrational Modes Framework
def analyze_vibrational_modes(flavor_groups):
    print("\nVibrational Mode Analysis (No-Friction Field):")
    print("-" * 50)
    print("In a frictionless Cl+(3,0,1) field, energetic perturbations do not decay.")
    print("They resolve into normal modes of the composite topological structure.\n")
    
    for f_str, data in flavor_groups.items():
        print(f"Analyzing {f_str} structure:")
        
        # Determine symmetry of the composite
        # If DDU, 2 are the same, 1 is different. 
        # If UUU, all 3 are the same.
        unique_flavors = set(f_str)
        
        if len(unique_flavors) == 1:
            sym = "Symmetric (Equilateral Tension)"
            self_modes = "3 degenerate self-breathing modes."
            mutual_modes = "2 degenerate mutual shear modes, 1 common rotational mode."
        elif len(unique_flavors) == 2:
            sym = "Asymmetric (Isosceles Tension)"
            majority_char = max(set(f_str), key=f_str.count)
            minority_char = min(set(f_str), key=f_str.count)
            self_modes = f"2 degenerate {majority_char}-breathing modes, 1 distinct {minority_char}-breathing mode."
            mutual_modes = "1 longitudinal (majority vs minority) wiggling mode, 1 transverse shear mode."
        else:
            sym = "Fully Asymmetric (Scalene Tension)"
            self_modes = "3 distinct self-breathing modes."
            mutual_modes = "3 distinct mutual spatial wiggling modes."
            
        print(f"  - Geometry: {sym}")
        print(f"  - Internal Self-Knot Ringing: {self_modes}")
        print(f"  - Mutual Knot Wiggling: {mutual_modes}")
        print("")

if __name__ == "__main__":
    groups = explore_3fold_states()
    analyze_vibrational_modes(groups)
