import re
import sys

def convert_to_rat(content):
    content = content.replace("abbrev ℝ := Float", "abbrev ℝ := Rat")
    content = content.replace("Float.abs", "Rat.abs")
    content = content.replace("Float", "Rat")
    
    # Simple regex to replace decimals with fractions or integers
    # Be careful not to replace inside strings randomly if not needed, but here it's mostly code
    def replacer(match):
        val_str = match.group(0)
        try:
            # if it's an integer like 100.0 -> 100
            f = float(val_str)
            if f.is_integer():
                return str(int(f))
            else:
                # convert to fraction
                from fractions import Fraction
                frac = Fraction(str(val_str))
                return f"({frac.numerator}/{frac.denominator})"
        except:
            return val_str

    content = re.sub(r'\b\d+\.\d+\b', replacer, content)
    content = content.replace("1e-6", "0")
    
    return content

for filename in ["StabilityFromAlgebra.lean", "PhaseC_Certificates.lean"]:
    with open(f"/home/d/code/scp/v59/furey_construction/lean/7D_Algebra/{filename}", "r") as f:
        content = f.read()
    
    content = convert_to_rat(content)
    
    with open(f"/home/d/code/scp/v59/furey_construction/lean/7D_Algebra/{filename}", "w") as f:
        f.write(content)

print("Refactored to Rat")
