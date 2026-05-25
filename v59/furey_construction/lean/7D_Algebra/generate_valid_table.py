fano_triples = [
    (1,2,3), (1,4,5), (1,7,6), (2,4,6), (2,5,7), (3,4,7), (3,6,5)
]

table = [[0]*8 for _ in range(8)]
for i in range(8): table[0][i] = (1, i)
for i in range(8): table[i][0] = (1, i)
for i in range(1, 8): table[i][i] = (-1, 0)

for (i, j, k) in fano_triples:
    table[i][j] = (1, k)
    table[j][i] = (-1, k)
    
    table[j][k] = (1, i)
    table[k][j] = (-1, i)
    
    table[k][i] = (1, j)
    table[i][k] = (-1, j)

# Format for Lean
lines = []
for i in range(8):
    row_strs = []
    for j in range(8):
        sgn, val = table[i][j]
        row_strs.append(f"({sgn},{val})")
    lines.append("  [" + ", ".join(row_strs) + "]")

print("[\n" + ",\n".join(lines) + "\n]")
