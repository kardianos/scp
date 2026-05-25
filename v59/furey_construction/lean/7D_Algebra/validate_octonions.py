octMultTable = [
  [(1,0), (1,1), (1,2), (1,3), (1,4), (1,5), (1,6), (1,7)],
  [(1,1), (-1,0), (1,3), (-1,2), (1,5), (-1,4), (1,7), (-1,6)],
  [(1,2), (-1,3), (-1,0), (1,1), (1,6), (1,7), (-1,4), (-1,5)],
  [(1,3), (1,2), (-1,1), (-1,0), (1,7), (-1,6), (1,5), (-1,4)],
  [(1,4), (-1,5), (-1,6), (-1,7), (-1,0), (1,1), (1,2), (1,3)],
  [(1,5), (1,4), (-1,7), (1,6), (-1,1), (-1,0), (-1,3), (1,2)],
  [(1,6), (1,7), (1,4), (-1,5), (-1,2), (1,3), (-1,0), (-1,1)],
  [(1,7), (-1,6), (1,5), (1,4), (-1,3), (-1,2), (1,1), (-1,0)]
]

def mul(a, b):
    sgn_a = 1 if a >= 0 else -1
    sgn_b = 1 if b >= 0 else -1
    abs_a = abs(a)
    abs_b = abs(b)
    
    sgn_res, val = octMultTable[abs_a][abs_b]
    return sgn_a * sgn_b * sgn_res * val

errors = 0
# Check quaternion subalgebras
for i in range(1, 8):
    for j in range(1, 8):
        if i == j: continue
        k = abs(mul(i, j))
        
        # i, j, k should form an associative subalgebra
        # e_i (e_j e_i) == (e_i e_j) e_i
        left = mul(i, mul(j, i))
        right = mul(mul(i, j), i)
        if left != right:
            print(f"Alternative failure: e{i}*(e{j}*e{i}) != (e{i}*e{j})*e{i}")
            errors += 1
            
        # check associativity inside the triple
        if mul(i, mul(j, k)) != mul(mul(i, j), k):
            print(f"Quat Assoc failure: e{i}*(e{j}*e{k}) != (e{i}*e{j})*e{k}")
            errors += 1

print(f"Total errors: {errors}")
