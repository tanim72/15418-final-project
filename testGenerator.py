n = 10 # change "N" as you see fit

# main diagonal
b = [(100 + i) if i % 2 == 0 else (i + 93) for i in range(n)]

# lower diagonal
a = [((i % 7) + 1) if i % 2 == 0 else -((i % 7) + 1) for i in range(n)]

# upper diagonal
cycle = [-2, 3, -4, 5]
c = [cycle[i % 4] for i in range(n)]

x = [1] * n

# rhs computation
d = [0] * n
d[0] = b[0] + c[0]  # since x[0]=1 and x[1]=1
for i in range(1, n-1):
    d[i] = a[i] + b[i] + c[i]
d[n-1] = a[n-1] + b[n-1]

output_filename = f"inputs/test{n}x{n}_input.txt" # change this to your desired output filename
with open(output_filename, "w") as f:
    f.write(str(n) + "\n")
    f.write(" ".join(map(str, b)) + "\n")
    f.write(" ".join(map(str, a)) + "\n")
    f.write(" ".join(map(str, c)) + "\n")
    f.write(" ".join(map(str, d)) + "\n")

print(f"Test case saved to {output_filename}")
