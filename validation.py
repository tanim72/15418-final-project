
def validate_tridiagonal_solution(input_file, solver_output_file, tolerance=1e-7):
    """
    Validate a tridiagonal solver output by comparing it with the seqeuntial solution.
    """
    tokens = []
    with open(input_file, 'r') as f:
        for line in f:
            if line.strip():
                tokens.extend(line.split())
    
    # Matrix size
    N = int(tokens[0])
    idx = 1
    
    # Main diagonal
    b = list(map(float, tokens[idx: idx + N]))
    idx += N
    
    # Super diagonal
    c = list(map(float, tokens[idx: idx + N - 1]))
    idx += (N - 1)
    c.append(0.0)
    
    # Sub diagonal
    a = [0.0]
    a.extend(list(map(float, tokens[idx: idx + N - 1])))
    idx += (N - 1)
    
    # RHS vector d
    d = list(map(float, tokens[idx: idx + N]))
    idx += N

    def thomas_solver(N, a, b, c, d):
        """
        Sequential Thomas algorithm as implemented in solverMatrix.cpp.
        """

        c_star = [0.0] * N
        d_star = [0.0] * N

        # Forward pass: first row
        c_star[0] = c[0] / b[0]
        d_star[0] = d[0] / b[0]

        # Process rows 1 to N-1
        for i in range(1, N):
            denom = b[i] - a[i] * c_star[i - 1]
            if i < N - 1:
                c_star[i] = c[i] / denom
            d_star[i] = (d[i] - a[i] * d_star[i - 1]) / denom

        # Backward substitution
        x = [0.0] * N
        x[N - 1] = d_star[N - 1]
        for i in range(N - 2, -1, -1):
            x[i] = d_star[i] - c_star[i] * x[i + 1]
        return x

    x_reference = thomas_solver(N, a, b, c, d)

    with open(solver_output_file, 'r') as f:
        x_solver = [float(line.strip()) for line in f if line.strip()]
    
    if len(x_solver) != N:
        print(f"Error: The output file contains {len(x_solver)} values, but {N} were expected.")
        return False

    max_diff = max(abs(x_reference[i] - x_solver[i]) for i in range(N))
    
    if max_diff < tolerance:
        print(f"PASS = {max_diff:.2e} < {tolerance}")
        return True
    else:
        print(f"FAIL = {max_diff:.2e} >= {tolerance}")
        return False

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python validation.py <input_test_file> <solver_output_file>")
        sys.exit(1)
        
    input_filename = sys.argv[1] # Input test file
    output_filename = sys.argv[2] # Output file from the solver
    
    validate_tridiagonal_solution(input_filename, output_filename)
