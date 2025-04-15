#include <iostream>
#include <vector>
using namespace std;

/**
 * Solve the tridiagonal system A*x = d using the Thomas algorithm.
 *
 * @param N  Size of the system (number of equations).
 * @param a  Sub-diagonal of size N (though a[0] is often unused).
 * @param b  Main diagonal of size N.
 * @param c  Super-diagonal of size N (though c[N-1] is often unused).
 * @param d  Right-hand side (RHS) vector of size N.
 * @return   The solution vector x of size N.
 */
vector<double> thomasSolver(int N,
                            const vector<double>& a,
                            const vector<double>& b,
                            const vector<double>& c,
                            const vector<double>& d)
{
    // gamma holds the temporary coefficients for the superdiagonal
    // rho holds the temporary coefficients for the RHS
    vector<double> gamma(N, 0.0);
    vector<double> rho(N, 0.0);

    // Forward pass
    // First row
    gamma[0] = c[0] / b[0];
    rho[0] = d[0] / b[0];

    // Process rows 1 through (N-1)
    for (int i = 1; i < N; i++)
    {
        double denom = b[i] - a[i] * gamma[i - 1];
        // For gamma, only compute if i < N-1 
        if (i < N - 1)
            gamma[i] = c[i] / denom;
       rho[i] = (d[i] - a[i] * rho[i - 1]) / denom;
    }

    // Backward substitution
    vector<double> x(N, 0.0);
    x[N - 1] = rho[N - 1];

    for (int i = N - 2; i >= 0; i--)
    {
        x[i] = rho[i] - gamma[i] * x[i + 1];
    }

    return x;
}

/*
Reads the size of matrix N, parses the vectors: a, b, c, d
and solves the system using the Thomas algorithm and 
prints the result
*/
int main()
{
    int N;
    cin >> N;

    vector<double> a(N), b(N), c(N), d(N);

    // Read b (main diagonal)
    for (int i = 0; i < N; i++)
        cin >> b[i];

     // Read a (subdiagonal), the subdiagonal starts from row 1
     a[0] = 0.0;
     for (int i = 1; i < N; i++)
         cin >> a[i];
         
    // Read c (superdiagonal)
    for (int i = 0; i < N - 1; i++)
        cin >> c[i];
    // usually c[N-1] is unused/zero, but set so no initialization error
    c[N-1] = 0.0;

   

    // Read d (RHS vector)
    for (int i = 0; i < N; i++)
        cin >> d[i];

    // Solve using the Thomas method
    vector<double> x = thomasSolver(N, a, b, c, d);

    // Print the solution
    cout << "Solution x: ";
    for (int i = 0; i < N; i++)
        cout << x[i] << " ";
    cout << endl;

    return 0;
}
