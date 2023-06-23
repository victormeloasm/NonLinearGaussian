#include <stdio.h>
#include <math.h>

#define N 3 // Number of equations in the system
#define MAX_ITER 1000 // Maximum number of iterations
#define EPSILON 0.000001 // Convergence criterion

// System of equations functions
void equations(double x[N], double f[N]) {
    // Enter your equations here
    f[0] = x[0] + 2 * x[1] - x[2] - 2;
    f[1] = x[0] * x[1] + x[1] * x[2] - 1;
    f[2] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 3;
}

// Jacobian matrix
void jacobian(double x[N], double J[N][N]) {
    // Enter the derivative of each equation with respect to each variable here
    J[0][0] = 1.0;
    J[0][1] = 2.0;
    J[0][2] = -1.0;
    J[1][0] = x[1];
    J[1][1] = x[0] + x[2];
    J[1][2] = x[1];
    J[2][0] = 2.0 * x[0];
    J[2][1] = 2.0 * x[1];
    J[2][2] = 2.0 * x[2];
}

// Gaussian elimination method
void gaussianElimination(double A[N][N + 1], double x[N]) {
    int i, j, k;

    for (i = 0; i < N; i++) {
        // Search for maximum in this column
        double maxEl = fabs(A[i][i]);
        int maxRow = i;
        for (k = i + 1; k < N; k++) {
            if (fabs(A[k][i]) > maxEl) {
                maxEl = fabs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row
        for (k = i; k < N + 1; k++) {
            double temp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = temp;
        }

        // Make all rows below this one 0 in current column
        for (k = i + 1; k < N; k++) {
            double factor = -A[k][i] / A[i][i];
            for (j = i; j < N + 1; j++) {
                if (i == j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += factor * A[i][j];
                }
            }
        }
    }

    // Back substitution
    for (i = N - 1; i >= 0; i--) {
        x[i] = A[i][N] / A[i][i];
        for (j = i - 1; j >= 0; j--) {
            A[j][N] -= A[j][i] * x[i];
        }
    }
}

// Newton's method combined with Gaussian elimination
void newtonMethod(double x[N]) {
    int i, j, k;
    double f[N], J[N][N], dx[N], A[N][N + 1], x_new[N];

    for (i = 0; i < MAX_ITER; i++) {
        equations(x, f);
        jacobian(x, J);

        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                A[j][k] = J[j][k];
            }
            A[j][N] = -f[j];
        }

        gaussianElimination(A, dx);

        double max_error = 0.0;
        for (j = 0; j < N; j++) {
            double error = fabs(dx[j]);
            if (error > max_error) {
                max_error = error;
            }
        }

        if (max_error < EPSILON) {
            printf("Newton's method converged after %d iterations.\n", i + 1);
            return;
        }

        for (j = 0; j < N; j++) {
            x_new[j] = x[j] + dx[j];
        }

        for (j = 0; j < N; j++) {
            x[j] = x_new[j];
        }
    }

    printf("Newton's method did not converge after %d iterations.\n", MAX_ITER);
}

int main() {
    double x[N] = {0.0, 0.0, 0.0};

    printf("Enter the equation structure:\n");
    printf("Example: Equation 1: f(x1, x2, x3) = x1 + 2*x2 - x3 - 2\n");
    printf("Enter each equation in the specified format:\n");

    for (int i = 0; i < N; i++) {
        printf("Equation %d: ", i + 1);
        scanf("%lf %lf %lf", &x[i], &x[i + 1], &x[i + 2]);
    }

    newtonMethod(x);

    printf("Solution by Newton's method: x1 = %lf, x2 = %lf, x3 = %lf\n", x[0], x[1], x[2]);

    return 0;
}
