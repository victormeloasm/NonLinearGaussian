# Non-Linear Equation Solver

This program implements a solver for non-linear equations using the Jacobian matrix, Gaussian elimination, and Newton's method. It provides an efficient and robust approach to finding the solutions to systems of non-linear equations.

## Features

- Solves systems of non-linear equations with any number of variables.
- Implements the Jacobian matrix and Gaussian elimination to solve the linear system.
- Utilizes Newton's method to iteratively converge to the solutions.
- Provides convergence criteria to control the accuracy of the solutions.
- Allows customization of the equation structure and the derivative calculations.

## Requirements

- C compiler (e.g., GCC)

## Usage

1. Clone the repository or download the source code files.
2. Compile the code using a C compiler. For example, using GCC:
   ```
   gcc -o solver solver.c -lm
   ```
3. Run the executable file:
   ```
   ./solver
   ```
4. Follow the prompts to enter the equation structure and initial values.
5. The program will output the solution using Newton's method.

## Customization

- **Equation Structure:** Modify the `equations` function to define the system of equations. Update the equations according to your requirements.
- **Derivative Calculation:** Modify the `jacobian` function to calculate the derivative of each equation with respect to each variable. Adjust the derivatives based on the equation structure.
- **Convergence Criteria:** Modify the `MAX_ITER` and `EPSILON` constants in the code to control the maximum number of iterations and the convergence criterion for Newton's method.

## Example

Here's an example of how to use the program:

```c
// System of equations
void equations(double x[N], double f[N]) {
    f[0] = x[0] + 2 * x[1] - x[2] - 2;
    f[1] = x[0] * x[1] + x[1] * x[2] - 1;
    f[2] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 3;
}

// Jacobian matrix
void jacobian(double x[N], double J[N][N]) {
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

int main() {
    // Initialize variables and inputs

    // Call the solver function
    newtonMethod(x);

    // Output the solution

    return 0;
}
```

## License

This project is licensed under the [MIT License](LICENSE).
