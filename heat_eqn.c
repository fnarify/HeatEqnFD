/**
 * Simple application that allows the user to choose to plot a solution of the heat equation u_t = alpha * u_xx
 * using a finite difference scheme.
 * Currently the only two options are either Dirichlet B.C. explicit or Neumann implicit.
 *
 * The results can sometimes be a little odd, a good choice of values are as such:
 * alpha = 1
 * spatial steps = 10
 * no. timesteps = 200
 * timestep = 0.0005
 * For either Dirichlet or Neumann.
 *
 * Author - Bao Smith
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "heat_eqn.h"

#define DEFAULT_SPACE_STEP 10 // Default space step size for printing points.
#define DEFAULT_TIME_STEP  5  // Default time step size for printing.
#define PLOT_CALL_SIZE 64     // Size of string for gnuplot system call.

/**
 * Initial condition from fourier series:
 * F(x) = 1/2 + sum(1, inf)[2/(n*pi) * (sin(3*n*pi/4) - sin(n*pi/4)) * cos(n*pi*x)]
 */
#define initcond(x) (((x) >= 0.25) && ((x) <= 0.75))

/**
 * Memory error from malloc
 */
#define memerror(s) do {fprintf(stderr, "Not enough memory.\n"); exit(s);} while(0)

/**
 * Plots a 3d graph using gnuplot, where the data plotted has dimension data[row * col].
 * With the x-coordinates defined by dx (spatial-step) and the y co-ordinates defined by dt (time-step).
 * If there are interior points (Dirichlet case) then we have different x-coords.
 * First the data is output to outputf and then script defines the settings for gnuplot.
 */
void plot(char *outputf, char *script, size_t col, size_t row, int interior, double dx, double dt, double *data)
{
    size_t i, j;

    // Print shortened output to respective file, then plot
    FILE *output = fopen(outputf, "w");
    if (output == NULL)
    {
        fprintf(stderr, "Cannot write to file outputd\n");
        exit(1);
    }

    for (i = 0; i < row; i += DEFAULT_TIME_STEP)
    {
        for (j = 0; j < col; j += DEFAULT_SPACE_STEP)
        {
            fprintf(output, "%15g %15g %15g\n", i * dt, interior ? ((j + 1) * dx) : (j * dx), data[i * col + j]);
        }
        // Just makes the Dirichlet plot look a bit nicer.
        if (interior) {fprintf(output, "%15g %15g %15d\n", i * dt, (j + 1) * dx, 0);}

        fprintf(output, "\n");
    }

    fclose(output);
    free(data); // Don't need this data anymore        

    // Requires gnuplot to be on your user PATH variable.
    // Otherwise, uncomment the line below and add the file path to the bin folder of a gnuplot installation.
    char call[PLOT_CALL_SIZE];
    snprintf(call, PLOT_CALL_SIZE, "gnuplot -persistent %s", script);
    //system("cd YOUR_DIRECTORY_PATH\gnuplot\bin");
    system(call);
}

/**
 * Approximates the 1-d heat eqn u_t = alpha * u_xx with Dirichlet boundary conditions:
 * u(0, t) = u(L, t) = 0 using an explicit finite difference scheme.
 * nx - number of points in x-direction.
 * nt - number of time steps. 
 * dt - time step width.
 *
 * For this approximation to be bounded as n -> inf we need |1 - 4 * nu| < 1 -> dt <= dx^2/(2*alpha)
 * e.g. for dx = 0.1 we need dt <= 0.005 to make our approximation stable.
 * The same also applies for the implicit scheme.
 */
void expheateqn_d(double alpha, size_t nx, size_t nt, double dt)
{
    size_t i, j; 
    double *u_app;                 // Solution matrix (COLUMN-MAJOR order).
    double *upper, *diag, *lower;  // Diagonals of finite difference matrix, as the matrix is tri-diagonal.
    double *xi;                    // Interior pts.
    double dx;                     // Spatial step. 
    double bound;                  // To check we have a numerically stable result.
    double nu;
    
    dx = alpha * 1 / (double) (nx + 1);
    if (!dt) {dt = 0.4 * dx * dx;} // CFL condition.

    // For numerical stability of the given scheme.
    bound = dx * dx / (double) (2 * alpha);
    if (dt > bound) {dt = bound;}

    xi = malloc(sizeof(double) * nx);
    if (!xi) {memerror(1);}   
    // Only work with interior points for Dirichlet problem.
    for (i = 1; i < nx + 1; i++) {xi[i - 1] = i * dx;}

    upper = malloc(sizeof(double) * (nx - 1));
    diag = malloc(sizeof(double) * nx);
    lower = malloc(sizeof(double) * (nx - 1));
    if (!upper || !diag || !lower) {memerror(2);}

    // Finite-difference matrix setup for forward difference scheme (FTCS).
    nu = dt / (dx * dx);
    for (i = 0; i < nx; ) {diag[i++] = 1 - 2 * nu;}
    for (i = 0; i < nx - 1; i++) {upper[i] = lower[i] = nu;}

    u_app = calloc(nx * (nt + 1), sizeof(double));
    if (!u_app) {memerror(3);}

    // Apply initial condition to first column.
    for (i = 0; i < nx; i++) {u_app[i] = initcond(xi[i]);}

    free(xi);
    
    // Compute the value at the next time step by multiplying 
    // the previous time step by the FD matrix.
    for (i = 0; i < nt; i++)
    {
        u_app[(i + 1) * nx] = diag[0] * u_app[i * nx] + upper[0] * u_app[i * nx + 1];
        for (j = 1; j < nx - 1; j++)
        {
            u_app[(i + 1) * nx + j] = lower[j - 1] * u_app[i * nx + j - 1] + diag[j] * u_app[i * nx + j]
                                      + upper[j] * u_app[i * nx + j + 1];
        }
        u_app[(i + 1) * nx + j] = lower[j - 1] * u_app[i * nx + j - 1] + diag[j] * u_app[i * nx + j];
    }

    free(upper); free(diag); free(lower);

    plot("outputd", "plotd", nx, nt + 1, 1, dx, dt, u_app);
}

/**
 * Approximates the heat eqn u_t = alpha * u_xx with Neumann boundary conditions
 * u_x(0, t) = u_x(L, t) = 0 using an implicit finite difference scheme.
 */
void impheateqn_n(double alpha, size_t nx, size_t nt, double dt)
{
    size_t i, j; 
    double *u_app;                 // Solution matrix (COLUMN-MAJOR order).
    double *lower, *diag, *upper;  // Diagonals of the tri-diagonal finite-difference matrix.
    double dx;                     // Spatial step.
    double bound;                  // To check we have a numerically stable result.
    double nu;
    double *pivotvals;             // The value to divide each pivot of fd by. 
    double *curcol;                // Current time step of u_app we are working on.
    
    dx = alpha * 1 / (double) (nx + 1);    
    if (!dt) {dt = 0.4 * dx * dx;} // CFL condition.

    // If we aren't given a stable time step for our spatial step, force one.
    bound = dx * dx / (double) (2 * alpha);
    if (dt > bound) {dt = bound;}

    // Finite-difference matrix setup for backward difference scheme (BTCS).    
    lower = malloc(sizeof(double) * (nx + 1));
    diag = malloc(sizeof(double) * (nx + 2));
    upper = malloc(sizeof(double) * (nx + 1));
    if (!lower || ! diag || !upper) {memerror(1);}  

    nu = alpha * dt / (dx * dx);
    for (i = 0; i < nx + 2; ) {diag[i++] = 1 + 2 * nu;}
    for (i = 0; i < nx + 1; i++) {upper[i] = lower[i] =-nu;}
    // Changes to FD matrix given by Neumann BC.
    lower[nx] = upper[0] = -2 * nu;

    u_app = calloc((nx + 2) * (nt + 1), sizeof(double));
    if (!u_app) {memerror(2);}

    // Apply initial condition.
    for (i = 0; i < nx + 2; i++) {u_app[i] = initcond(i * dx);}

    /**
     * Elsewise given the system u_m = fd * u_m+1, we solve the system directly.
     * Since the matrix is tri-diagonal, and is clearly strictly diagonally dominant
     * we can apply Thomas' algorithm, and attain stable results.
     * Where the matrix is transformed such that it is an upper bi-diagonal matrix with
     * 1's down the main diagonal.
     */
    if ((nu < 0) && !(1 + 2 * nu)) 
    {
        fprintf(stderr, "choice of dt or dx creates a zero pivot, exiting program\n");
        exit(-1);
    }

    /* Thomas' algorithm. */

    pivotvals = malloc(sizeof(double) * (nx + 1));
    if (!pivotvals) {memerror(3);}   
    curcol = malloc(sizeof(double) * (nx + 2));
    if (!curcol) {memerror(4);}    

    // We can calculate the change in the upper diagonal beforehand as it does not depends on the u_m.
    pivotvals[0] = upper[0] / diag[0];
    for (i = 1; i < nx + 1; i++)
    {
        pivotvals[i] = upper[i] / (diag[i] - (lower[i - 1] * pivotvals[i - 1]));
    }

    free(upper);

    for (i = 0; i < nt; i++)
    {
        for (j = 0; j < nx + 2; j++) {curcol[j] = u_app[i * (nx + 2) + j];}

        // Calculate the change by our algorithm to the vector u_m.
        curcol[0] /= diag[0];            
        for (j = 1; j < nx + 2; j++)
        {
            curcol[j] = (curcol[j] - (lower[j - 1] * curcol[j - 1])) / (diag[j] - (lower[j - 1] * pivotvals[j - 1]));            
        }

        // Back substitution.
        j--;
        u_app[(i + 1) * (nx + 2) + j] = curcol[j];
        while (j-- > 0)
        {
            u_app[(i + 1) * (nx + 2) + j] = curcol[j];          
            u_app[(i + 1) * (nx + 2) + j] -= (pivotvals[j] * u_app[(i + 1) * (nx + 2) + j + 1]);
        }
    }

    free(lower); free(diag);
    free(pivotvals); free(curcol);

    plot("outputn", "plotn", nx + 2, nt + 1, 0, dx, dt, u_app);
}

int main()
{
    char scheme[2];
    size_t nx, nt;
    double dt;

    printf("This program solves the heat equation \nu_t = u_xx with either Dirichlet or Neumann boundary conditions.\n");
    printf("Do you want to solve for Dirichlet (D) or Neumann (N) scheme: ");
    scanf("%s", scheme);
    printf("Enter number of spatial steps (positive integer): ");
    scanf("%" PRSIZET, &nx);
    printf("Enter number of timesteps (positive integer): ");
    scanf("%" PRSIZET, &nt);
    printf("Enter timestep (positive rational number).\nThe value of the timestep is bounded by |1 - 4 *  dt / dx^2| : ");
    scanf("%lf", &dt);

    if (!strcmp(scheme, "D") || !strcmp(scheme, "d"))
        expheateqn_d(1, nx * DEFAULT_SPACE_STEP, nt * DEFAULT_TIME_STEP, dt);
    else if (!strcmp(scheme, "N") || !strcmp(scheme, "n"))
        impheateqn_n(1, nx * DEFAULT_SPACE_STEP, nt * DEFAULT_TIME_STEP, dt);

    return 0;
}

