/**
 * Simple application that allows the user to choose to plot a solution of the heat equation u_t = alpha * u_xx
 * using a finite difference scheme.
 * Currently the only two options are either Dirichlet B.C. explicit or Neumann implicit.
 *
 * This version differs from the original, by being more focused towards solving the problem exactly,
 * and making some of the details of the computation even less obvious than before.
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
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "heat_eqn.h"

#define DEFAULT_SPACE_STEP 10 // Default space step size for printing points.
#define DEFAULT_TIME_STEP  5  // Default time step size for printing.

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
 * Read error on scanf.
 */
#define readerror(s) do {fprintf(stderr, "Could not read input,\n"); exit(s);} while(0)

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

    errno = 0;
    // gnuplot should be on your user PATH variable, otherwise replace the first argument with the path to gnuplot/bin.
    execlp("gnuplot", "gnuplot", "-persistent", script, (char *) 0);
    if (errno) {fprintf(stderr, "%s", strerror(errno)); exit(-errno);}
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
    size_t i, j, colsize, rowsize; 
    double *u_app;                 // Solution matrix (COLUMN-MAJOR order).
    double diag, outdiag;          // Diagonals of finite difference matrix, as the matrix is tri-diagonal.
                                   // Since the diagonals of each matrix are composed of the same value, we can save storage this way.
    double dx;                     // Spatial step. 
    double bound;                  // To check we have a numerically stable result.
    double nu;                     // dt / dx^2
    
    dx = alpha * 1 / (double) (nx + 1);
    if (!dt) {dt = 0.4 * dx * dx;} // CFL condition.

    // For numerical stability of the given scheme.
    bound = dx * dx / (double) (2 * alpha);
    if (dt > bound) {dt = bound;}

    // Finite-difference matrix setup for forward difference scheme (FTCS).
    nu = dt / (dx * dx);
    diag = 1 - 2 * nu;
    outdiag = nu;

    // Initial setup for result matrix.
    colsize = nx; rowsize = nt + 1;
    u_app = calloc(colsize * rowsize, sizeof(double));
    if (!u_app) {memerror(1);}

    // Apply initial condition to first column.
    for (i = 0; i < colsize; i++) {u_app[i] = initcond((i + 1) * dx);}

    // Compute the value at the next time step by multiplying 
    // the previous time step by the FD matrix.
    for (i = 0; i < rowsize - 1; i++)
    {
        u_app[(i + 1) * colsize] = diag * u_app[i * colsize] + outdiag * u_app[i * colsize + 1];
        for (j = 1; j < colsize - 1; j++)
        {
            u_app[(i + 1) * colsize + j] = outdiag * u_app[i * colsize + j - 1] + diag * u_app[i * colsize + j]
                                      + outdiag * u_app[i * colsize + j + 1];
        }
        u_app[(i + 1) * colsize + j] = outdiag * u_app[i * colsize + j - 1] + diag * u_app[i * colsize + j];
    }

    plot("outputd", "plotd", colsize, rowsize, 1, dx, dt, u_app);
}

/**
 * Approximates the heat eqn u_t = alpha * u_xx with Neumann boundary conditions
 * u_x(0, t) = u_x(L, t) = 0 using an implicit finite difference scheme.
 */
void impheateqn_n(double alpha, size_t nx, size_t nt, double dt)
{
    size_t i, j, colsize, rowsize; 
    double *u_app;                 // Solution matrix (COLUMN-MAJOR order).
    double diag, *lower;           // Diagonals of the tri-diagonal finite-difference matrix.
                                   // Similar to the explicit scheme, we also store the change given by the B.C.
                                   // And since the upper diag is the the lower diag flipped.
    double dx;                     // Spatial step.
    double bound;                  // To check we have a numerically stable result.
    double nu;                     // dt / dx^2
    double *pivotvals;             // The value to divide each pivot of fd by. 
    double *curcol;                // Current time step of u_app we are working on.

    colsize = nx + 2; rowsize = nt + 1;    
    
    dx = alpha * 1 / (double) (nx + 1);    
    if (!dt) {dt = 0.4 * dx * dx;} // CFL condition.

    // If we aren't given a stable time step for our spatial step, force one.
    bound = dx * dx / (double) (2 * alpha);
    if (dt > bound) {dt = bound;}

    // Finite-difference matrix setup for backward difference scheme (BTCS).        
    lower = malloc(sizeof(double) * (colsize - 1));
    if (!lower) {memerror(1);}

    nu = alpha * dt / (dx * dx);
    diag = 1 + 2 * nu;
    for (i = 0; i < colsize - 2; i++) {lower[i] = -nu;}    
    // Changes to FD matrix given by Neumann BC.
    lower[i] = -2 * nu;

    u_app = calloc(colsize * rowsize, sizeof(double));
    if (!u_app) {memerror(2);}

    // Apply initial condition.
    for (i = 0; i < colsize; i++) {u_app[i] = initcond(i * dx);}

    /**
     * Elsewise given the system u_m = fd * u_m+1, we solve the system directly.
     * Since the matrix is tri-diagonal, and is clearly strictly diagonally dominant
     * we can apply Thomas' algorithm, and attain stable results.
     * Where the matrix is transformed such that it is an upper bi-diagonal matrix with
     * 1's down the main diagonal.
     */
    if ((nu < 0) && !diag) 
    {
        fprintf(stderr, "choice of dt or dx creates a zero pivot, exiting program\n");
        exit(-1);
    }

    /* Thomas' algorithm. */

    pivotvals = malloc(sizeof(double) * (colsize - 1));
    if (!pivotvals) {memerror(3);}   
    curcol = malloc(sizeof(double) * colsize);
    if (!curcol) {memerror(4);}    

    // We can calculate the change in the upper diagonal beforehand as it does not depends on the u_m.
    pivotvals[0] = lower[colsize - 2] / diag;
    for (i = 1; i < colsize - 1; i++)
    {
        pivotvals[i] = lower[colsize - 2 - i] / (diag - (lower[i - 1] * pivotvals[i - 1]));
    }

    for (i = 0; i < rowsize - 1; i++)
    {
        for (j = 0; j < colsize; j++) {curcol[j] = u_app[i * colsize + j];}

        // Calculate the change by our algorithm to the vector u_m.
        curcol[0] /= diag;            
        for (j = 1; j < colsize; j++)
        {
            curcol[j] = (curcol[j] - (lower[j - 1] * curcol[j - 1])) / (diag - (lower[j - 1] * pivotvals[j - 1]));            
        }

        // Back substitution.
        j--;
        u_app[(i + 1) * colsize + j] = curcol[j];
        while (j-- > 0)
        {
            u_app[(i + 1) * colsize + j] = curcol[j];          
            u_app[(i + 1) * colsize + j] -= (pivotvals[j] * u_app[(i + 1) * colsize + j + 1]);
        }
    }

    free(lower);
    free(pivotvals); free(curcol);

    plot("outputn", "plotn", colsize, rowsize, 0, dx, dt, u_app);
}

int main()
{
    char scheme[2];
    size_t nx, nt;
    double dt;

    printf("This program solves the heat equation \nu_t = u_xx with either Dirichlet or Neumann boundary conditions.\n");
    printf("Do you want to solve for Dirichlet (D) or Neumann (N) scheme: ");
    if (scanf("%1s", scheme) != 1) {readerror(-1);}
    scheme[1] = '\0';
    fflush(stdin); printf("Enter number of spatial steps (positive integer): ");
    if (scanf("%" PRSIZET, &nx) != 1) {readerror(-1);}
    fflush(stdin); printf("Enter number of timesteps (positive integer): ");
    if (scanf("%" PRSIZET, &nt) != 1) {readerror(-1);}
    fflush(stdin); printf("Enter timestep (positive rational number).\nThe value of the timestep is bounded by |1 - 4 *  dt / dx^2| : ");
    if (scanf("%lf", &dt) != 1) {readerror(-1);}

    if (!strcmp(scheme, "D") || !strcmp(scheme, "d"))
        expheateqn_d(1, nx * DEFAULT_SPACE_STEP, nt * DEFAULT_TIME_STEP, dt);
    else if (!strcmp(scheme, "N") || !strcmp(scheme, "n"))
        impheateqn_n(1, nx * DEFAULT_SPACE_STEP, nt * DEFAULT_TIME_STEP, dt);

    return 0;
}

