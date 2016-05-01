#define DEFAULT_STEP 5 // Step size for printing points.

/**
 * for printing our matrix.
 */
#define prettyprint(x) do {fprintf(stdout, "%7.5g ", (x));} while(0)

/**
 * Prints a given n x m matrix.
 */
void printmat(size_t n, size_t m, double A[][m])
{
    size_t i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
             prettyprint(A[i][j]);
        }
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");    
}

/**
 * Prints results from finite difference approximation with time and spatial steps printed too.
 */
void printres(int interior, size_t nx, size_t nt, double dx, double dt, double A[][nt])
{
    size_t i, j;
    fprintf(stdout, "\ndx | dt ");
    for (i = 0; i < nt; i++) {prettyprint(i * dt);}

    for (i = 0; i < nx; i++)
    {
        fprintf(stdout, "\n");
        // If we only have interior points to print (saves making an array for them).
        prettyprint(interior ? (i * dx + dx) : (i * dx));
        for (j = 0; j < nt; j++)
        {
            prettyprint(A[i][j]);
        }
    }
}

/**
 * Prints only a few time steps
 */
void printfew(int interior, size_t printn, size_t nx, size_t nt, double dx, double dt, double A[][nt])
{
    size_t i, j, stepsize;
    stepsize = (printn && printn < nt) ? printn : DEFAULT_STEP;
   
    fprintf(stdout, "\ndx | dt ");
    for (i = 0; i < nt; i += stepsize) {prettyprint(i * dt);}

    for (i = 0; i < nx; i++)
    {
        fprintf(stdout, "\n");
        prettyprint(interior ? (i * dx + dx) : (i * dx));
        for (j = 0; j < nt; j += stepsize)
        {
            prettyprint(A[i][j]);
        }
    }
}

