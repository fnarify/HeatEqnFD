/**
 * Performs A * x where A and x are matrices.
 * Does this by multiplying each element of the first row of A by the respective element in the first
 * column of x, and again with the next column of x (Do this for each row of A).
 */
void matmul(size_t rA, size_t cA, size_t rx, size_t cx, double A[][cA], double x[][cx], double Ax[][cx])
{
    size_t i, j, k;
    double sum;

    if (cA != rx)
    {
        fprintf(stderr, "Column size of A must equal row size of x\n");
        return;
    }

    for (i = 0; i < rA; i++)
    {
        for (j = 0; j < cx; j++)
        {
            sum = 0.0;
            for (k = 0; k < cA; k++)
            {
                sum += A[i][k] * x[k][j];
            }
            Ax[i][j] = sum;
        }
    }
}

/**
 * Calculates determinent of a matrix m (assumed to be square).
 * Using the Laplace method.
 */
double det(size_t n, double A[][n])
{
    size_t i, j, k, subi, subj;
    double detA = 0;
    double submat[n - 1][n - 1];

    if (n == 2) 
    {
        return (A[0][0] * A[1][1] - A[1][0] * A[0][1]);
    }
    else if (n >= 3)
    {
        for (i = 0; i < n; i++)
        {
            subi = 0;

            // Collect correct points from A into our sub matrix.
            for (j = 1; j < n; j++)
            {
               subj = 0;
                for (k = 0; k < n; k++)
                {
                    if (k == i) {continue;}
                    submat[subi][subj++] = A[j][k];
                }
                subi++;
            }

            // Det calculation for each element in the first row.
            detA += A[0][i] * pow(-1, i + 2) * det(n - 1, submat);            
        }
    }
    else
    {
        return A[0][0];
    }

    return detA;
}

/**
 * Finds the co-factor matrix of A, assumes the matrix is square.
 */
void cofac(size_t n, double A[][n], double coA[][n])
{
    size_t i, j, k, l, coAi, coAj;    
    double submat[n - 1][n - 1];

    // Find the cofactor of each element in A.
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            coAi = 0;
            for (k = 0; k < n; k++)
            {
                coAj = 0;
                if (k == i) {continue;}
                for (l = 0; l < n; l++)
                {
                    if (l == j) {continue;}
                    submat[coAi][coAj++] = A[k][l];
                }
                coAi++;
            }
            coA[i][j] = pow(-1, i + j + 2) * det(n - 1, submat);
        }
    }
}

/**
 * Calculates the transpose in-place of a square matrix A.
 */
void transpose(size_t n, double A[][n])
{
    size_t i, j;
    double temp;

    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            temp = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = temp;
        }
    }
}

/**
 * Inverts an n x n matrix A, that is assumed to be pre-defined.
 * Finds the co-factor matrix then takes the transpose to obtain the adjoint matrix
 * and divides each element by the determinant of A.
 */
void inv(size_t n, double A[][n], double invA[][n])
{
    double detA;
    size_t i, j;
    // Sets invA to the the adjoint matrix of A.
    cofac(n, A, invA);
    transpose(n, invA);

    detA = det(n, A);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        { 
            invA[i][j] /= detA;
        }
    }
}

