#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double r8vec_dot ( int n, double a1[], double a2[] )
/*
    r8vec_dot computes the dot product of a pair of R8VEC's.
*/
{
    int i;
    double value;

    value = 0.0;
    for ( i = 0; i < n; i++ ) value = value + a1[i] * a2[i];
    return value;
}

void ax_st ( int n, int nz_num, int ia[], int ja[], double a[], double x[], double w[] )
/*
    ax_st computes A*x for a matrix stored in sparse triplet form.
*/
{
    int i;
    int j;
    int k;

    for ( i = 0; i < n; i++ ) w[i] = 0.0;

    for ( k = 0; k < nz_num; k++ )
    {
        i = ia[k];
        j = ja[k];
        w[i] = w[i] + a[k] * x[j];
    }

    return;
}

double **dmatrix( int nrl, int nrh, int ncl, int nch )
/*
    dmatrix allocates a double matrix.
*/
{
    int i;
    double **m;
    int nrow = nrh - nrl + 1;
    int ncol = nch - ncl + 1;

    /* 
    Allocate pointers to the rows.
    */
    m = ( double ** ) malloc ( (size_t) ( ( nrow + 1 ) * sizeof ( double* ) ) );
    m = m + 1;
    m = m - nrl;

    m[nrl] = ( double * ) malloc ( (size_t) ( ( nrow * ncol + 1 ) * sizeof ( double ) ) );
    m[nrl] = m[nrl] + 1;
    m[nrl] = m[nrl] - ncl;

    for ( i = nrl + 1; i <= nrh; i++ ) m[i] = m[i-1] + ncol;

    return m;
}

void mult_givens ( double c, double s, int k, double *g )
/*
    mult_givens applies a Givens rotation to two vector elements.
*/
{
  double g1;
  double g2;

  g1 = c * g[k] - s * g[k+1];
  g2 = s * g[k] + c * g[k+1];

  g[k]   = g1;
  g[k+1] = g2;

  return;
}

void free_dmatrix ( double **m, int nrl, int nrh, int ncl, int nch )
/*
    free_dmatrix frees a double matrix allocated by DMATRIX.
*/
{
  free ( ( char * ) ( m[nrl] + ncl - 1 ) );
  free ( ( char * ) ( m + nrl - 1 ) );

  return;
}

void mgmres_st (int n, int nz_num, int ia[], int ja[], double a[], double x[], double rhs[], int itr_max, int mr, double tol_abs, double tol_rel)
/*
    mgmres_st applies the restarted GMRES algorithm.
*/
{
    double av;
    double *c;
    double delta = 1.0e-03;
    double *g;
    double **h;
    double htmp;
    int i;
    int itr;
    int itr_used;
    int j;
    int k;
    int k_copy;
    double mu;
    double *r;
    double rho;
    double rho_tol;
    double *s;
    double **v;
    int verbose = 1;
    double *y;

    itr_used = 0;

    c = ( double * ) malloc ( mr * sizeof ( double ) );
    g = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );
    h = dmatrix(0,mr,0,mr-1);
    r = ( double * ) malloc ( n * sizeof ( double ) );
    s = ( double * ) malloc ( mr * sizeof ( double ) );
    v = dmatrix(0,mr,0,n-1);
    y = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) );

    for ( itr = 0; itr < itr_max; itr++ ) 
    {
        ax_st ( n, nz_num, ia, ja, a, x, r );

        for ( i = 0; i < n; i++ ) r[i] = rhs[i] - r[i];

        rho = sqrt ( r8vec_dot ( n, r, r ) );

        if ( verbose ) printf ( "  ITR = %8d  Residual = %e\n", itr, rho );

        if ( itr == 0 ) rho_tol = rho * tol_rel;

        for ( i = 0; i < n; i++ ) v[0][i] = r[i] / rho;

        g[0] = rho;
        for ( i = 1; i < mr + 1; i++ ) g[i] = 0.0;

        for ( i = 0; i < mr + 1; i++ ) for ( j = 0; j < mr; j++ ) h[i][j] = 0.0;

        for ( k = 0; k < mr; k++ )
        {
            k_copy = k;

            ax_st ( n, nz_num, ia, ja, a, v[k], v[k+1] );

            av = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );

            for ( j = 0; j < k+1; j++ )
            {
                h[j][k] = r8vec_dot ( n, v[k+1], v[j] );
                for ( i = 0; i < n; i++ ) v[k+1][i] = v[k+1][i] - h[j][k] * v[j][i];
            }

            h[k+1][k] = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );

            if ( ( av + delta * h[k+1][k] ) == av )
            {
                for ( j = 0; j < k+1; j++ )
                {
                    htmp = r8vec_dot ( n, v[k+1], v[j] );
                    h[j][k] = h[j][k] + htmp;
                    for ( i = 0; i < n; i++ ) v[k+1][i] = v[k+1][i] - htmp * v[j][i];
                }
                h[k+1][k] = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );
            }

            if ( h[k+1][k] != 0.0 )
            {
                for ( i = 0; i < n; i++ ) v[k+1][i] = v[k+1][i] / h[k+1][k];
            }

            if ( 0 < k )
            {
                for ( i = 0; i < k + 2; i++ ) y[i] = h[i][k];
                for ( j = 0; j < k; j++ ) mult_givens ( c[j], s[j], j, y );
                for ( i = 0; i < k + 2; i++ ) h[i][k] = y[i];
            }

            mu = sqrt ( h[k][k] * h[k][k] + h[k+1][k] * h[k+1][k] );
            c[k] = h[k][k] / mu;
            s[k] = -h[k+1][k] / mu;
            h[k][k] = c[k] * h[k][k] - s[k] * h[k+1][k];
            h[k+1][k] = 0.0;
            mult_givens ( c[k], s[k], k, g );

            rho = fabs ( g[k+1] );

            itr_used = itr_used + 1;

            if ( verbose ) printf ( "  K =   %8d  Residual = %e\n", k, rho );

            if ( rho <= rho_tol && rho <= tol_abs ) break;
        }

        k = k_copy;

        y[k] = g[k] / h[k][k];
        for ( i = k - 1; 0 <= i; i-- )
        {
            y[i] = g[i];
            for ( j = i+1; j < k + 1; j++ ) y[i] = y[i] - h[i][j] * y[j];
            y[i] = y[i] / h[i][i];
        }

        for ( i = 0; i < n; i++ ) for ( j = 0; j < k + 1; j++ ) x[i] = x[i] + v[j][i] * y[j];

        if ( rho <= rho_tol && rho <= tol_abs ) break;
    }

    if ( verbose )
    {
        printf ( "\n" );
        printf ( "MGMRES_ST:\n" );
        printf ( "  Iterations = %d\n", itr_used );
        printf ( "  Final residual = %e\n", rho );
    }
    /*
    Free memory.
    */
    free ( c );
    free ( g );
    free_dmatrix ( h, 0, mr, 0, mr - 1 );
    free ( r );
    free ( s );
    free_dmatrix ( v, 0, mr, 0, n - 1 );
    free ( y );

    return;
}

void solveGMRES(int n, int na, double *a, int *ia, int *ja, double *rhs, double *x)
{
    int mr = 5000;
    int inter_max = 10;

    mgmres_st(n, na, ia, ja, a, x, rhs, inter_max, mr, 1e-8, 1e-8);
}