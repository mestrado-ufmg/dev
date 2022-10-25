#include <math.h>
#include <stdio.h>
#include <lapacke.h>

/*
#####################################################
    CONSTANTS
#####################################################
*/
const double ZERO_ERROR = 1e-8;
const double PI = 3.14159265359;
const double FACTOR = 1 / (4 * PI);
const int LAYERS = 300;
const double CTAU_CRIT = 1e-1;
const int LAMINAR_FLOW = 0;

/*
#####################################################
    STRUCTS
#####################################################
*/
struct Point
{
    double x;
    double y;
    double z;
};

struct VerticeConnection
{
    int n;
    double *coeffs;
    int *faces;
};

struct FacesConnection
{
    int n;
    int *faces;
};

struct ProfileParameters
{
    int n;
    double *eta;
    double *U;
    double *W;
    double *dU_deta;
    double *dW_deta;
    double *S;
    double *T;
    double *R;
};

struct FreestreamParameters
{
    double velocity;
    double density;
    double viscosity;
    double mach;
};

struct IntegralThicknessParameters
{
    double delta_1_ast;
    double delta_2_ast;
    double phi_11;
    double phi_12;
    double phi_21;
    double phi_22;
    double phi_1_ast;
    double phi_2_ast;
    double delta_1_line;
    double delta_2_line;
    double delta_q;
    double delta_q_o;
    double theta_1_o;
    double theta_2_o;
    double delta_1_o;
    double delta_2_o;
    double C_D;
    double C_D_x;
    double C_D_o;
    double C_f_1;
    double C_f_2;
    double theta_11;
    double theta_22;
};

struct IntegralDefectParameters
{
    double M_x;
    double M_y;
    double J_xx;
    double J_xy;
    double J_yx;
    double J_yy;
    double E_x;
    double E_y;
    double K_o_x;
    double K_o_y;
    double Q_x;
    double Q_y;
    double Q_o_x;
    double Q_o_y;
    double tau_w_x;
    double tau_w_y;
    double D;
    double D_x;
    double D_o;
    double K_tau_xx;
    double K_tau_xy;
    double K_tau_yx;
    double K_tau_yy;
    double S_tau_x;
    double S_tau_y;
};

struct EquationsParameters
{
    double M_x;
    double M_y;
    double J_xx;
    double J_xy;
    double J_yx;
    double J_yy;
    double E_x;
    double E_y;
    double K_o_x;
    double K_o_y;
    double Q_x;
    double Q_y;
    double Q_o_x;
    double Q_o_y;
    double tau_w_x;
    double tau_w_y;
    double D;
    double D_x;
    double D_o;
    double K_tau_xx;
    double K_tau_xy;
    double K_tau_yx;
    double K_tau_yy;
    double S_tau_x;
    double S_tau_y;
    double grad_q2_x, grad_q2_y;
    double grad_phi_x, grad_phi_y;
    double div_M;
    double div_J_x;
    double div_J_y;
    double div_E;
    double div_K_o;
    double div_K_tau_x;
    double div_K_tau_y;
    double vel;
    double density;
};

struct DivergentParameters
{
    double M;
    double J_x;
    double J_y;
    double E;
    double K_o;
    double K_tau_x;
    double K_tau_y;
};

struct GradientParameters
{
    double q2_x;
    double q2_y;
    double phi_x;
    double phi_y;
};

struct BoundaryLayerEquations
{
    double momentum_x;
    double momentum_y;
    double kinetic_energy;
    double lateral_curvature;
    double shear_stress_x;
    double shear_stress_y;
};

/*
#####################################################
    MATH FUNCTIONS
#####################################################
*/
double division(double a,
                double b)
{
    if ((-ZERO_ERROR < b) && (b < 0))
    {
        return -a / ZERO_ERROR;
    }
    else if ((0 < b) && (b < ZERO_ERROR))
    {
        return a / ZERO_ERROR;
    }
    else
    {
        return a / b;
    }
}

double norm(struct Point p)
{
    return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}

struct Point cross(struct Point p1,
                   struct Point p2)
{
    struct Point p = {p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z, p1.x * p2.y - p1.y * p2.x};
    return p;
}

double dot(struct Point p1,
           struct Point p2)
{
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

double angleBetweenVectors(struct Point p1,
                           struct Point p2)
{

    double norm1 = norm(p1);
    double norm2 = norm(p2);
    double dot12 = dot(p1, p2);

    return acos(dot12 / (norm1 * norm2));
}

double absValue(double a)
{
    if (a < 0)
    {
        return -a;
    }
    else
    {
        return a;
    }
}

void integrate_trap(int n,
                    double *x,
                    double *y,
                    double *out,
                    double mult)
{

    int i;

    *out = 0.0;

    for (i = 0; i < n - 1; i++)
    {
        *out = *out + 0.5 * (x[i + 1] - x[i]) * (y[i + 1] + y[i]);
    }

    *out = *out * mult;
}

/*
#####################################################
    HELPER FUNCTIONS
#####################################################
*/
void calculateVerticesConnection(int nv,
                                 int nf,
                                 double *vertices,
                                 int *faces,
                                 struct VerticeConnection *connection)
{

    /* Parameters */
    int i, j, k, n, faceLine, verticeLine1, verticeLine2, verticeLine3;
    int *facesIds;
    double *angles;
    struct Point point1, point2;
    double angle;
    double sum;

    /* Initialize */
    facesIds = (int *)malloc(50 * sizeof(int));
    angles = (double *)malloc(50 * sizeof(double));

    /* Loop over vertices */
    for (i = 0; i < nv; i++)
    {

        // Reset the number of faces and angles
        n = 0;
        sum = 0.0;

        /* Loop over faces */
        for (j = 0; j < nf; j++)
        {

            faceLine = j * 3;

            /* Check if the face contain the vertice */
            if ((faces[faceLine] == i) || (faces[faceLine + 1] == i) || (faces[faceLine + 2] == i))
            {

                /* Calculate the angle */
                if (faces[faceLine] == i)
                {
                    verticeLine1 = 3 * faces[faceLine];
                    verticeLine2 = 3 * faces[faceLine + 1];
                    verticeLine3 = 3 * faces[faceLine + 2];
                }
                else if (faces[faceLine + 1] == i)
                {
                    verticeLine3 = 3 * faces[faceLine];
                    verticeLine1 = 3 * faces[faceLine + 1];
                    verticeLine2 = 3 * faces[faceLine + 2];
                }
                else
                {
                    verticeLine2 = 3 * faces[faceLine];
                    verticeLine3 = 3 * faces[faceLine + 1];
                    verticeLine1 = 3 * faces[faceLine + 2];
                }

                point1.x = vertices[verticeLine2] - vertices[verticeLine1];
                point1.y = vertices[verticeLine2 + 1] - vertices[verticeLine1 + 1];
                point1.z = vertices[verticeLine2 + 2] - vertices[verticeLine1 + 2];

                point2.x = vertices[verticeLine3] - vertices[verticeLine1];
                point2.y = vertices[verticeLine3 + 1] - vertices[verticeLine1 + 1];
                point2.z = vertices[verticeLine3 + 2] - vertices[verticeLine1 + 2];

                angle = angleBetweenVectors(point1, point2);

                /* Increase the number of faces and save the angle */
                angles[n] = angle;
                facesIds[n] = j;
                sum = sum + angle;

                // Increase the number of faces
                n++;
            }
        }

        /* Save the number of faces, faces id's and angles in connection */
        connection[i].n = n;
        connection[i].coeffs = (double *)malloc(n * sizeof(double));
        connection[i].faces = (int *)malloc(n * sizeof(int));

        for (k = 0; k < n; k++)
        {
            connection[i].coeffs[k] = angles[k] / sum;
            connection[i].faces[k] = facesIds[k];
        }
    }

    free(facesIds);
    free(angles);
}

void calculateFacesConnection(int nv,
                              int nf,
                              int *faces,
                              struct VerticeConnection *verticesConnection,
                              struct FacesConnection *facesConnection)
{

    /* Parameters */
    int i, j, k;
    int check;
    int index;
    int *connectedFaces;
    int n;

    /* Initialize */
    connectedFaces = (int *)malloc(50 * sizeof(int));

    for (i = 0; i < nf; i++)
    {

        n = 0;

        // First index
        index = faces[3 * i];

        for (j = 0; j < verticesConnection[index].n; j++)
        {
            if (verticesConnection[index].faces[j] != i)
            {
                connectedFaces[n] = verticesConnection[index].faces[j];
                n = n + 1;
            }
        }

        // Second index
        index = faces[3 * i + 1];

        for (j = 0; j < verticesConnection[index].n; j++)
        {

            check = 0;

            for (k = 0; k < n; k++)
            {
                if (connectedFaces[k] == verticesConnection[index].faces[j])
                {
                    check = 1;
                    break;
                }
            }

            if ((verticesConnection[index].faces[j] != i) && (check == 0))
            {
                connectedFaces[n] = verticesConnection[index].faces[j];
                n = n + 1;
            }
        }

        // Second index
        index = faces[3 * i + 2];

        for (j = 0; j < verticesConnection[index].n; j++)
        {

            check = 0;

            for (k = 0; k < n; k++)
            {
                if (connectedFaces[k] == verticesConnection[index].faces[j])
                {
                    check = 1;
                    break;
                }
            }

            if ((verticesConnection[index].faces[j] != i) && (check == 0))
            {
                connectedFaces[n] = verticesConnection[index].faces[j];
                n = n + 1;
            }
        }

        facesConnection[i].n = n;
        facesConnection[i].faces = (int *)malloc(n * sizeof(int));

        for (j = 0; j < n; j++)
        {
            facesConnection[i].faces[j] = connectedFaces[j];
        }
    }

    free(connectedFaces);
}

void calculateVerticesValues(int nv,
                             int nf,
                             struct VerticeConnection *connection,
                             double *in,
                             double *out)
{

    int i, j;

    for (i = 0; i < nv; i++)
    {
        out[i] = 0;
        for (j = 0; j < connection[i].n; j++)
        {
            out[i] = out[i] + in[connection[i].faces[j]] * connection[i].coeffs[j];
        }
    }
}

void find_exp_ratio(double S,
                    double a0,
                    double n,
                    double *r0)
{

    /* Parameters */
    int i;
    int calc;
    double step;
    double tol;
    double aux, faux;
    double dfdr, dx;
    double a, fa, b, fb;

    calc = 0;

    // Convergence
    tol = 1e-8;

    // Check convergence
    if ((1 - tol < *r0) && (*r0 < 1 + tol))
        *r0 = 1 + tol;

    aux = a0 * (1 - pow(*r0, n)) / (1 - *r0) - S;

    if ((-tol < aux) && (aux < tol))
        calc = 1;

    if (calc == 0)
    {

        // Find edges
        step = 1e-2;

        if (aux < 0)
        {
            fa = aux;
            a = *r0;
            while (aux <= 0)
            {
                *r0 = *r0 + step;
                if ((1 - tol < *r0) && (*r0 < 1 + tol))
                    *r0 = 1 + tol;
                aux = a0 * (1 - pow(*r0, n)) / (1 - *r0) - S;
                if (aux < 0)
                {
                    fa = aux;
                    a = *r0;
                }
            }
            b = *r0;
            fb = aux;
        }
        else
        {
            fb = aux;
            b = *r0;
            while (aux >= 0)
            {
                *r0 = *r0 - step;
                if ((1 - tol < *r0) && (*r0 < 1 + tol))
                    *r0 = 1 + tol;
                aux = a0 * (1 - pow(*r0, n)) / (1 - *r0) - S;
                if (aux > 0)
                {
                    fb = aux;
                    b = *r0;
                }
            }
            a = *r0;
            fa = aux;
        }

        // Find root
        for (i = 0; i < 500; i++)
        {

            aux = 0.5 * (a + b);
            if ((1 - tol < aux) && (aux < 1 + tol))
                aux = 1 + tol;
            faux = a0 * (1 - pow(aux, n)) / (1 - aux) - S;

            *r0 = aux;

            if ((-tol < faux) && (faux < tol))
            {
                break;
            }

            if (faux < 0)
            {
                a = aux;
                fa = faux;
            }
            else
            {
                b = aux;
                fb = faux;
            }
        }
    }
}

struct Point gradient(double centerValue,
                      struct Point p1,
                      struct Point p2,
                      struct Point p3,
                      struct Point e1,
                      struct Point e2,
                      struct Point e3,
                      struct Point vel,
                      double transpiration)
{

    struct Point p0;

    struct Point p01;
    struct Point p02;
    struct Point p03;

    struct Point n1;
    struct Point n2;
    struct Point n3;
    struct Point n;

    struct Point grad;

    // Gradient in plane system
    p0.x = (p1.x + p2.x + p3.x) / 3;
    p0.y = (p1.y + p2.y + p3.y) / 3;
    p0.z = centerValue;

    p01.x = p1.x - p0.x;
    p01.y = p1.y - p0.y;
    p01.z = p1.z - p0.z;

    p02.x = p2.x - p0.x;
    p02.y = p2.y - p0.y;
    p02.z = p2.z - p0.z;

    p03.x = p3.x - p0.x;
    p03.y = p3.y - p0.y;
    p03.z = p3.z - p0.z;

    n1.x = p01.y * p02.z - p01.z * p02.y;
    n1.y = p01.z * p02.x - p01.x * p02.z;
    n1.z = p01.x * p02.y - p01.y * p02.x;

    n2.x = p02.y * p03.z - p02.z * p03.y;
    n2.y = p02.z * p03.x - p02.x * p03.z;
    n2.z = p02.x * p03.y - p02.y * p03.x;

    n3.x = p03.y * p01.z - p03.z * p01.y;
    n3.y = p03.z * p01.x - p03.x * p01.z;
    n3.z = p03.x * p01.y - p03.y * p01.x;

    n.x = n1.x + n2.x + n3.x;
    n.y = n1.y + n2.y + n3.y;
    n.z = n1.z + n2.z + n3.z;

    grad.x = -n.x / n.z;
    grad.y = -n.y / n.z;
    grad.z = 0.0;

    // Convert to streamline system
    struct Point dir1;
    struct Point dir2;
    double norm;
    double s1;
    double s2;

    norm = sqrt(pow(vel.x - transpiration * e3.x, 2) + pow(vel.y - transpiration * e3.y, 2) + pow(vel.z - transpiration * e3.z, 2));

    dir1.x = (vel.x - transpiration * e3.x) / norm;
    dir1.y = (vel.y - transpiration * e3.y) / norm;
    dir1.z = (vel.z - transpiration * e3.z) / norm;

    dir2.x = e3.y * dir1.z - e3.z * dir1.y;
    dir2.y = e3.z * dir1.x - e3.x * dir1.z;
    dir2.z = e3.x * dir1.y - e3.y * dir1.x;

    s1 = grad.x * dir1.x + grad.y * dir1.y + grad.z * dir1.z;
    s2 = grad.x * dir2.x + grad.y * dir2.y + grad.z * dir2.z;

    // Output
    struct Point out;

    out.x = s1;
    out.y = s2;
    out.z = 0.0;

    return out;
}

double divergence(struct Point p1,
                  struct Point p2,
                  struct Point p3,
                  struct Point v1,
                  struct Point v2,
                  struct Point v3,
                  double area)
{

    // Sides
    double norm1;
    double norm2;
    double norm3;

    struct Point line1;
    struct Point line2;
    struct Point line3;

    struct Point vec1;
    struct Point vec2;
    struct Point vec3;

    norm1 = sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
    norm2 = sqrt(pow(p3.x - p2.x, 2) + pow(p3.y - p2.y, 2));
    norm3 = sqrt(pow(p1.x - p3.x, 2) + pow(p1.y - p3.y, 2));

    line1.x = (p2.x - p1.x) / norm1;
    line1.y = (p2.y - p1.y) / norm1;
    line1.z = 0.0;
    line2.x = (p3.x - p2.x) / norm2;
    line2.y = (p3.y - p2.y) / norm2;
    line2.z = 0.0;
    line3.x = (p1.x - p3.x) / norm3;
    line3.y = (p1.y - p3.y) / norm3;
    line3.z = 0.0;

    vec1.x = line1.y;
    vec1.y = -line1.x;
    vec1.z = 0.0;

    vec2.x = line2.y;
    vec2.y = -line2.x;
    vec2.z = 0.0;

    vec3.x = line3.y;
    vec3.y = -line3.x;
    vec3.z = 0.0;

    // Flux
    double a, b;
    double integral;

    // Line 1
    a = vec1.x * v1.x + vec1.y * v1.y;
    b = vec1.x * v2.x + vec1.y * v2.y;
    integral = 0.5 * norm1 * (a + b);

    // Line 2
    a = vec2.x * v2.x + vec2.y * v2.y;
    b = vec2.x * v3.x + vec2.y * v3.y;
    integral = integral + 0.5 * norm2 * (a + b);

    // Line 3
    a = vec3.x * v3.x + vec3.y * v3.y;
    b = vec3.x * v1.x + vec3.y * v1.y;
    integral = integral + 0.5 * norm3 * (a + b);

    return integral / area;
}

void calculateGradient(double centerValue,
                       struct Point p1,
                       struct Point p2,
                       struct Point p3,
                       struct Point e1,
                       struct Point e2,
                       struct Point e3,
                       struct Point vel,
                       double transpiration,
                       double *x, double *y)
{

    struct Point p0;

    struct Point p01;
    struct Point p02;
    struct Point p03;

    struct Point n1;
    struct Point n2;
    struct Point n3;
    struct Point n;

    struct Point grad;

    // Gradient in plane system
    p0.x = (p1.x + p2.x + p3.x) / 3;
    p0.y = (p1.y + p2.y + p3.y) / 3;
    p0.z = centerValue;

    p01.x = p1.x - p0.x;
    p01.y = p1.y - p0.y;
    p01.z = p1.z - p0.z;

    p02.x = p2.x - p0.x;
    p02.y = p2.y - p0.y;
    p02.z = p2.z - p0.z;

    p03.x = p3.x - p0.x;
    p03.y = p3.y - p0.y;
    p03.z = p3.z - p0.z;

    n1.x = p01.y * p02.z - p01.z * p02.y;
    n1.y = p01.z * p02.x - p01.x * p02.z;
    n1.z = p01.x * p02.y - p01.y * p02.x;

    n2.x = p02.y * p03.z - p02.z * p03.y;
    n2.y = p02.z * p03.x - p02.x * p03.z;
    n2.z = p02.x * p03.y - p02.y * p03.x;

    n3.x = p03.y * p01.z - p03.z * p01.y;
    n3.y = p03.z * p01.x - p03.x * p01.z;
    n3.z = p03.x * p01.y - p03.y * p01.x;

    n.x = n1.x + n2.x + n3.x;
    n.y = n1.y + n2.y + n3.y;
    n.z = n1.z + n2.z + n3.z;

    grad.x = -n.x / n.z;
    grad.y = -n.y / n.z;
    grad.z = 0.0;

    // Convert to streamline system
    struct Point dir1;
    struct Point dir2;
    double norm;
    double s1;
    double s2;

    norm = sqrt(pow(vel.x - transpiration * e3.x, 2) + pow(vel.y - transpiration * e3.y, 2) + pow(vel.z - transpiration * e3.z, 2));

    dir1.x = (vel.x - transpiration * e3.x) / norm;
    dir1.y = (vel.y - transpiration * e3.y) / norm;
    dir1.z = (vel.z - transpiration * e3.z) / norm;

    dir2.x = e3.y * dir1.z - e3.z * dir1.y;
    dir2.y = e3.z * dir1.x - e3.x * dir1.z;
    dir2.z = e3.x * dir1.y - e3.y * dir1.x;

    s1 = grad.x * dir1.x + grad.y * dir1.y + grad.z * dir1.z;
    s2 = grad.x * dir2.x + grad.y * dir2.y + grad.z * dir2.z;

    // Output
    *x = s1;
    *x = s2;
}

void calculateDivergence(struct Point p1,
                         struct Point p2,
                         struct Point p3,
                         struct Point v1,
                         struct Point v2,
                         struct Point v3,
                         double area,
                         double *out) {

    // Sides
    double norm1;
    double norm2;
    double norm3;

    struct Point line1;
    struct Point line2;
    struct Point line3;

    struct Point vec1;
    struct Point vec2;
    struct Point vec3;

    norm1 = sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
    norm2 = sqrt(pow(p3.x - p2.x, 2) + pow(p3.y - p2.y, 2));
    norm3 = sqrt(pow(p1.x - p3.x, 2) + pow(p1.y - p3.y, 2));

    line1.x = (p2.x - p1.x) / norm1;
    line1.y = (p2.y - p1.y) / norm1;
    line1.z = 0.0;
    line2.x = (p3.x - p2.x) / norm2;
    line2.y = (p3.y - p2.y) / norm2;
    line2.z = 0.0;
    line3.x = (p1.x - p3.x) / norm3;
    line3.y = (p1.y - p3.y) / norm3;
    line3.z = 0.0;

    vec1.x = line1.y;
    vec1.y = -line1.x;
    vec1.z = 0.0;

    vec2.x = line2.y;
    vec2.y = -line2.x;
    vec2.z = 0.0;

    vec3.x = line3.y;
    vec3.y = -line3.x;
    vec3.z = 0.0;

    // Flux
    double a, b;
    double integral;

    // Line 1
    a = vec1.x * v1.x + vec1.y * v1.y;
    b = vec1.x * v2.x + vec1.y * v2.y;
    integral = 0.5 * norm1 * (a + b);

    // Line 2
    a = vec2.x * v2.x + vec2.y * v2.y;
    b = vec2.x * v3.x + vec2.y * v3.y;
    integral = integral + 0.5 * norm2 * (a + b);

    // Line 3
    a = vec3.x * v3.x + vec3.y * v3.y;
    b = vec3.x * v1.x + vec3.y * v1.y;
    integral = integral + 0.5 * norm3 * (a + b);

    *out = integral / area;
}

/*
#####################################################
    LINEAR SYSTEM
#####################################################
*/
void solveLinearSystem(lapack_int n,
                       double *A,
                       double *b)
{

    /* Locals */
    lapack_int info;

    /* Local arrays */
    lapack_int *ipiv;

    /* Initialization */
    ipiv = (lapack_int *)malloc(n * sizeof(lapack_int));

    /* Solve the equations A*X = B */
    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 1, A, n, ipiv, b, 1);
}

/*
#####################################################
    POTENTIAL FLOW
#####################################################
*/
void sourceFunc(struct Point p,
                struct Point p1,
                struct Point p2,
                struct Point p3,
                struct Point e1,
                struct Point e2,
                struct Point e3,
                double area,
                double maxDistance,
                double *vel)
{

    double u, v, w;
    double distance = norm(p);

    if (distance > maxDistance)
    {

        double pNorm3 = pow(norm(p), 3);

        u = FACTOR * area * p.x / pNorm3;
        v = FACTOR * area * p.y / pNorm3;
        w = FACTOR * area * p.z / pNorm3;
    }
    else
    {

        // Velocity parameters
        double r1, r2, r3;
        double l1, l2, l3;
        double h1, h2, h3;

        double d12, d23, d31;
        double m12, m23, m31;
        double ln12, ln23, ln31;

        // Calculate
        r1 = sqrt(pow(p.x - p1.x, 2) + pow(p.y - p1.y, 2) + pow(p.z, 2));
        r2 = sqrt(pow(p.x - p2.x, 2) + pow(p.y - p2.y, 2) + pow(p.z, 2));
        r3 = sqrt(pow(p.x - p3.x, 2) + pow(p.y - p3.y, 2) + pow(p.z, 2));

        l1 = pow(p.x - p1.x, 2) + pow(p.z, 2);
        l2 = pow(p.x - p2.x, 2) + pow(p.z, 2);
        l3 = pow(p.x - p3.x, 2) + pow(p.z, 2);

        h1 = (p.x - p1.x) * (p.y - p1.y);
        h2 = (p.x - p2.x) * (p.y - p2.y);
        h3 = (p.x - p3.x) * (p.y - p3.y);

        d12 = sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
        m12 = division(p2.y - p1.y, p2.x - p1.x);

        d23 = sqrt(pow(p3.x - p2.x, 2) + pow(p3.y - p2.y, 2));
        m23 = division(p3.y - p2.y, p3.x - p2.x);

        d31 = sqrt(pow(p1.x - p3.x, 2) + pow(p1.y - p3.y, 2));
        m31 = division(p1.y - p3.y, p1.x - p3.x);

        ln12 = log(division(r1 + r2 - d12, r1 + r2 + d12));
        ln23 = log(division(r2 + r3 - d23, r2 + r3 + d23));
        ln31 = log(division(r3 + r1 - d31, r3 + r1 + d31));

        u = -FACTOR * (division((p2.y - p1.y), d12) * ln12 + division((p3.y - p2.y), d23) * ln23 + division((p1.y - p3.y), d31) * ln31);
        v = FACTOR * (division((p2.x - p1.x), d12) * ln12 + division((p3.x - p2.x), d23) * ln23 + division((p1.x - p3.x), d31) * ln31);
        w = -FACTOR * (atan(division(m12 * l1 - h1, p.z * r1)) - atan(division(m12 * l2 - h2, p.z * r2)) + atan(division(m23 * l2 - h2, p.z * r2)) - atan(division(m23 * l3 - h3, p.z * r3)) + atan(division(m31 * l3 - h3, p.z * r3)) - atan(division(m31 * l1 - h1, p.z * r1)));
    }

    vel[0] = u * e1.x + v * e2.x + w * e3.x;
    vel[1] = u * e1.y + v * e2.y + w * e3.y;
    vel[2] = u * e1.z + v * e2.z + w * e3.z;
}

void lineFunc(struct Point p,
              struct Point p1,
              struct Point p2,
              double *vel)
{

    struct Point r1 = {p1.x - p.x, p1.y - p.y, p1.z - p.z};
    struct Point r2 = {p2.x - p.x, p2.y - p.y, p2.z - p.z};

    struct Point r1xr2 = cross(r1, r2);

    double r1Norm = norm(r1);
    double r2Norm = norm(r2);

    double r1xr2Norm2 = pow(norm(r1xr2), 2);

    double dot = (1 / r1xr2Norm2) * ((r1.x - r2.x) * (r1.x / r1Norm - r2.x / r2Norm) + (r1.y - r2.y) * (r1.y / r1Norm - r2.y / r2Norm) + (r1.z - r2.z) * (r1.z / r1Norm - r2.z / r2Norm));

    vel[0] = FACTOR * r1xr2.x * dot;
    vel[1] = FACTOR * r1xr2.y * dot;
    vel[2] = FACTOR * r1xr2.z * dot;
}

void doubletFunc(struct Point p,
                 struct Point p1,
                 struct Point p2,
                 struct Point p3,
                 struct Point e1,
                 struct Point e2,
                 struct Point e3,
                 double area,
                 double maxDistance,
                 double *vel)
{

    double distance = norm(p);

    if (distance > maxDistance)
    {

        double pxLocal = p.x * e1.x + p.y * e1.y + p.z * e1.z;
        double pyLocal = p.x * e2.x + p.y * e2.y + p.z * e2.z;
        double pzLocal = p.x * e3.x + p.y * e3.y + p.z * e3.z;
        double den = pow(pxLocal * pxLocal + pyLocal * pyLocal + pzLocal * pzLocal, 2.5);

        double u = 0.75 * FACTOR * area * pzLocal * pxLocal / den;
        double v = 0.75 * FACTOR * area * pzLocal * pyLocal / den;
        double w = -FACTOR * area * (pxLocal * pxLocal + pyLocal * pyLocal - 2 * pzLocal * pzLocal) / den;

        vel[0] = u * e1.x + v * e2.x + w * e3.x;
        vel[1] = u * e1.y + v * e2.y + w * e3.y;
        vel[2] = u * e1.z + v * e2.z + w * e3.z;
    }
    else
    {

        double u, v, w;

        double *vel1 = (double *)malloc(3 * sizeof(double *));
        double *vel2 = (double *)malloc(3 * sizeof(double *));
        double *vel3 = (double *)malloc(3 * sizeof(double *));

        lineFunc(p, p1, p2, vel1);
        lineFunc(p, p2, p3, vel2);
        lineFunc(p, p3, p1, vel3);

        u = vel1[0] + vel2[0] + vel3[0];
        v = vel1[1] + vel2[1] + vel3[1];
        w = vel1[2] + vel2[2] + vel3[2];

        vel[0] = u * e1.x + v * e2.x + w * e3.x;
        vel[1] = u * e1.y + v * e2.y + w * e3.y;
        vel[2] = u * e1.z + v * e2.z + w * e3.z;

        free(vel1);
        free(vel2);
        free(vel3);
    }
}

void createLinearSystem(int n,
                        double *facesAreas,
                        double *facesMaxDistance,
                        double *facesCenter,
                        double *controlPoints,
                        double *p1,
                        double *p2,
                        double *p3,
                        double *e1,
                        double *e2,
                        double *e3,
                        double *freestream,
                        double *sigma,
                        int nSpanLeftWing,
                        int nWakeLeftWing,
                        int *leftWingGrid,
                        double *leftWingVertices,
                        int *leftWingFaces,
                        int nSpanRightWing,
                        int nWakeRightWing,
                        int *rightWingGrid,
                        double *rightWingVertices,
                        int *rightWingFaces,
                        int nSpanTail,
                        int nWakeTail,
                        int *tailGrid,
                        double *tailVertices,
                        int *tailFaces,
                        double *matrix,
                        double *array,
                        double *matrixVelx,
                        double *matrixVely,
                        double *matrixVelz,
                        double *arrayVel)
{

    ///---------------------------------------///
    /// Parameters

    // Loops
    int i, j, k, l;

    // Point
    int i3D1, i3D2, i3D3, j3D1, j3D2, j3D3;
    int j2D1, j2D2;
    int indexLine1, indexLine2;
    struct Point p;
    struct Point p1Line;
    struct Point p2Line;
    struct Point pLocal;
    struct Point p1Local;
    p1Local.z = 0.0;
    struct Point p2Local;
    p2Local.z = 0.0;
    struct Point p3Local;
    p3Local.z = 0.0;

    // Base vectors
    struct Point e3iPoint;

    struct Point e1jPoint;
    struct Point e2jPoint;
    struct Point e3jPoint;

    // Velocities
    double *sourceVel = (double *)malloc(3 * sizeof(double));
    double *doubletVel = (double *)malloc(3 * sizeof(double));
    double *lineVel = (double *)malloc(3 * sizeof(double));

    ///---------------------------------------///
    /// Linear system

    // Create
    for (i = 0; i < n; i++)
    {

        i3D1 = i * 3;
        i3D2 = i3D1 + 1;
        i3D3 = i3D1 + 2;

        e3iPoint.x = e3[i3D1];
        e3iPoint.y = e3[i3D2];
        e3iPoint.z = e3[i3D3];

        array[i] = 0;
        arrayVel[i * 3] = 0;
        arrayVel[i * 3 + 1] = 0;
        arrayVel[i * 3 + 2] = 0;

        // Surface
        // Effect of j on i
        for (j = 0; j < n; j++)
        {

            j3D1 = j * 3;
            j3D2 = j3D1 + 1;
            j3D3 = j3D1 + 2;

            j2D1 = j * 2;
            j2D2 = j2D1 + 1;

            // Points
            e1jPoint.x = e1[j3D1];
            e1jPoint.y = e1[j3D2];
            e1jPoint.z = e1[j3D3];
            e2jPoint.x = e2[j3D1];
            e2jPoint.y = e2[j3D2];
            e2jPoint.z = e2[j3D3];
            e3jPoint.x = e3[j3D1];
            e3jPoint.y = e3[j3D2];
            e3jPoint.z = e3[j3D3];

            p.x = controlPoints[i3D1] - facesCenter[j3D1];
            p.y = controlPoints[i3D2] - facesCenter[j3D2];
            p.z = controlPoints[i3D3] - facesCenter[j3D3];

            pLocal.x = p.x * e1jPoint.x + p.y * e1jPoint.y + p.z * e1jPoint.z;
            pLocal.y = p.x * e2jPoint.x + p.y * e2jPoint.y + p.z * e2jPoint.z;
            pLocal.z = p.x * e3jPoint.x + p.y * e3jPoint.y + p.z * e3jPoint.z;

            p1Local.x = p1[j2D1];
            p1Local.y = p1[j2D2];
            p2Local.x = p2[j2D1];
            p2Local.y = p2[j2D2];
            p3Local.x = p3[j2D1];
            p3Local.y = p3[j2D2];

            sourceFunc(pLocal, p1Local, p2Local, p3Local, e1jPoint, e2jPoint, e3jPoint, facesAreas[j], facesMaxDistance[j], sourceVel);
            doubletFunc(pLocal, p1Local, p2Local, p3Local, e1jPoint, e2jPoint, e3jPoint, facesAreas[j], facesMaxDistance[j], doubletVel);

            matrix[i * n + j] = doubletVel[0] * e3iPoint.x + doubletVel[1] * e3iPoint.y + doubletVel[2] * e3iPoint.z;
            array[i] = array[i] - sigma[j] * (sourceVel[0] * e3iPoint.x + sourceVel[1] * e3iPoint.y + sourceVel[2] * e3iPoint.z);

            matrixVelx[i * n + j] = doubletVel[0];
            matrixVely[i * n + j] = doubletVel[1];
            matrixVelz[i * n + j] = doubletVel[2];

            arrayVel[i * 3] = arrayVel[i * 3] + sigma[j] * sourceVel[0];
            arrayVel[i * 3 + 1] = arrayVel[i * 3 + 1] + sigma[j] * sourceVel[1];
            arrayVel[i * 3 + 2] = arrayVel[i * 3 + 2] + sigma[j] * sourceVel[2];
        }

        array[i] = array[i] - (freestream[0] * e3iPoint.x + freestream[1] * e3iPoint.y + freestream[2] * e3iPoint.z);

        arrayVel[i * 3] = arrayVel[i * 3] + freestream[0];
        arrayVel[i * 3 + 1] = arrayVel[i * 3 + 1] + freestream[1];
        arrayVel[i * 3 + 2] = arrayVel[i * 3 + 2] + freestream[2];

        // Wake
        p.x = controlPoints[i3D1];
        p.y = controlPoints[i3D2];
        p.z = controlPoints[i3D3];

        // Left wing
        for (k = 0; k < nSpanLeftWing; k++)
        {

            for (l = 0; l < nWakeLeftWing - 1; l++)
            {

                indexLine1 = leftWingGrid[k * nWakeLeftWing + l] * 3;
                p1Line.x = leftWingVertices[indexLine1];
                p1Line.y = leftWingVertices[indexLine1 + 1];
                p1Line.z = leftWingVertices[indexLine1 + 2];

                indexLine2 = leftWingGrid[k * nWakeLeftWing + l + 1] * 3;
                p2Line.x = leftWingVertices[indexLine2];
                p2Line.y = leftWingVertices[indexLine2 + 1];
                p2Line.z = leftWingVertices[indexLine2 + 2];

                lineFunc(p, p1Line, p2Line, lineVel);

                if (k == 0)
                {

                    matrix[i * n + leftWingFaces[k * 2]] = matrix[i * n + leftWingFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[k * 2 + 1]] = matrix[i * n + leftWingFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + leftWingFaces[k * 2]] = matrixVelx[i * n + leftWingFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + leftWingFaces[k * 2]] = matrixVely[i * n + leftWingFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + leftWingFaces[k * 2]] = matrixVelz[i * n + leftWingFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + leftWingFaces[k * 2 + 1]] = matrixVelx[i * n + leftWingFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + leftWingFaces[k * 2 + 1]] = matrixVely[i * n + leftWingFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + leftWingFaces[k * 2 + 1]] = matrixVelz[i * n + leftWingFaces[k * 2 + 1]] + lineVel[2];
                }
                else if (k == nSpanLeftWing - 1)
                {

                    matrix[i * n + leftWingFaces[(k - 1) * 2]] = matrix[i * n + leftWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + leftWingFaces[(k - 1) * 2]] = matrixVelx[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + leftWingFaces[(k - 1) * 2]] = matrixVely[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + leftWingFaces[(k - 1) * 2]] = matrixVelz[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[2];
                }
                else
                {

                    matrix[i * n + leftWingFaces[(k - 1) * 2]] = matrix[i * n + leftWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + leftWingFaces[(k - 1) * 2]] = matrixVelx[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + leftWingFaces[(k - 1) * 2]] = matrixVely[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + leftWingFaces[(k - 1) * 2]] = matrixVelz[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[2];

                    matrix[i * n + leftWingFaces[k * 2]] = matrix[i * n + leftWingFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[k * 2 + 1]] = matrix[i * n + leftWingFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + leftWingFaces[k * 2]] = matrixVelx[i * n + leftWingFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + leftWingFaces[k * 2]] = matrixVely[i * n + leftWingFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + leftWingFaces[k * 2]] = matrixVelz[i * n + leftWingFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + leftWingFaces[k * 2 + 1]] = matrixVelx[i * n + leftWingFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + leftWingFaces[k * 2 + 1]] = matrixVely[i * n + leftWingFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + leftWingFaces[k * 2 + 1]] = matrixVelz[i * n + leftWingFaces[k * 2 + 1]] + lineVel[2];

                    indexLine1 = leftWingGrid[(k - 1) * nWakeLeftWing] * 3;
                    p1Line.x = leftWingVertices[indexLine1];
                    p1Line.y = leftWingVertices[indexLine1 + 1];
                    p1Line.z = leftWingVertices[indexLine1 + 2];

                    indexLine2 = leftWingGrid[k * nWakeLeftWing] * 3;
                    p2Line.x = leftWingVertices[indexLine2];
                    p2Line.y = leftWingVertices[indexLine2 + 1];
                    p2Line.z = leftWingVertices[indexLine2 + 2];

                    lineFunc(p, p1Line, p2Line, lineVel);

                    matrix[i * n + leftWingFaces[(k - 1) * 2]] = matrix[i * n + leftWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + leftWingFaces[(k - 1) * 2]] = matrixVelx[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + leftWingFaces[(k - 1) * 2]] = matrixVely[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + leftWingFaces[(k - 1) * 2]] = matrixVelz[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[2];
                }
            }
        }

        // Right wing
        for (k = 0; k < nSpanRightWing; k++)
        {

            for (l = 0; l < nWakeRightWing - 1; l++)
            {

                indexLine1 = rightWingGrid[k * nWakeRightWing + l] * 3;
                p1Line.x = rightWingVertices[indexLine1];
                p1Line.y = rightWingVertices[indexLine1 + 1];
                p1Line.z = rightWingVertices[indexLine1 + 2];

                indexLine2 = rightWingGrid[k * nWakeRightWing + l + 1] * 3;
                p2Line.x = rightWingVertices[indexLine2];
                p2Line.y = rightWingVertices[indexLine2 + 1];
                p2Line.z = rightWingVertices[indexLine2 + 2];

                lineFunc(p, p1Line, p2Line, lineVel);

                if (k == 0)
                {

                    matrix[i * n + rightWingFaces[k * 2]] = matrix[i * n + rightWingFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[k * 2 + 1]] = matrix[i * n + rightWingFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + rightWingFaces[k * 2]] = matrixVelx[i * n + rightWingFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + rightWingFaces[k * 2]] = matrixVely[i * n + rightWingFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + rightWingFaces[k * 2]] = matrixVelz[i * n + rightWingFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + rightWingFaces[k * 2 + 1]] = matrixVelx[i * n + rightWingFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + rightWingFaces[k * 2 + 1]] = matrixVely[i * n + rightWingFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + rightWingFaces[k * 2 + 1]] = matrixVelz[i * n + rightWingFaces[k * 2 + 1]] + lineVel[2];
                }
                else if (k == nSpanRightWing - 1)
                {

                    matrix[i * n + rightWingFaces[(k - 1) * 2]] = matrix[i * n + rightWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + rightWingFaces[(k - 1) * 2]] = matrixVelx[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + rightWingFaces[(k - 1) * 2]] = matrixVely[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + rightWingFaces[(k - 1) * 2]] = matrixVelz[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[2];
                }
                else
                {

                    matrix[i * n + rightWingFaces[(k - 1) * 2]] = matrix[i * n + rightWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + rightWingFaces[(k - 1) * 2]] = matrixVelx[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + rightWingFaces[(k - 1) * 2]] = matrixVely[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + rightWingFaces[(k - 1) * 2]] = matrixVelz[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[2];

                    matrix[i * n + rightWingFaces[k * 2]] = matrix[i * n + rightWingFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[k * 2 + 1]] = matrix[i * n + rightWingFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + rightWingFaces[k * 2]] = matrixVelx[i * n + rightWingFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + rightWingFaces[k * 2]] = matrixVely[i * n + rightWingFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + rightWingFaces[k * 2]] = matrixVelz[i * n + rightWingFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + rightWingFaces[k * 2 + 1]] = matrixVelx[i * n + rightWingFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + rightWingFaces[k * 2 + 1]] = matrixVely[i * n + rightWingFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + rightWingFaces[k * 2 + 1]] = matrixVelz[i * n + rightWingFaces[k * 2 + 1]] + lineVel[2];

                    indexLine1 = rightWingGrid[(k - 1) * nWakeRightWing] * 3;
                    p1Line.x = rightWingVertices[indexLine1];
                    p1Line.y = rightWingVertices[indexLine1 + 1];
                    p1Line.z = rightWingVertices[indexLine1 + 2];

                    indexLine2 = rightWingGrid[k * nWakeRightWing] * 3;
                    p2Line.x = rightWingVertices[indexLine2];
                    p2Line.y = rightWingVertices[indexLine2 + 1];
                    p2Line.z = rightWingVertices[indexLine2 + 2];

                    lineFunc(p, p1Line, p2Line, lineVel);

                    matrix[i * n + rightWingFaces[(k - 1) * 2]] = matrix[i * n + rightWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + rightWingFaces[(k - 1) * 2]] = matrixVelx[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + rightWingFaces[(k - 1) * 2]] = matrixVely[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + rightWingFaces[(k - 1) * 2]] = matrixVelz[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[2];
                }
            }
        }

        // Tail
        for (k = 0; k < nSpanTail; k++)
        {

            for (l = 0; l < nWakeTail - 1; l++)
            {

                indexLine1 = tailGrid[k * nWakeTail + l] * 3;
                p1Line.x = tailVertices[indexLine1];
                p1Line.y = tailVertices[indexLine1 + 1];
                p1Line.z = tailVertices[indexLine1 + 2];

                indexLine2 = tailGrid[k * nWakeTail + l + 1] * 3;
                p2Line.x = tailVertices[indexLine2];
                p2Line.y = tailVertices[indexLine2 + 1];
                p2Line.z = tailVertices[indexLine2 + 2];

                lineFunc(p, p1Line, p2Line, lineVel);

                if (k == 0)
                {

                    matrix[i * n + tailFaces[k * 2]] = matrix[i * n + tailFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[k * 2 + 1]] = matrix[i * n + tailFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + tailFaces[k * 2]] = matrixVelx[i * n + tailFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + tailFaces[k * 2]] = matrixVely[i * n + tailFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + tailFaces[k * 2]] = matrixVelz[i * n + tailFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + tailFaces[k * 2 + 1]] = matrixVelx[i * n + tailFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + tailFaces[k * 2 + 1]] = matrixVely[i * n + tailFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + tailFaces[k * 2 + 1]] = matrixVelz[i * n + tailFaces[k * 2 + 1]] + lineVel[2];
                }
                else if (k == nSpanTail - 1)
                {

                    matrix[i * n + tailFaces[(k - 1) * 2]] = matrix[i * n + tailFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[(k - 1) * 2 + 1]] = matrix[i * n + tailFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + tailFaces[(k - 1) * 2]] = matrixVelx[i * n + tailFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + tailFaces[(k - 1) * 2]] = matrixVely[i * n + tailFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + tailFaces[(k - 1) * 2]] = matrixVelz[i * n + tailFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[2];
                }
                else
                {

                    matrix[i * n + tailFaces[(k - 1) * 2]] = matrix[i * n + tailFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[(k - 1) * 2 + 1]] = matrix[i * n + tailFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + tailFaces[(k - 1) * 2]] = matrixVelx[i * n + tailFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + tailFaces[(k - 1) * 2]] = matrixVely[i * n + tailFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + tailFaces[(k - 1) * 2]] = matrixVelz[i * n + tailFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[2];

                    matrix[i * n + tailFaces[k * 2]] = matrix[i * n + tailFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[k * 2 + 1]] = matrix[i * n + tailFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + tailFaces[k * 2]] = matrixVelx[i * n + tailFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + tailFaces[k * 2]] = matrixVely[i * n + tailFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + tailFaces[k * 2]] = matrixVelz[i * n + tailFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + tailFaces[k * 2 + 1]] = matrixVelx[i * n + tailFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + tailFaces[k * 2 + 1]] = matrixVely[i * n + tailFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + tailFaces[k * 2 + 1]] = matrixVelz[i * n + tailFaces[k * 2 + 1]] + lineVel[2];

                    indexLine1 = tailGrid[(k - 1) * nWakeTail] * 3;
                    p1Line.x = tailVertices[indexLine1];
                    p1Line.y = tailVertices[indexLine1 + 1];
                    p1Line.z = tailVertices[indexLine1 + 2];

                    indexLine2 = tailGrid[k * nWakeTail] * 3;
                    p2Line.x = tailVertices[indexLine2];
                    p2Line.y = tailVertices[indexLine2 + 1];
                    p2Line.z = tailVertices[indexLine2 + 2];

                    lineFunc(p, p1Line, p2Line, lineVel);

                    matrix[i * n + tailFaces[(k - 1) * 2]] = matrix[i * n + tailFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[(k - 1) * 2 + 1]] = matrix[i * n + tailFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + tailFaces[(k - 1) * 2]] = matrixVelx[i * n + tailFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + tailFaces[(k - 1) * 2]] = matrixVely[i * n + tailFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + tailFaces[(k - 1) * 2]] = matrixVelz[i * n + tailFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[2];
                }
            }
        }
    }

    free(sourceVel);
    free(doubletVel);
    free(lineVel);
}

void calculateDoubletDistribution(int n,
                                  double *A,
                                  double *b,
                                  double *transpiration,
                                  double *sol)
{

    /* Loop parameter */
    int i;

    /* Copy the input */
    double *aux = (double *)malloc(n * n * sizeof(double));
    for (i = 0; i < n * n; i++)
        aux[i] = A[i];
    for (i = 0; i < n; i++)
        sol[i] = b[i] + transpiration[i];

    /* Solve the linear system */
    solveLinearSystem(n, aux, sol);

    free(aux);
}

void calculateSurfaceParameters(int n,
                                double *matrixVelx,
                                double *matrixVely,
                                double *matrixVelz,
                                double *arrayVel,
                                double *doublet,
                                double freestream,
                                double *velx, double *vely, double *velz,
                                double *velNorm,
                                double *cp,
                                double *mach,
                                double sound_speed)
{

    /* Loop parameter */
    int i, j, line1, line2, point;

    /* Calculate parameters */
    for (i = 0; i < n; i++)
    {

        line1 = i * 3;

        // Set velocity equal 0
        velx[i] = arrayVel[line1];
        vely[i] = arrayVel[line1 + 1];
        velz[i] = arrayVel[line1 + 2];

        // Current line
        line2 = i * n;

        for (j = 0; j < n; j++)
        {

            point = line2 + j;

            velx[i] = velx[i] + matrixVelx[point] * doublet[j];
            vely[i] = vely[i] + matrixVely[point] * doublet[j];
            velz[i] = velz[i] + matrixVelz[point] * doublet[j];
        }

        // Velocity norm
        velNorm[i] = sqrt(velx[i] * velx[i] + vely[i] * vely[i] + velz[i] * velz[i]);

        // Pressure coefficient
        cp[i] = 1 - pow(velNorm[i] / freestream, 2);

        // Mach
        mach[i] = velNorm[i] / sound_speed;
    }
}

/*
#####################################################
    BOUNDARY LAYER
#####################################################
*/
void calculateProfiles(double delta,
                       double A,
                       double B,
                       double Psi,
                       double Ctau1,
                       double Ctau2,
                       struct FreestreamParameters *freestream,
                       struct ProfileParameters *profiles)
{

    /* Parameters */
    int i;                             // loop
    int flow_type;                     // laminar or turbulent
    double delta_eta;                  // eta step or first step
    double Re_delta;                   // reynolds number
    double f0, f1, f2, f3;             // laminar curves
    double mu_mui;                     // viscosity ratio
    double h_hi;                       // enthalpy ratio
    double epsilon_line;               // thermodynamic aux parameter
    double Utau, Wtau, qtau;           // turbulent shear velocities
    double y_plus_1;                   // first element height
    double u_plus_max;                 // u_plus(delta_plus)
    double g0;                         // outer layer profile
    double dg0deta;                    // outer layer profile derivative
    double delta_plus;                 // dimensionless turbulent boundary layer height
    double Upsilon, K;                 // turbulent outer layer parameters
    double *u_plus, *y_plus;           // turbulent law of the wall profiles
    double k, C;                       // law of the wall parameters
    double u_min, y_min, u_max, y_max; // buffer region interpolation limits
    double log_y_min, log_y_max;       // log of the interpolation limits
    double a, b, c, t;                 // interpolation parameters
    double exp_ratio;                  // expansion ratio
    double eta2, eta3, eta4, eta5;     // power of eta

    /* Initilize */
    if (sqrt(pow(Ctau1, 2) + pow(Ctau2, 2)) > CTAU_CRIT)
    {
        flow_type = 1;
    }
    else
    {
        flow_type = 0;
    };
    Re_delta = freestream->velocity * freestream->density * delta / freestream->viscosity;
    u_plus = (double *)malloc(profiles->n * sizeof(double));
    y_plus = (double *)malloc(profiles->n * sizeof(double));
    k = 0.41;
    C = 5.0;
    u_min = 5.0;
    y_min = 5.0;
    u_max = 17.922725284263503;
    y_max = 200;
    log_y_min = log10(y_min);
    log_y_max = log10(y_max);
    a = u_min + 10 * log_y_min * log(10) * 0.26957378;
    b = 14.2135593;
    c = u_max - (1 / (k * log10(M_E))) * 0.51958278;
    y_plus_1 = 0.1;
    exp_ratio = 0.95;
    Utau = A / ((pow(pow(A, 2) + pow(B, 2), 0.25) * sqrt(Re_delta)) + 1e-10);
    Wtau = B / ((pow(pow(A, 2) + pow(B, 2), 0.25) * sqrt(Re_delta)) + 1e-10);
    epsilon_line = 0.2 * pow(freestream->mach, 2);

    /* Define eta */
    if (flow_type == LAMINAR_FLOW)
    {

        delta_eta = 1 / ((double)profiles->n - 1);

        for (i = 0; i < profiles->n; i++)
        {
            profiles->eta[i] = i * delta_eta;
        }
    }
    else
    {

        h_hi = 1 + epsilon_line;
        mu_mui = pow(h_hi, 1.5) * 2 / (h_hi + 1);
        delta_plus = sqrt(Re_delta) * (1 / mu_mui) * (1 / h_hi) * pow(pow(A, 2) + pow(B, 2), 0.25);

        if (delta_plus / ((double)profiles->n - 1) > y_plus_1)
        { // Geoemtric distribution

            // Find expansion ratio
            find_exp_ratio(delta_plus, y_plus_1, (double)profiles->n - 1, &exp_ratio);

            for (i = 0; i < profiles->n; i++)
            {
                if (i == 0)
                {
                    y_plus[i] = 0.0;
                }
                else if (i == profiles->n - 1)
                {
                    y_plus[i] = delta_plus;
                }
                else
                {
                    y_plus[i] = y_plus[i - 1] + y_plus_1 * pow(exp_ratio, i - 1);
                }
                profiles->eta[i] = y_plus[i] / delta_plus;
            }
        }
        else
        { // Linear distribution

            delta_eta = 1 / ((double)profiles->n - 1);

            for (i = 0; i < profiles->n; i++)
            {
                profiles->eta[i] = i * delta_eta;
                y_plus[i] = delta_plus * profiles->eta[i];
            }
        }
    }

    /* Create profiles */
    if (flow_type == LAMINAR_FLOW)
    {

        // Velocities
        for (i = 0; i < profiles->n; i++)
        {

            eta2 = pow(profiles->eta[i], 2);
            eta3 = pow(profiles->eta[i], 3);
            eta4 = pow(profiles->eta[i], 4);
            eta5 = pow(profiles->eta[i], 5);

            f0 = 6 * eta2 - 8 * eta3 + 3 * eta4;
            f1 = profiles->eta[i] - 3 * eta2 + 3 * eta3 - eta4;
            f2 = (profiles->eta[i] - 4 * eta2 + 6 * eta3 - 4 * eta4 + eta5) * pow(1 - profiles->eta[i], 2);
            f3 = (eta2 - 3 * eta3 + 3 * eta4 - eta5) * pow(1 - profiles->eta[i], 2);

            profiles->U[i] = A * (1 - 0.6 * (A - 3) * eta3) * f1 + f0;
            profiles->W[i] = B * f2 + Psi * f3;

            profiles->R[i] = 1 / (1 + epsilon_line * (1 - pow(profiles->U[i], 2) - pow(profiles->W[i], 2)));
        }

        // Gradient and Shear stress
        for (i = 0; i < profiles->n; i++)
        {

            f0 = 12 * profiles->eta[i] - 24 * pow(profiles->eta[i], 2) + 12 * pow(profiles->eta[i], 3);
            f1 = 1 - 6 * profiles->eta[i] + 9 * pow(profiles->eta[i], 2) - 4 * pow(profiles->eta[i], 3);
            f2 = (1 - 8 * profiles->eta[i] + 18 * pow(profiles->eta[i], 2) - 16 * pow(profiles->eta[i], 3) + 5 * pow(profiles->eta[i], 4)) * pow(1 - profiles->eta[i], 2) - 2 * (1 - profiles->eta[i]) * (profiles->eta[i] - 3 * pow(profiles->eta[i], 2) + 3 * pow(profiles->eta[i], 3) - pow(profiles->eta[i], 4));
            f3 = (2 * profiles->eta[i] - 9 * pow(profiles->eta[i], 2) + 12 * pow(profiles->eta[i], 3) - 5 * pow(profiles->eta[i], 4)) * pow(1 - profiles->eta[i], 2) - 2 * (1 - profiles->eta[i]) * (pow(profiles->eta[i], 2) - 3 * pow(profiles->eta[i], 3) + 3 * pow(profiles->eta[i], 4) - pow(profiles->eta[i], 5));

            profiles->dU_deta[i] = -1.8 * A * (A - 3) * pow(profiles->eta[i], 2) * (profiles->eta[i] - 3 * pow(profiles->eta[i], 2) + 3 * pow(profiles->eta[i], 3) - pow(profiles->eta[i], 4)) + A * (1 - 0.6 * (A - 3) * pow(profiles->eta[i], 3)) * f1 + f0;
            profiles->dW_deta[i] = B * f2 + Psi * f3;

            mu_mui = (1 / profiles->R[i], 1.5) * 2 / (1 / profiles->R[i] + 1);

            profiles->S[i] = (mu_mui / Re_delta) * profiles->dU_deta[i];
            profiles->T[i] = (mu_mui / Re_delta) * profiles->dW_deta[i];
        }
    }
    else
    {

        // u_plus_max
        if (delta_plus <= y_min)
        {
            u_plus_max = delta_plus;
        }
        else if ((y_min < delta_plus) && (delta_plus < y_max))
        {
            t = (log10(delta_plus) - log_y_min) / (log_y_max - log_y_min);
            u_plus_max = pow(1 - t, 4) * u_min + 4 * pow(1 - t, 3) * t * a + 6 * pow(1 - t, 2) * pow(t, 2) * b + 4 * (1 - t) * pow(t, 3) * c + pow(t, 4) * u_max;
        }
        else
        {
            u_plus_max = (1 / k) * log(delta_plus) + C;
        }

        K = sqrt(pow(Wtau * u_plus_max, 2) + pow(1 - Utau * u_plus_max, 2));
        Upsilon = atan(Wtau * u_plus_max / (1 - Utau * u_plus_max));

        // Velocities
        for (i = 0; i < profiles->n; i++)
        {

            // u_plus
            if (y_plus[i] < y_min)
            {
                u_plus[i] = y_plus[i];
            }
            else if ((y_min <= y_plus[i]) && (y_plus[i] <= y_max))
            {
                t = (log10(y_plus[i]) - log_y_min) / (log_y_max - log_y_min);
                u_plus[i] = pow(1 - t, 4) * u_min + 4 * pow(1 - t, 3) * t * a + 6 * pow(1 - t, 2) * pow(t, 2) * b + 4 * (1 - t) * pow(t, 3) * c + pow(t, 4) * u_max;
            }
            else
            {
                u_plus[i] = (1 / k) * log(y_plus[i]) + C;
            }

            g0 = 3 * pow(profiles->eta[i], 2) - 2 * pow(profiles->eta[i], 3);

            profiles->U[i] = Utau * u_plus[i] + K * cos(Upsilon - Psi * pow(1 - profiles->eta[i], 2)) * g0;
            profiles->W[i] = Wtau * u_plus[i] - K * sin(Upsilon - Psi * pow(1 - profiles->eta[i], 2)) * g0;

            profiles->R[i] = 1 / (1 + epsilon_line * (1 - pow(profiles->U[i], 2) - pow(profiles->W[i], 2)));
        }

        // Gradient and Shear stress
        for (i = 0; i < profiles->n; i++)
        {

            g0 = 3 * pow(profiles->eta[i], 2) - 2 * pow(profiles->eta[i], 3);

            profiles->dU_deta[i] = Utau * delta_plus * (1 / (1 + exp(-k * C) * (k * exp(k * u_plus[i]) - k - pow(k, 2) * u_plus[i] - 0.5 * k * pow(k * u_plus[i], 2)))) + 2 * Psi * (1 - profiles->eta[i]) * K * sin(Upsilon - Psi * pow(1 - profiles->eta[i], 2)) * g0 + K * cos(Upsilon - Psi * pow(1 - profiles->eta[i], 2)) * (6 * profiles->eta[i] - 6 * pow(profiles->eta[i], 2));
            profiles->dW_deta[i] = Wtau * delta_plus * (1 / (1 + exp(-k * C) * (k * exp(k * u_plus[i]) - k - pow(k, 2) * u_plus[i] - 0.5 * k * pow(k * u_plus[i], 2)))) + 2 * Psi * (1 - profiles->eta[i]) * K * cos(Upsilon - Psi * pow(1 - profiles->eta[i], 2)) * g0 - K * sin(Upsilon - Psi * pow(1 - profiles->eta[i], 2)) * (6 * profiles->eta[i] - 6 * pow(profiles->eta[i], 2));

            mu_mui = (1 / profiles->R[i], 1.5) * 2 / (1 / profiles->R[i] + 1);

            g0 = 3 * pow(profiles->eta[i], 2) - 2 * pow(profiles->eta[i], 3);
            dg0deta = 6 * profiles->eta[i] - 6 * pow(profiles->eta[i], 2);

            profiles->S[i] = profiles->R[i] * Utau * sqrt(pow(Utau, 2) + pow(Wtau, 2)) * (1 - g0) + profiles->R[i] * Ctau1 * K * cos(Upsilon - Psi * pow(1 - profiles->eta[i], 2)) * dg0deta;
            profiles->T[i] = profiles->R[i] * Wtau * sqrt(pow(Utau, 2) + pow(Wtau, 2)) * (1 - g0) + profiles->R[i] * Ctau2 * K * sin(Upsilon - Psi * pow(1 - profiles->eta[i], 2)) * dg0deta;
        }
    }

    /* Free arrays */
    free(u_plus);
    free(y_plus);
}

void calculateIntegralThickness(struct ProfileParameters *profiles,
                                struct IntegralThicknessParameters *integralThickness,
                                double delta,
                                double Psi)
{

    /* Parameters */
    int i;
    double *func;

    /* Initialize */
    func = (double *)malloc(profiles->n * sizeof(double));

    /* Calculate integral */
    for (i = 0; i < profiles->n; i++)
        func[i] = 1 - profiles->R[i] * profiles->U[i];
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->delta_1_ast, delta);

    for (i = 0; i < profiles->n; i++)
        func[i] = -profiles->R[i] * profiles->W[i];
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->delta_2_ast, delta);

    for (i = 0; i < profiles->n; i++)
        func[i] = 1 - profiles->R[i] * pow(profiles->U[i], 2);
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->phi_11, delta);

    for (i = 0; i < profiles->n; i++)
        func[i] = -profiles->R[i] * profiles->U[i] * profiles->W[i];
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->phi_12, delta);
    integralThickness->phi_21 = integralThickness->phi_12;

    for (i = 0; i < profiles->n; i++)
        func[i] = -profiles->R[i] * pow(profiles->W[i], 2);
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->phi_22, delta);

    for (i = 0; i < profiles->n; i++)
        func[i] = 1 - profiles->R[i] * profiles->U[i] * (pow(profiles->U[i], 2) + pow(profiles->W[i], 2));
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->phi_1_ast, delta);

    for (i = 0; i < profiles->n; i++)
        func[i] = -profiles->R[i] * profiles->W[i] * (pow(profiles->U[i], 2) + pow(profiles->W[i], 2));
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->phi_2_ast, delta);

    for (i = 0; i < profiles->n; i++)
        func[i] = 1 - profiles->U[i];
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->delta_1_line, delta);

    for (i = 0; i < profiles->n; i++)
        func[i] = -profiles->W[i];
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->delta_2_line, delta);

    integralThickness->delta_q = integralThickness->phi_11 + integralThickness->phi_22;

    for (i = 0; i < profiles->n; i++)
        func[i] = -Psi * profiles->R[i] * (pow(profiles->U[i], 2) + pow(profiles->W[i], 2));
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->delta_q_o, delta);

    for (i = 0; i < profiles->n; i++)
        func[i] = -Psi * profiles->R[i] * profiles->U[i] * (pow(profiles->U[i], 2) + pow(profiles->W[i], 2));
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->theta_1_o, delta);

    for (i = 0; i < profiles->n; i++)
        func[i] = -Psi * profiles->R[i] * profiles->W[i] * (pow(profiles->U[i], 2) + pow(profiles->W[i], 2));
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->theta_2_o, delta);

    for (i = 0; i < profiles->n; i++)
        func[i] = -Psi * profiles->U[i];
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->delta_1_o, delta);

    for (i = 0; i < profiles->n; i++)
        func[i] = -Psi * profiles->W[i];
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->delta_2_o, delta);

    for (i = 0; i < profiles->n; i++)
        func[i] = profiles->S[i] * profiles->dU_deta[i] + profiles->T[i] * profiles->dW_deta[i];
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->C_D, 1.0);

    for (i = 0; i < profiles->n; i++)
        func[i] = profiles->S[i] * profiles->dW_deta[i] - profiles->T[i] * profiles->dU_deta[i];
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->C_D_x, 1.0);

    for (i = 0; i < profiles->n; i++)
        func[i] = Psi * (profiles->S[i] * profiles->dW_deta[i] - profiles->T[i] * profiles->dU_deta[i]);
    integrate_trap(profiles->n, profiles->eta, func, &integralThickness->C_D_o, 1.0);

    integralThickness->C_f_1 = 2 * profiles->S[0];
    integralThickness->C_f_2 = 2 * profiles->T[0];

    integralThickness->theta_11 = integralThickness->phi_11 - integralThickness->delta_1_line;
    integralThickness->theta_22 = integralThickness->phi_22 - integralThickness->delta_2_line;

    /* Free arrays */
    free(func);
}

void calculateIntegralDefect(struct ProfileParameters *profiles,
                             struct IntegralThicknessParameters *integralThickness,
                             struct FreestreamParameters *freestream,
                             struct IntegralDefectParameters *integralDefect,
                             double delta,
                             double A, double B,
                             double Ctau1, double Ctau2) {

    /* Aux */
    double aux_1, aux_2, aux_3;

    /* Initialize */
    aux_1 = freestream->density * freestream->velocity;
    aux_2 = freestream->density * freestream->velocity * freestream->velocity;
    aux_3 = freestream->density * freestream->velocity * freestream->velocity * freestream->velocity;

    /* Calculate defect parameters */
    integralDefect->M_x = aux_1 * integralThickness->delta_1_ast;
    integralDefect->M_y = aux_1 * integralThickness->delta_2_ast;

    integralDefect->J_xx = aux_2 * integralThickness->phi_11;
    integralDefect->J_xy = aux_2 * integralThickness->phi_12;
    integralDefect->J_yx = aux_2 * integralThickness->phi_21;
    integralDefect->J_yy = aux_2 * integralThickness->phi_22;

    integralDefect->E_x = aux_3 * integralThickness->phi_1_ast;
    integralDefect->E_y = aux_3 * integralThickness->phi_2_ast;

    integralDefect->K_o_x = aux_3 * integralThickness->theta_1_o;
    integralDefect->K_o_y = aux_3 * integralThickness->theta_2_o;

    integralDefect->Q_x = freestream->velocity * integralThickness->delta_1_line;
    integralDefect->Q_y = freestream->velocity * integralThickness->delta_2_line;

    integralDefect->Q_o_x = freestream->velocity * integralThickness->theta_1_o;
    integralDefect->Q_o_y = freestream->velocity * integralThickness->theta_2_o;

    integralDefect->tau_w_x = 0.5 * aux_2 * integralThickness->C_f_1;
    integralDefect->tau_w_y = 0.5 * aux_2 * integralThickness->C_f_2;

    integralDefect->D = aux_3 * integralThickness->C_D;
    integralDefect->D_x = aux_3 * integralThickness->C_D_x;
    integralDefect->D_o = aux_3 * integralThickness->C_D_o;

    double mod_Ctau = sqrt(pow(Ctau1, 2) + pow(Ctau2, 2));

    if (mod_Ctau <= CTAU_CRIT) {

        double H_k_1, H_k_2, Re_theta_11, Re_theta_22, f1, f2;

        if (absValue(integralThickness->delta_1_ast / (integralThickness->theta_11 + 1e-8)) > 10 || absValue(integralThickness->delta_1_ast / (integralThickness->theta_11 + 1e-8)) <= 1.5) {
            integralDefect->S_tau_x = 0.0;
        } else {
            H_k_1 = (absValue(integralThickness->delta_1_ast / (integralThickness->theta_11 + 1e-8)) - 0.29 * freestream->mach * freestream->mach) / (1 + 0.113 * freestream->mach * freestream->mach);
            Re_theta_11 = freestream->velocity * freestream->density * absValue(integralThickness->theta_11) / freestream->viscosity;
            f1 = 0.01 * sqrt(pow(2.4 * H_k_1 - 3.7 + 2.5 * tanh(1.5 * H_k_1 - 4.65), 2) + 0.25) * (Re_theta_11 - pow(10, (1.415 / (H_k_1 - 1) - 0.489) * tanh(20 / (H_k_1 - 1) - 12.9)));
            integralDefect->S_tau_x = f1 * freestream->velocity * mod_Ctau / integralThickness->theta_11;
        }

        if (absValue(integralThickness->delta_2_ast / (integralThickness->theta_22 + 1e-8)) > 10 || absValue(integralThickness->delta_2_ast / (integralThickness->theta_22 + 1e-8)) <= 1.5) {
            integralDefect->S_tau_y = 0.0;
        } else {
            H_k_2 = (absValue(integralThickness->delta_2_ast / (integralThickness->theta_22 + 1e-8)) - 0.29 * freestream->mach * freestream->mach) / (1 + 0.113 * freestream->mach * freestream->mach);
            Re_theta_22 = freestream->velocity * freestream->density * absValue(integralThickness->theta_22) / freestream->viscosity;
            f2 = 0.01 * sqrt(pow(2.4 * H_k_2 - 3.7 + 2.5 * tanh(1.5 * H_k_2 - 4.65), 2) + 0.25) * (Re_theta_22 - pow(10, (1.415 / (H_k_2 - 1) - 0.489) * tanh(20 / (H_k_2 - 1) - 12.9)));
            integralDefect->S_tau_y = f2 * freestream->velocity * mod_Ctau / integralThickness->theta_11;
        }

    } else {

        double P_tau_x;
        double P_tau_y;
        double D_tau_x;
        double D_tau_y;

        double *P_tau_x_func = (double *)malloc(profiles->n * sizeof(double));
        double *P_tau_y_func = (double *)malloc(profiles->n * sizeof(double));
        double *D_tau_x_func = (double *)malloc(profiles->n * sizeof(double));
        double *D_tau_y_func = (double *)malloc(profiles->n * sizeof(double));

        double tau_x, tau_y;

        for (int i = 0; i < profiles->n; i++) {

            P_tau_x_func[i] = freestream->density * profiles->R[i] * (freestream->velocity * freestream->velocity / profiles->R[i]) * sqrt(pow(profiles->S[i], 2) + pow(profiles->T[i], 2)) * freestream->velocity * profiles->dU_deta[i];
            P_tau_y_func[i] = freestream->density * profiles->R[i] * (freestream->velocity * freestream->velocity / profiles->R[i]) * sqrt(pow(profiles->S[i], 2) + pow(profiles->T[i], 2)) * freestream->velocity * profiles->dW_deta[i];

            tau_x = freestream->velocity * freestream->velocity * profiles->S[i] / profiles->R[i];
            tau_y = freestream->velocity * freestream->velocity * profiles->T[i] / profiles->R[i];

            D_tau_x_func[i] = 2 * freestream->density * profiles->R[i] * pow(pow(tau_x, 2) + pow(tau_y, 2), 0.25) * tau_x;
            D_tau_y_func[i] = 2 * freestream->density * profiles->R[i] * pow(pow(tau_x, 2) + pow(tau_y, 2), 0.25) * tau_y;
        }

        integrate_trap(profiles->n, profiles->eta, P_tau_x_func, &P_tau_x, 1.0);
        integrate_trap(profiles->n, profiles->eta, P_tau_y_func, &P_tau_y, 1.0);
        integrate_trap(profiles->n, profiles->eta, D_tau_x_func, &D_tau_x, delta);
        integrate_trap(profiles->n, profiles->eta, D_tau_y_func, &D_tau_y, delta);

        integralDefect->S_tau_x = 0.30 * (P_tau_x - D_tau_x);
        integralDefect->S_tau_y = 0.30 * (P_tau_y - D_tau_y);

        free(P_tau_x_func);
        free(P_tau_y_func);
        free(D_tau_x_func);
        free(D_tau_y_func);

    }

    double *K_tau_xx_func = (double *)malloc(profiles->n * sizeof(double));
    double *K_tau_xy_func = (double *)malloc(profiles->n * sizeof(double));
    double *K_tau_yx_func = (double *)malloc(profiles->n * sizeof(double));
    double *K_tau_yy_func = (double *)malloc(profiles->n * sizeof(double));

    double tau_x, tau_y;

    for (int i = 0; i < profiles->n; i++) {

        tau_x = freestream->velocity * freestream->velocity * profiles->S[i] / profiles->R[i];
        tau_y = freestream->velocity * freestream->velocity * profiles->T[i] / profiles->R[i];

        K_tau_xx_func[i] = profiles->R[i] * freestream->density * tau_x * freestream->velocity * profiles->U[i];
        K_tau_xy_func[i] = profiles->R[i] * freestream->density * tau_x * freestream->velocity * profiles->W[i];
        K_tau_yx_func[i] = profiles->R[i] * freestream->density * tau_y * freestream->velocity * profiles->U[i];
        K_tau_yy_func[i] = profiles->R[i] * freestream->density * tau_y * freestream->velocity * profiles->W[i];
    }

    integrate_trap(profiles->n, profiles->eta, K_tau_xx_func, &integralDefect->K_tau_xx, delta);
    integrate_trap(profiles->n, profiles->eta, K_tau_xy_func, &integralDefect->K_tau_xy, delta);
    integrate_trap(profiles->n, profiles->eta, K_tau_yx_func, &integralDefect->K_tau_yx, delta);
    integrate_trap(profiles->n, profiles->eta, K_tau_yy_func, &integralDefect->K_tau_yy, delta);

    free(K_tau_xx_func);
    free(K_tau_xy_func);
    free(K_tau_yx_func);
    free(K_tau_yy_func);
}

void calculateEquationsParams(double delta,
                              double A,
                              double B,
                              double Psi,
                              double Ctau1,
                              double Ctau2,
                              struct FreestreamParameters *freestream,
                              struct ProfileParameters *profiles,
                              struct IntegralThicknessParameters *integralThickness,
                              struct IntegralDefectParameters *integralDefect,
                              struct EquationsParameters *params)
{

    /* Profiles */
    calculateProfiles(delta, A, B, Psi, Ctau1, Ctau2, freestream, profiles);

    /* Integral thickness */
    calculateIntegralThickness(profiles, integralThickness, delta, Psi);

    /* Integral defect */
    calculateIntegralDefect(profiles, integralThickness, freestream, integralDefect, delta, A, B, Ctau1, Ctau2);

    /* Assing params */
    params->D = integralDefect->D;
    params->D_o = integralDefect->D_o;
    params->D_x = integralDefect->D_x;
    params->E_x = integralDefect->E_x;
    params->E_y = integralDefect->E_y;
    params->J_xx = integralDefect->J_xx;
    params->J_xy = integralDefect->J_xy;
    params->J_yx = integralDefect->J_yx;
    params->J_yy = integralDefect->J_yy;
    params->K_o_x = integralDefect->K_o_x;
    params->K_o_y = integralDefect->K_o_y;
    params->K_tau_xx = integralDefect->K_tau_xx;
    params->K_tau_xy = integralDefect->K_tau_xy;
    params->K_tau_yx = integralDefect->K_tau_yx;
    params->K_tau_yy = integralDefect->K_tau_yy;
    params->M_x = integralDefect->M_x;
    params->M_y = integralDefect->M_y;
    params->Q_o_x = integralDefect->Q_o_x;
    params->Q_o_y = integralDefect->Q_o_y;
    params->Q_x = integralDefect->Q_x;
    params->Q_y = integralDefect->Q_y;
    params->S_tau_x = integralDefect->S_tau_x;
    params->S_tau_y = integralDefect->S_tau_y;
    params->tau_w_x = integralDefect->tau_w_x;
    params->tau_w_y = integralDefect->tau_w_y;
    params->vel = freestream->velocity;
    params->density = freestream->density;
}

void calculateDivergents(int face,
                         int *faces,
                         struct VerticeConnection *vertices_connection,
                         struct EquationsParameters *params,
                         double area,
                         double *p1, double *p2, double *p3) {

    /* Parameters */
    int k; // Counters

    struct Point point_1, point_2, point_3; // Face corners
    struct Point v_1, v_2, v_3;             // Face corners values
    struct Point point_lateral;             // Face edge orthogonal vector
    int index_1, index_2, index_3;          // Vertices indexes

    /* Initialize */
    point_1.x = p1[face];
    point_1.y = p1[face + 1];
    point_1.z = 0.0;
    point_2.x = p2[face];
    point_2.y = p2[face + 1];
    point_2.z = 0.0;
    point_3.x = p3[face];
    point_3.y = p3[face + 1];
    point_3.z = 0.0;

    v_1.z = 0.0;
    v_2.z = 0.0;
    v_3.z = 0.0;

    index_1 = faces[3 * face];
    index_2 = faces[3 * face + 1];
    index_3 = faces[3 * face + 2];

    /* Div(M) */
    v_1.x = 0.0;
    v_1.y = 0.0;
    v_2.x = 0.0;
    v_2.y = 0.0;
    v_3.x = 0.0;
    v_3.y = 0.0;

    for (k = 0; k < vertices_connection[index_1].n; k++)
    {
        v_1.x = v_1.x + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].M_x;
        v_1.y = v_1.y + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].M_y;
    }

    for (k = 0; k < vertices_connection[index_2].n; k++)
    {
        v_2.x = v_2.x + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].M_x;
        v_2.y = v_2.y + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].M_y;
    }

    for (k = 0; k < vertices_connection[index_3].n; k++)
    {
        v_3.x = v_3.x + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].M_x;
        v_3.y = v_3.y + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].M_y;
    }

    calculateDivergence(point_1, point_2, point_3, v_1, v_2, v_3, area, &params[face].div_M);

    /* Div(Jx) */
    v_1.x = 0.0;
    v_1.y = 0.0;
    v_2.x = 0.0;
    v_2.y = 0.0;
    v_3.x = 0.0;
    v_3.y = 0.0;

    for (k = 0; k < vertices_connection[index_1].n; k++)
    {
        v_1.x = v_1.x + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].J_xx;
        v_1.y = v_1.y + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].J_xy;
    }

    for (k = 0; k < vertices_connection[index_2].n; k++)
    {
        v_2.x = v_2.x + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].J_xx;
        v_2.y = v_2.y + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].J_xy;
    }

    for (k = 0; k < vertices_connection[index_3].n; k++)
    {
        v_3.x = v_3.x + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].J_xx;
        v_3.y = v_3.y + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].J_xy;
    }

    calculateDivergence(point_1, point_2, point_3, v_1, v_2, v_3, area, &params[face].div_J_x);

    /* Div(Jy) */
    v_1.x = 0.0;
    v_1.y = 0.0;
    v_2.x = 0.0;
    v_2.y = 0.0;
    v_3.x = 0.0;
    v_3.y = 0.0;

    for (k = 0; k < vertices_connection[index_1].n; k++)
    {
        v_1.x = v_1.x + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].J_yx;
        v_1.y = v_1.y + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].J_yy;
    }

    for (k = 0; k < vertices_connection[index_2].n; k++)
    {
        v_2.x = v_2.x + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].J_yx;
        v_2.y = v_2.y + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].J_yy;
    }

    for (k = 0; k < vertices_connection[index_3].n; k++)
    {
        v_3.x = v_3.x + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].J_yx;
        v_3.y = v_3.y + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].J_yy;
    }

    calculateDivergence(point_1, point_2, point_3, v_1, v_2, v_3, area, &params[face].div_J_y);

    /* Div(E) */
    v_1.x = 0.0;
    v_1.y = 0.0;
    v_2.x = 0.0;
    v_2.y = 0.0;
    v_3.x = 0.0;
    v_3.y = 0.0;

    for (k = 0; k < vertices_connection[index_1].n; k++)
    {
        v_1.x = v_1.x + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].E_x;
        v_1.y = v_1.y + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].E_y;
    }

    for (k = 0; k < vertices_connection[index_2].n; k++)
    {
        v_2.x = v_2.x + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].E_x;
        v_2.y = v_2.y + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].E_y;
    }

    for (k = 0; k < vertices_connection[index_3].n; k++)
    {
        v_3.x = v_3.x + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].E_x;
        v_3.y = v_3.y + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].E_y;
    }

    calculateDivergence(point_1, point_2, point_3, v_1, v_2, v_3, area, &params[face].div_E);

    /* Div(ko) */
    v_1.x = 0.0;
    v_1.y = 0.0;
    v_2.x = 0.0;
    v_2.y = 0.0;
    v_3.x = 0.0;
    v_3.y = 0.0;

    for (k = 0; k < vertices_connection[index_1].n; k++)
    {
        v_1.x = v_1.x + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].K_o_x;
        v_1.y = v_1.y + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].K_o_y;
    }

    for (k = 0; k < vertices_connection[index_2].n; k++)
    {
        v_2.x = v_2.x + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].K_o_x;
        v_2.y = v_2.y + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].K_o_y;
    }

    for (k = 0; k < vertices_connection[index_3].n; k++)
    {
        v_3.x = v_3.x + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].K_o_x;
        v_3.y = v_3.y + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].K_o_y;
    }

    calculateDivergence(point_1, point_2, point_3, v_1, v_2, v_3, area, &params[face].div_K_o);

    /* Div(ktaux) */
    v_1.x = 0.0;
    v_1.y = 0.0;
    v_2.x = 0.0;
    v_2.y = 0.0;
    v_3.x = 0.0;
    v_3.y = 0.0;

    for (k = 0; k < vertices_connection[index_1].n; k++)
    {
        v_1.x = v_1.x + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].K_tau_xx;
        v_1.y = v_1.y + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].K_tau_xy;
    }

    for (k = 0; k < vertices_connection[index_2].n; k++)
    {
        v_2.x = v_2.x + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].K_tau_xx;
        v_2.y = v_2.y + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].K_tau_xy;
    }

    for (k = 0; k < vertices_connection[index_3].n; k++)
    {
        v_3.x = v_3.x + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].K_tau_xx;
        v_3.y = v_3.y + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].K_tau_xy;
    }

    calculateDivergence(point_1, point_2, point_3, v_1, v_2, v_3, area, &params[face].div_K_tau_x);

    /* Div(ktauy) */
    v_1.x = 0.0;
    v_1.y = 0.0;
    v_2.x = 0.0;
    v_2.y = 0.0;
    v_3.x = 0.0;
    v_3.y = 0.0;

    for (k = 0; k < vertices_connection[index_1].n; k++)
    {
        v_1.x = v_1.x + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].K_tau_yx;
        v_1.y = v_1.y + vertices_connection[index_1].coeffs[k] * params[vertices_connection[index_1].faces[k]].K_tau_yy;
    }

    for (k = 0; k < vertices_connection[index_2].n; k++)
    {
        v_2.x = v_2.x + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].K_tau_yx;
        v_2.y = v_2.y + vertices_connection[index_2].coeffs[k] * params[vertices_connection[index_2].faces[k]].K_tau_yy;
    }

    for (k = 0; k < vertices_connection[index_3].n; k++)
    {
        v_3.x = v_3.x + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].K_tau_yx;
        v_3.y = v_3.y + vertices_connection[index_3].coeffs[k] * params[vertices_connection[index_3].faces[k]].K_tau_yy;
    }

    calculateDivergence(point_1, point_2, point_3, v_1, v_2, v_3, area, &params[face].div_K_tau_y);
}

void calculateGradients(int face,
                        int *faces,
                        struct VerticeConnection *vertices_connection,
                        struct EquationsParameters *params,
                        double *e1, double *e2, double *e3,
                        double *p1, double *p2, double *p3,
                        double *velNorm,
                        double *velx, double *vely, double *velz,
                        double *transpiration)
{

    /* Parameters */
    int k; // Counters

    struct Point point_1, point_2, point_3;    // Face corners
    struct Point v_1, v_2, v_3;                // Face corners values
    struct Point e1_point, e2_point, e3_point; // Face base vectors
    struct Point vel_point;                    // Face velocity
    struct Point point_lateral;                // Face edge orthogonal vector
    int index_1, index_2, index_3;             // Vertices indexes
    double aux_1, aux_2;

    /* Initialize */
    e1_point.x = e1[3 * face];
    e1_point.y = e1[3 * face + 1];
    e1_point.z = e1[3 * face + 2];
    e2_point.x = e2[3 * face];
    e2_point.y = e2[3 * face + 1];
    e2_point.z = e2[3 * face + 2];
    e3_point.x = e3[3 * face];
    e3_point.y = e3[3 * face + 1];
    e3_point.z = e3[3 * face + 2];

    point_1.x = p1[2 * face];
    point_1.y = p1[2 * face + 1];
    point_1.z = 0.0;
    point_2.x = p2[2 * face];
    point_2.y = p2[2 * face + 1];
    point_2.z = 0.0;
    point_3.x = p3[2 * face];
    point_3.y = p3[2 * face + 1];
    point_3.z = 0.0;

    v_1.z = 0.0;
    v_2.z = 0.0;
    v_3.z = 0.0;

    index_1 = faces[3 * face];
    index_2 = faces[3 * face + 1];
    index_3 = faces[3 * face + 2];

    vel_point.x = velx[face];
    vel_point.y = vely[face];
    vel_point.z = velz[face];

    // pow(velNorm, 2)
    point_1.z = 0.0;
    point_2.z = 0.0;
    point_3.z = 0.0;

    for (k = 0; k < vertices_connection[index_1].n; k++)
    {
        point_1.z = point_1.z + vertices_connection[index_1].coeffs[k] * velNorm[vertices_connection[index_1].faces[k]] * velNorm[vertices_connection[index_1].faces[k]];
    }

    for (k = 0; k < vertices_connection[index_2].n; k++)
    {
        point_2.z = point_2.z + vertices_connection[index_2].coeffs[k] * velNorm[vertices_connection[index_2].faces[k]] * velNorm[vertices_connection[index_2].faces[k]];
    }

    for (k = 0; k < vertices_connection[index_3].n; k++)
    {
        point_3.z = point_3.z + vertices_connection[index_3].coeffs[k] * velNorm[vertices_connection[index_3].faces[k]] * velNorm[vertices_connection[index_3].faces[k]];
    }

    calculateGradient(velNorm[face] * velNorm[face], point_1, point_2, point_3, e1_point, e2_point, e3_point, vel_point, transpiration[face], &params[face].grad_q2_x, &params[face].grad_q2_y);

    // Phi
    v_1.x = 0.0;
    v_1.y = 0.0;
    v_1.z = 0.0;
    v_2.x = 0.0;
    v_2.y = 0.0;
    v_2.z = 0.0;
    v_3.x = 0.0;
    v_3.y = 0.0;
    v_3.z = 0.0;

    for (k = 0; k < vertices_connection[index_1].n; k++)
    {
        v_1.x = v_1.x + vertices_connection[index_1].coeffs[k] * velx[vertices_connection[index_1].faces[k]];
        v_1.y = v_1.y + vertices_connection[index_1].coeffs[k] * vely[vertices_connection[index_1].faces[k]];
        v_1.z = v_1.z + vertices_connection[index_1].coeffs[k] * velz[vertices_connection[index_1].faces[k]];
    }

    for (k = 0; k < vertices_connection[index_2].n; k++)
    {
        v_2.x = v_2.x + vertices_connection[index_2].coeffs[k] * velx[vertices_connection[index_2].faces[k]];
        v_2.y = v_2.y + vertices_connection[index_2].coeffs[k] * vely[vertices_connection[index_2].faces[k]];
        v_2.z = v_2.z + vertices_connection[index_2].coeffs[k] * velz[vertices_connection[index_2].faces[k]];
    }

    for (k = 0; k < vertices_connection[index_3].n; k++)
    {
        v_3.x = v_3.x + vertices_connection[index_3].coeffs[k] * velx[vertices_connection[index_3].faces[k]];
        v_3.y = v_3.y + vertices_connection[index_3].coeffs[k] * vely[vertices_connection[index_3].faces[k]];
        v_3.z = v_3.z + vertices_connection[index_3].coeffs[k] * velz[vertices_connection[index_3].faces[k]];
    }

    // Remove vel in e3 direction
    aux_1 = e3_point.x * vel_point.x + e3_point.y * vel_point.y + e3_point.z * vel_point.z;
    vel_point.x = vel_point.x - e3_point.x * aux_1;
    vel_point.y = vel_point.y - e3_point.y * aux_1;
    vel_point.z = vel_point.z - e3_point.z * aux_1;

    aux_1 = e3_point.x * v_1.x + e3_point.y * v_1.y + e3_point.z * v_1.z;
    v_1.x = v_1.x - e3_point.x * aux_1;
    v_1.y = v_1.y - e3_point.y * aux_1;
    v_1.z = v_1.z - e3_point.z * aux_1;

    aux_1 = e3_point.x * v_2.x + e3_point.y * v_2.y + e3_point.z * v_2.z;
    v_2.x = v_2.x - e3_point.x * aux_1;
    v_2.y = v_2.y - e3_point.y * aux_1;
    v_2.z = v_2.z - e3_point.z * aux_1;

    aux_1 = e3_point.x * v_3.x + e3_point.y * v_3.y + e3_point.z * v_3.z;
    v_3.x = v_3.x - e3_point.x * aux_1;
    v_3.y = v_3.y - e3_point.y * aux_1;
    v_3.z = v_3.z - e3_point.z * aux_1;

    // Unary vectors
    aux_1 = norm(vel_point);
    vel_point.x = vel_point.x / aux_1;
    vel_point.y = vel_point.y / aux_1;
    vel_point.z = vel_point.z / aux_1;
    aux_1 = norm(v_1);
    v_1.x = v_1.x / aux_1;
    v_1.y = v_1.y / aux_1;
    v_1.z = v_1.z / aux_1;
    aux_1 = norm(v_2);
    v_2.x = v_2.x / aux_1;
    v_2.y = v_2.y / aux_1;
    v_2.z = v_2.z / aux_1;
    aux_1 = norm(v_3);
    v_3.x = v_3.x / aux_1;
    v_3.y = v_3.y / aux_1;
    v_3.z = v_3.z / aux_1;

    point_lateral = cross(vel_point, e3_point);

    aux_1 = point_lateral.x * v_1.x + point_lateral.y * v_1.y + point_lateral.z * v_1.z;
    aux_2 = vel_point.x * v_1.x + vel_point.y * v_1.y + vel_point.z * v_1.z;
    point_1.z = atan(aux_1 / aux_2);

    aux_1 = point_lateral.x * v_2.x + point_lateral.y * v_2.y + point_lateral.z * v_2.z;
    aux_2 = vel_point.x * v_2.x + vel_point.y * v_2.y + vel_point.z * v_2.z;
    point_2.z = atan(aux_1 / aux_2);

    aux_1 = point_lateral.x * v_3.x + point_lateral.y * v_3.y + point_lateral.z * v_3.z;
    aux_2 = vel_point.x * v_3.x + vel_point.y * v_3.y + vel_point.z * v_3.z;
    point_3.z = atan(aux_1 / aux_2);

    vel_point.x = velx[face];
    vel_point.y = vely[face];
    vel_point.z = velz[face];

    calculateGradient(0.0, point_1, point_2, point_3, e1_point, e2_point, e3_point, vel_point, transpiration[face], &params[face].grad_phi_x, &params[face].grad_phi_y);
}

void calculateObjectiveFunction(struct EquationsParameters params,
                                double *obj,
                                double velocity,
                                double *momentum_x,
                                double *momentum_y,
                                double *kinetic_energy,
                                double *lateral_curvature,
                                double *shear_stress_x,
                                double *shear_stress_y) {

    /* Face integral equations */
    *momentum_x = params.div_J_x - params.vel * params.div_M - params.tau_w_x;
    *momentum_y = params.div_J_y - params.tau_w_y;
    *kinetic_energy = params.div_E - params.vel * params.vel * params.div_M - params.density * (params.Q_x * params.grad_q2_x + params.Q_y * params.grad_q2_y) - 2 * params.D;
    *lateral_curvature = params.div_K_o + (params.E_x * params.grad_phi_x + params.E_y * params.grad_phi_y) + 0.5 * params.density * (params.Q_x * params.grad_q2_y - params.Q_y * params.grad_q2_x) - params.density * (params.Q_o_x * params.grad_q2_x + params.Q_o_y * params.grad_q2_y) + params.D_x - 2 * params.D_o;
    *shear_stress_x = params.div_K_tau_x - params.S_tau_x;
    *shear_stress_y = params.div_K_tau_y - params.S_tau_y;

    // printf("%.3e %.3e %.3e %.3e\n", params.div_K_tau_x, params.S_tau_x, params.div_K_tau_y, params.S_tau_y);
    
    /* Objective function */
    // *obj = pow(*momentum_x / (velocity), 2) + pow(*momentum_y / (velocity), 2) + pow(*kinetic_energy * 1e-1 / (velocity * velocity), 2) + pow(*lateral_curvature * 1e-3 / (velocity), 2) + pow(*shear_stress_x * 1e-6, 2) + pow(*shear_stress_y * 1e-6, 2);
    // *obj = pow(*momentum_x / (velocity), 2) + pow(*momentum_y / (velocity), 2) + pow(*kinetic_energy / (velocity * velocity), 2) + pow(*lateral_curvature / (velocity), 2) + pow(*shear_stress_x, 2) + pow(*shear_stress_y, 2);
    *obj = absValue(*momentum_x) + absValue(*momentum_y); // + pow(*kinetic_energy / (velocity * velocity), 2) + pow(*lateral_curvature / (velocity), 2) + pow(*shear_stress_x, 2) + pow(*shear_stress_y, 2);

    if (*obj > 1e10) printf("%.2e %.2e %.2e %.2e\n", *obj, *shear_stress_y, params.div_K_tau_y, params.S_tau_y);

}

void calculateObjFuncGradient(int nf,
                              int n,
                              double *x,
                              double *grad,
                              int *faces_local,
                              int *faces,
                              struct VerticeConnection *vertices_connection,
                              double velocity,
                              double *facesArea,
                              double *p1, double *p2, double *p3,
                              double *e1, double *e2, double *e3,
                              double *velNorm,
                              double *velx, double *vely, double *velz,
                              double *mach,
                              double density,
                              double viscosity,
                              double *transpiration,
                              double *norm_array,
                              double *div_M, double *tau_w_x, double *tau_w_y,
                              double *max_error,
                              double *momentum_x,
                              double *momentum_y,
                              double *kinetic_energy,
                              double *lateral_curvature,
                              double *shear_stress_x,
                              double *shear_stress_y) {
    
    // Parameters
    int i, j;
    double *obj;
    struct FreestreamParameters freestream;
    struct IntegralThicknessParameters integralThickness;
    struct IntegralDefectParameters integralDefect;
    struct EquationsParameters params_aux;
    double eps;
    struct EquationsParameters *params;
    struct ProfileParameters profiles;

    // Initialize
    obj = (double*)malloc((6 * n + 1) * sizeof(double));
    eps = 1e-8;

    params = (struct EquationsParameters*)malloc(nf * sizeof(struct EquationsParameters));

    profiles.n = LAYERS;
    profiles.eta = (double*)malloc(LAYERS * sizeof(double));
    profiles.U = (double*)malloc(LAYERS * sizeof(double));
    profiles.W = (double*)malloc(LAYERS * sizeof(double));
    profiles.R = (double*)malloc(LAYERS * sizeof(double));
    profiles.S = (double*)malloc(LAYERS * sizeof(double));
    profiles.T = (double*)malloc(LAYERS * sizeof(double));
    profiles.dU_deta = (double*)malloc(LAYERS * sizeof(double));
    profiles.dW_deta = (double*)malloc(LAYERS * sizeof(double));

    freestream.density = density;
    freestream.viscosity = viscosity;

    // Reference value
    for (i = 0; i < n; i++) {
        freestream.velocity = velNorm[faces_local[i]];
        freestream.mach = mach[faces_local[i]];
        calculateEquationsParams(norm_array[0] * x[6 * i], norm_array[1] * x[6 * i + 1], norm_array[2] * x[6 * i + 2], norm_array[3] * x[6 * i + 3], norm_array[4] * x[6 * i + 4], norm_array[5] * x[6 * i + 5], &freestream, &profiles, &integralThickness, &integralDefect, &params[faces_local[i]]);
    }
    calculateDivergents(faces_local[0], faces, vertices_connection, params, facesArea[faces_local[0]], p1, p2, p3);
    calculateGradients(faces_local[0], faces, vertices_connection, params, e1, e2, e3, p1, p2, p3, velNorm, velx, vely, velz, transpiration);
    calculateObjectiveFunction(params[faces_local[0]], &obj[0], velocity, momentum_x, momentum_y, kinetic_energy, lateral_curvature, shear_stress_x, shear_stress_y);
    *max_error = obj[0];
    
    *div_M = params[faces_local[0]].div_M;
    *tau_w_x = params[faces_local[0]].tau_w_x;
    *tau_w_y = params[faces_local[0]].tau_w_y;

    for (j = 0; j < n; j++) {

        params_aux = params[faces_local[j]];

        freestream.velocity = velNorm[faces_local[j]];
        freestream.mach = mach[faces_local[j]];

        // delta
        calculateEquationsParams(norm_array[0] * (x[6 * j] + eps), norm_array[1] * x[6 * j + 1], norm_array[2] * x[6 * j + 2], norm_array[3] * x[6 * j + 3], norm_array[4] * x[6 * j + 4], norm_array[5] * x[6 * j + 5], &freestream, &profiles, &integralThickness, &integralDefect, &params[faces_local[j]]);
        calculateDivergents(faces_local[0], faces, vertices_connection, params, facesArea[faces_local[0]], p1, p2, p3);
        calculateGradients(faces_local[0], faces, vertices_connection, params, e1, e2, e3, p1, p2, p3, velNorm, velx, vely, velz, transpiration);
        calculateObjectiveFunction(params[faces_local[0]], &obj[1 + 6 * j], velocity, momentum_x, momentum_y, kinetic_energy, lateral_curvature, shear_stress_x, shear_stress_y);

        // A
        calculateEquationsParams(norm_array[0] * x[6 * j], norm_array[1] * (x[6 * j + 1] + eps), norm_array[2] * x[6 * j + 2], norm_array[3] * x[6 * j + 3], norm_array[4] * x[6 * j + 4], norm_array[5] * x[6 * j + 5], &freestream, &profiles, &integralThickness, &integralDefect, &params[faces_local[j]]);
        calculateDivergents(faces_local[0], faces, vertices_connection, params, facesArea[faces_local[0]], p1, p2, p3);
        calculateGradients(faces_local[0], faces, vertices_connection, params, e1, e2, e3, p1, p2, p3, velNorm, velx, vely, velz, transpiration);
        calculateObjectiveFunction(params[faces_local[0]], &obj[1 + 6 * j + 1], velocity, momentum_x, momentum_y, kinetic_energy, lateral_curvature, shear_stress_x, shear_stress_y);
        
        // B
        calculateEquationsParams(norm_array[0] * x[6 * j], norm_array[1] * x[6 * j + 1], norm_array[2] * (x[6 * j + 2] + eps), norm_array[3] * x[6 * j + 3], norm_array[4] * x[6 * j + 4], norm_array[5] * x[6 * j + 5], &freestream, &profiles, &integralThickness, &integralDefect, &params[faces_local[j]]);
        calculateDivergents(faces_local[0], faces, vertices_connection, params, facesArea[faces_local[0]], p1, p2, p3);
        calculateGradients(faces_local[0], faces, vertices_connection, params, e1, e2, e3, p1, p2, p3, velNorm, velx, vely, velz, transpiration);
        calculateObjectiveFunction(params[faces_local[0]], &obj[1 + 6 * j + 2], velocity, momentum_x, momentum_y, kinetic_energy, lateral_curvature, shear_stress_x, shear_stress_y);

        // Psi
        calculateEquationsParams(norm_array[0] * x[6 * j], norm_array[1] * x[6 * j + 1], norm_array[2] * x[6 * j + 2], norm_array[3] * (x[6 * j + 3] + eps), norm_array[4] * x[6 * j + 4], norm_array[5] * x[6 * j + 5], &freestream, &profiles, &integralThickness, &integralDefect, &params[faces_local[j]]);
        calculateDivergents(faces_local[0], faces, vertices_connection, params, facesArea[faces_local[0]], p1, p2, p3);
        calculateGradients(faces_local[0], faces, vertices_connection, params, e1, e2, e3, p1, p2, p3, velNorm, velx, vely, velz, transpiration);
        calculateObjectiveFunction(params[faces_local[0]], &obj[1 + 6 * j + 3], velocity, momentum_x, momentum_y, kinetic_energy, lateral_curvature, shear_stress_x, shear_stress_y);

        // Ctau1
        calculateEquationsParams(norm_array[0] * x[6 * j], norm_array[1] * x[6 * j + 1], norm_array[2] * x[6 * j + 2], norm_array[3] * x[6 * j + 3], norm_array[4] * (x[6 * j + 4] + eps), norm_array[5] * x[6 * j + 5], &freestream, &profiles, &integralThickness, &integralDefect, &params[faces_local[j]]);
        calculateDivergents(faces_local[0], faces, vertices_connection, params, facesArea[faces_local[0]], p1, p2, p3);
        calculateGradients(faces_local[0], faces, vertices_connection, params, e1, e2, e3, p1, p2, p3, velNorm, velx, vely, velz, transpiration);
        calculateObjectiveFunction(params[faces_local[0]], &obj[1 + 6 * j + 4], velocity, momentum_x, momentum_y, kinetic_energy, lateral_curvature, shear_stress_x, shear_stress_y);

        // Ctau2
        calculateEquationsParams(norm_array[0] * x[6 * j], norm_array[1] * x[6 * j + 1], norm_array[2] * x[6 * j + 2], norm_array[3] * x[6 * j + 3], norm_array[4] * (x[6 * j + 4] + eps), norm_array[5] * (x[6 * j + 5] + eps), &freestream, &profiles, &integralThickness, &integralDefect, &params[faces_local[j]]);
        calculateDivergents(faces_local[0], faces, vertices_connection, params, facesArea[faces_local[0]], p1, p2, p3);
        calculateGradients(faces_local[0], faces, vertices_connection, params, e1, e2, e3, p1, p2, p3, velNorm, velx, vely, velz, transpiration);
        calculateObjectiveFunction(params[faces_local[0]], &obj[1 + 6 * j + 5], velocity, momentum_x, momentum_y, kinetic_energy, lateral_curvature, shear_stress_x, shear_stress_y);

        params[faces_local[j]] = params_aux;

    }

    // Gradient
    for (i = 0; i < n; i++) grad[i] = (obj[1 + i] - obj[0]) / eps;

    // Free
    free(obj);
    free(params);
    free(profiles.eta);
    free(profiles.U);
    free(profiles.W);
    free(profiles.S);
    free(profiles.T);
    free(profiles.R);
    free(profiles.dU_deta);
    free(profiles.dW_deta);

}

void calculateObjFunc(int nf,
                      int n,
                      double *x,
                      double *obj,
                      int *faces_local,
                      int *faces,
                      struct VerticeConnection *vertices_connection,
                      double velocity,
                      double *facesArea,
                      double *p1, double *p2, double *p3,
                      double *e1, double *e2, double *e3,
                      double *velNorm,
                      double *velx, double *vely, double *velz,
                      double *mach,
                      double density,
                      double viscosity,
                      double *transpiration,
                      double *norm_array) {
    
    // Parameters
    int i;
    struct FreestreamParameters freestream;
    struct IntegralThicknessParameters integralThickness;
    struct IntegralDefectParameters integralDefect;
    struct EquationsParameters params_aux;
    double eps;
    struct EquationsParameters *params;
    struct ProfileParameters profiles;
    double momentum_x;
    double momentum_y;
    double kinetic_energy;
    double lateral_curvature;
    double shear_stress_x;
    double shear_stress_y;

    // Initialize
    params = (struct EquationsParameters*)malloc(nf * sizeof(struct EquationsParameters));

    profiles.n = LAYERS;
    profiles.eta = (double*)malloc(LAYERS * sizeof(double));
    profiles.U = (double*)malloc(LAYERS * sizeof(double));
    profiles.W = (double*)malloc(LAYERS * sizeof(double));
    profiles.R = (double*)malloc(LAYERS * sizeof(double));
    profiles.S = (double*)malloc(LAYERS * sizeof(double));
    profiles.T = (double*)malloc(LAYERS * sizeof(double));
    profiles.dU_deta = (double*)malloc(LAYERS * sizeof(double));
    profiles.dW_deta = (double*)malloc(LAYERS * sizeof(double));

    freestream.density = density;
    freestream.viscosity = viscosity;

    // Reference value
    for (i = 0; i < n; i++) {
        freestream.velocity = velNorm[faces_local[i]];
        freestream.mach = mach[faces_local[i]];
        calculateEquationsParams(norm_array[0] * x[6 * i], norm_array[1] * x[6 * i + 1], norm_array[2] * x[6 * i + 2], norm_array[3] * x[6 * i + 3], norm_array[4] * x[6 * i + 4], norm_array[5] * x[6 * i + 5], &freestream, &profiles, &integralThickness, &integralDefect, &params[faces_local[i]]);
    }
    calculateDivergents(faces_local[0], faces, vertices_connection, params, facesArea[faces_local[0]], p1, p2, p3);
    calculateGradients(faces_local[0], faces, vertices_connection, params, e1, e2, e3, p1, p2, p3, velNorm, velx, vely, velz, transpiration);
    calculateObjectiveFunction(params[faces_local[0]], obj, velocity, &momentum_x, &momentum_y, &kinetic_energy, &lateral_curvature, &shear_stress_x, &shear_stress_y);

    // Free
    free(params);
    free(profiles.eta);
    free(profiles.U);
    free(profiles.W);
    free(profiles.S);
    free(profiles.T);
    free(profiles.R);
    free(profiles.dU_deta);
    free(profiles.dW_deta);

}

double normVec(int n, double *a) {
    double out = 0.0;
    for (int i = 0; i < n; i++) out = out + a[i] * a[i];
    return sqrt(out);
}

void solveBoundaryLayer(int nf,
                        int nv,
                        struct VerticeConnection *vertices_connection,
                        double *vertices,
                        int *faces,
                        double *facesCenter,
                        double *facesArea,
                        double *e1, double *e2, double *e3,
                        double *p1, double *p2, double *p3,
                        double *transpiration,
                        double *delta,
                        double *A,
                        double *B,
                        double *Psi,
                        double *Ctau1,
                        double *Ctau2,
                        double *tau_x,
                        double *tau_y,
                        double *velNorm,
                        double *velx, double *vely, double *velz,
                        double *mach,
                        double density,
                        double viscosity,
                        double *cp,
                        double sound_speed,
                        double *matrix, double *array,
                        double *matrixVelx, double *matrixVely, double *matrixVelz, double *arrayVel,
                        double *doublet,
                        double freestreamNorm) {

    /* Parameters */
    double *x;
    double *x_new;
    int *div_x;
    double max_step;
    double *steps;
    double *max_ratio;
    double *max_inc;
    double *norm_array;                                 // value to normalize delta
    int n_local_faces;                                  // number of local faces
    double *x_local;                                    // parameters of local faces
    double *x_local_new;                                // parameters of local faces
    int *faces_local;                                   // local faces ids
    double *grad_local;                                 // local gradient
    double *grad_local_new;                             // local gradient

    struct FacesConnection *faces_connection;           // faces connection

    double max_error_aux;                               // maximum error
    double max_momentum_x_aux;                          // maximum momentum x error
    double max_momentum_y_aux;                          // maximum momentum x error
    double max_kinetic_energy_aux;                      // maximum kinetic energy error
    double max_lateral_curvature_aux;                   // maximum lateral curvature error
    double max_shear_stress_x_aux;                      // maximum shear stress x error
    double max_shear_stress_y_aux;                      // maximum shear stress y error
    double max_error;                                   // maximum error
    double max_momentum_x;                              // maximum momentum x error
    double max_momentum_y;                              // maximum momentum x error
    double max_kinetic_energy;                          // maximum kinetic energy error
    double max_lateral_curvature;                       // maximum lateral curvature error
    double max_shear_stress_x;                          // maximum shear stress x error
    double max_shear_stress_y;                          // maximum shear stress y error

    double max_x_value;                                 // maximum x geometric value

    int i, j, k, l;                                     // loop

    int int_max;                                        // maximum interactions

    double *div_M;                                      // mass divergence
    double *tau_w_x;                                    // x wall shear stress
    double *tau_w_y;                                    // y wall shear stress

    double alpha;

    //------------------------------------//
    /* Contant */
    const double ga_max_inter = 100;
    const double ga_ratio = (sqrt(5) + 1) / 2;

    /* Interactions */
    int ga_interaction;
    double ga_error;
    int ga_max_inter_ga;

    /* Search */
    double *ga_grad;
    double *ga_diff;
    double ga_fa, ga_fb, ga_fc, ga_fd;
    double *ga_a, *ga_b, *ga_c, *ga_d;

    /* Loop */
    int ga_j;
    //------------------------------------//
    /* Initialize */
    int_max = 5;

    alpha = 1.0;

    x = (double*)malloc(6 * nf * sizeof(double));
    x_new = (double*)malloc(6 * nf * sizeof(double));
    div_x = (int*)malloc(6 * nf * sizeof(int));

    steps = (double*)malloc(6 * sizeof(double));
    max_ratio = (double*)malloc(6 * sizeof(double));
    max_inc = (double*)malloc(6 * sizeof(double));

    max_inc[0] = 0.5; max_inc[1] = 0.1; max_inc[2] = 0.1; max_inc[3] = 0.1; max_inc[4] = 0.1; max_inc[5] = 0.1;

    max_x_value = -10.0;
    for (i = 0; i < nv; i++) if (vertices[3 * i] > max_x_value) max_x_value = vertices[3 * i];
    
    norm_array = (double*)malloc(6 * sizeof(double));
    norm_array[0] = 1e-3;
    norm_array[1] = 1.0;
    norm_array[2] = 1.0;
    norm_array[3] = 1.0;
    norm_array[4] = 1e-5;
    norm_array[5] = 1e-5;

    for (i = 0; i < nf; i++) {
        x[6 * i] = (1 / norm_array[0]) * (0.001 + 5.0 * (max_x_value - facesCenter[3 * i]) / sqrt(density * freestreamNorm * (max_x_value - facesCenter[3 * i]) / viscosity));
        x[6 * i + 1] = 1.0;
        x[6 * i + 2] = 0.0001;
        x[6 * i + 3] = 0.0001;
        x[6 * i + 4] = 0.0001;
        x[6 * i + 5] = 0.0001;
    }

    div_M = (double*)malloc(nf * sizeof(double));
    tau_w_x = (double*)malloc(nf * sizeof(double));
    tau_w_y = (double*)malloc(nf * sizeof(double));

    faces_connection = (struct FacesConnection *)malloc(nf * sizeof(struct FacesConnection));

    /* Faces connection */
    calculateFacesConnection(nv, nf, faces, vertices_connection, faces_connection);

    /* Print Interactions */
    printf("\n      Interaction   Max. error      Momentum x       Momentum y    Kinetic Energy    Lateral Curv.   Shear Stress x   Shear Stress y\n");

    /* Interaction loop */
    for (i = 1; i <= int_max; i++) {

        /* Initial error */
        max_error = -1.0;

        for (k = 0; k < 6 * nf; k++) {
            x_new[k] = 0.0;
            div_x[k] = 0;
        }

        /* Faces loop */
        for (j = 0; j < nf; j++) {

            // Local faces parameters
            n_local_faces = 1 + faces_connection[j].n;
            x_local = (double*)malloc(6 * n_local_faces * sizeof(double));
            grad_local = (double*)malloc(6 * n_local_faces * sizeof(double));
            grad_local_new = (double*)malloc(6 * n_local_faces * sizeof(double));
            faces_local = (int*)malloc(n_local_faces * sizeof(int));

            for (k = 0; k < 6; k++) x_local[k] = x[6 * j + k];
            faces_local[0] = j;
            for (k = 0; k < faces_connection[j].n; k++) {
                for (l = 0; l < 6; l++) x_local[6 * (1 + k) + l] = x[6 * faces_connection[j].faces[k] + l];
                faces_local[k + 1] = faces_connection[j].faces[k];
            }

            // Calculate gradient
            calculateObjFuncGradient(nf, n_local_faces, x_local, grad_local, faces_local, faces, vertices_connection, freestreamNorm, facesArea, p1, p2, p3, e1, e2, e3, velNorm, velx, vely, velz, mach, density, viscosity, transpiration, norm_array, &div_M[j], &tau_w_x[j], &tau_w_y[j], &max_error_aux, &max_momentum_x_aux, &max_momentum_y_aux, &max_kinetic_energy_aux, &max_lateral_curvature_aux, &max_shear_stress_x_aux, &max_shear_stress_y_aux);

            // Error
            if (max_error < max_error_aux) {
                max_error = max_error_aux;
                max_momentum_x = max_momentum_x_aux;
                max_momentum_y = max_momentum_y_aux;
                max_kinetic_energy = max_kinetic_energy_aux;
                max_lateral_curvature = max_lateral_curvature_aux;
                max_shear_stress_x = max_shear_stress_x_aux;
                max_shear_stress_y = max_shear_stress_y_aux;
            }

            // Find max. ratio components
            for (k = 0; k < n_local_faces; k++) {
                if (k == 0) {
                    max_ratio[0] = absValue(max_inc[0] / (grad_local[6 * k] + 1e-8));
                    max_ratio[1] = absValue(max_inc[1] / (grad_local[6 * k + 1] + 1e-8));
                    max_ratio[2] = absValue(max_inc[2] / (grad_local[6 * k + 2] + 1e-8));
                    max_ratio[3] = absValue(max_inc[3] / (grad_local[6 * k + 3] + 1e-8));
                    max_ratio[4] = absValue(max_inc[4] / (grad_local[6 * k + 4] + 1e-8));
                    max_ratio[5] = absValue(max_inc[5] / (grad_local[6 * k + 5] + 1e-8));
                } else {
                    if (absValue(max_inc[0] / (grad_local[6 * k] + 1e-8)) < max_ratio[0]) max_ratio[0] = absValue(max_inc[0] / (grad_local[6 * k] + 1e-8));
                    if (absValue(max_inc[1] / (grad_local[6 * k + 1] + 1e-8)) < max_ratio[0]) max_ratio[1] = absValue(max_inc[1] / (grad_local[6 * k + 1] + 1e-8));
                    if (absValue(max_inc[1] / (grad_local[6 * k + 2] + 1e-8)) < max_ratio[0]) max_ratio[2] = absValue(max_inc[2] / (grad_local[6 * k + 2] + 1e-8));
                    if (absValue(max_inc[1] / (grad_local[6 * k + 3] + 1e-8)) < max_ratio[0]) max_ratio[3] = absValue(max_inc[3] / (grad_local[6 * k + 3] + 1e-8));
                    if (absValue(max_inc[1] / (grad_local[6 * k + 4] + 1e-8)) < max_ratio[0]) max_ratio[4] = absValue(max_inc[4] / (grad_local[6 * k + 4] + 1e-8));
                    if (absValue(max_inc[1] / (grad_local[6 * k + 5] + 1e-8)) < max_ratio[0]) max_ratio[5] = absValue(max_inc[5] / (grad_local[6 * k + 5] + 1e-8));
                }
                // printf("grad: %.4e %.4e %.4e %.4e %.4e %.4e\n", grad_local[6 * k + 0], grad_local[6 * k + 1], grad_local[6 * k + 2], grad_local[6 * k + 3], grad_local[6 * k + 4], grad_local[6 * k + 5]);
            }

            // Max. step
            max_step = max_ratio[0];
            for (k = 0; k < 6; k++) if (max_step > max_ratio[k]) max_step = max_ratio[k];
            for (k = 0; k < n_local_faces; k++) if ((x_local[6 * k] - max_step * grad_local[6 * k]) < 0) max_step = 0.5 * absValue(x_local[6 * k] / (grad_local[6 * k] + 1e-8));

            //--------------------------------------------//
            // Search

            ga_diff = (double*)malloc(6 * n_local_faces * sizeof(double));
            ga_a = (double*)malloc(6 * n_local_faces * sizeof(double));
            ga_b = (double*)malloc(6 * n_local_faces * sizeof(double));
            ga_c = (double*)malloc(6 * n_local_faces * sizeof(double));
            ga_d = (double*)malloc(6 * n_local_faces * sizeof(double));
            ga_grad = (double*)malloc(6 * n_local_faces * sizeof(double));

            for (k = 0; k < 6 * n_local_faces; k++) ga_grad[k] = max_step * grad_local[k];
            ga_max_inter_ga = (int) ceil(log(1e-8 / normVec(6 * n_local_faces, ga_grad)) / log((sqrt(5) - 1) / 2));

            // First points
            for (k = 0; k < 6 * n_local_faces; k++) {
                ga_a[k] = x_local[k];
                ga_b[k] = ga_a[k] - ga_grad[k];
                ga_c[k] = ga_b[k] - (ga_b[k] - ga_a[k]) / ga_ratio;
                ga_d[k] = ga_a[k] + (ga_b[k] - ga_a[k]) / ga_ratio;
            }

            // Unidirecional search
            for (k = 0; k < ga_max_inter_ga; k++) {

                calculateObjFunc(nf, n_local_faces, ga_a, &ga_fa, faces_local, faces, vertices_connection, freestreamNorm, facesArea, p1, p2, p3, e1, e2, e3, velNorm, velx, vely, velz, mach, density, viscosity, transpiration, norm_array);
                calculateObjFunc(nf, n_local_faces, ga_b, &ga_fb, faces_local, faces, vertices_connection, freestreamNorm, facesArea, p1, p2, p3, e1, e2, e3, velNorm, velx, vely, velz, mach, density, viscosity, transpiration, norm_array);
                calculateObjFunc(nf, n_local_faces, ga_c, &ga_fc, faces_local, faces, vertices_connection, freestreamNorm, facesArea, p1, p2, p3, e1, e2, e3, velNorm, velx, vely, velz, mach, density, viscosity, transpiration, norm_array);
                calculateObjFunc(nf, n_local_faces, ga_d, &ga_fd, faces_local, faces, vertices_connection, freestreamNorm, facesArea, p1, p2, p3, e1, e2, e3, velNorm, velx, vely, velz, mach, density, viscosity, transpiration, norm_array);

                if ((ga_fa < ga_fc && ga_fa < ga_fd && ga_fa < ga_fb) || (ga_fc < ga_fa && ga_fc < ga_fd && ga_fc < ga_fb)) {
                    for (l = 0; l < 6 * n_local_faces; l++) ga_b[l] = ga_d[l];
                } else {
                    for (l = 0; l < 6 * n_local_faces; l++) ga_a[l] = ga_c[l];
                }

                // for (l = 0; l < 6 * n_local_faces; l++) ga_diff[l] = ga_b[l] - ga_a[l];
                // printf("   Inner: %d; Error.: %.4e | Func.: %.4e %.4e %.4e %.4e\n", k, normVec(6 * n_local_faces, ga_diff), ga_fa, ga_fc, ga_fd, ga_fb);

                for (l = 0; l < 6 * n_local_faces; l++) ga_c[l] = ga_b[l] - (ga_b[l] - ga_a[l]) / ga_ratio;
                for (l = 0; l < 6 * n_local_faces; l++) ga_d[l] = ga_a[l] + (ga_b[l] - ga_a[l]) / ga_ratio;

            }

            for (k = 0; k < 6; k++) {

                x_new[6 * faces_local[k]] = x_new[6 * faces_local[k]] + 0.5 * (ga_c[6 * k] + ga_d[6 * k]);
                x_new[6 * faces_local[k] + 1] = x_new[6 * faces_local[k] + 1] + 0.5 * (ga_c[6 * k + 1] + ga_d[6 * k + 1]);
                x_new[6 * faces_local[k] + 2] = x_new[6 * faces_local[k] + 2] + 0.5 * (ga_c[6 * k + 2] + ga_d[6 * k + 2]);
                x_new[6 * faces_local[k] + 3] = x_new[6 * faces_local[k] + 3] + 0.5 * (ga_c[6 * k + 3] + ga_d[6 * k + 3]);
                x_new[6 * faces_local[k] + 4] = x_new[6 * faces_local[k] + 4] + 0.5 * (ga_c[6 * k + 4] + ga_d[6 * k + 4]);
                x_new[6 * faces_local[k] + 5] = x_new[6 * faces_local[k] + 5] + 0.5 * (ga_c[6 * k + 5] + ga_d[6 * k + 5]);

                div_x[6 * faces_local[k]] = div_x[6 * faces_local[k]] + 1;
                div_x[6 * faces_local[k] + 1] = div_x[6 * faces_local[k] + 1] + 1;
                div_x[6 * faces_local[k] + 2] = div_x[6 * faces_local[k] + 2] + 1;
                div_x[6 * faces_local[k] + 3] = div_x[6 * faces_local[k] + 3] + 1;
                div_x[6 * faces_local[k] + 4] = div_x[6 * faces_local[k] + 4] + 1;
                div_x[6 * faces_local[k] + 5] = div_x[6 * faces_local[k] + 5] + 1;
            }

            free(ga_diff);
            free(ga_a);
            free(ga_b);
            free(ga_c);
            free(ga_d);
            free(ga_grad);
            //--------------------------------------------//

            // Free
            free(x_local);
            free(grad_local);
            free(grad_local_new);
            free(faces_local);
        }

        // Assing values
        for (j = 0; j < nf; j++) {
            x[6 * j] = (1 - alpha) * x[6 * j] + alpha * x_new[6 * j] / div_x[6 * j];
            x[6 * j + 1] = (1 - alpha) * x[6 * j + 1] + alpha * x_new[6 * j + 1] / div_x[6 * j + 1];
            x[6 * j + 2] = (1 - alpha) * x[6 * j + 2] + alpha * x_new[6 * j + 2] / div_x[6 * j + 2];
            x[6 * j + 3] = (1 - alpha) * x[6 * j + 3] + alpha * x_new[6 * j + 3] / div_x[6 * j + 3];
            x[6 * j + 4] = (1 - alpha) * x[6 * j + 4] + alpha * x_new[6 * j + 4] / div_x[6 * j + 4];
            x[6 * j + 5] = (1 - alpha) * x[6 * j + 5] + alpha * x_new[6 * j + 5] / div_x[6 * j + 5];
        }

        /* Print error */
        printf("           %d        %.4e      %.4e       %.4e      %.4e       %.4e      %.4e      %.4e\n", i, max_error, max_momentum_x, max_momentum_y, max_kinetic_energy, max_lateral_curvature, max_shear_stress_x, max_shear_stress_y);
        if (max_error < 1e-8)  break;

        /* Calculate inviscid parameters */
        if (i > 3000) {

            /* Transpiration */
            for (j = 0; j < nf; j++) {
                transpiration[j] = absValue(div_M[j]) / density;
                tau_x[j] = tau_w_x[j];
                tau_y[j] = tau_w_y[j];
            }

            /* Solve linear system with zero transpiration */
            calculateDoubletDistribution(nf, matrix, array, transpiration, doublet);

            /* Calculate potential surface parameters */
            calculateSurfaceParameters(nf, matrixVelx, matrixVely, matrixVelz, arrayVel, doublet, freestreamNorm, velx, vely, velz, velNorm, cp, mach, sound_speed);
        }
    }

    /* Transpiration */
    for (i = 0; i < nf; i++) {
        transpiration[j] = div_M[j] / density;
        tau_x[i] = tau_w_x[i];
        tau_y[i] = tau_w_y[i];
    }

    /* Assing values */
    for (i = 0; i < nf; i++) {
        delta[i] = norm_array[0] * x[6 * i];
        A[i] = norm_array[1] * x[6 * i + 1];
        B[i] = norm_array[2] * x[6 * i + 2];
        Psi[i] = norm_array[3] * x[6 * i + 3];
        Ctau1[i] = norm_array[4] * x[6 * i + 4];
        Ctau2[i] = norm_array[5] * x[6 * i + 5];
    }

}

/*
#####################################################
    SOLVER
#####################################################
*/
void solve(int type,
           int nv,
           int nf,
           double *vertices,
           int *faces,
           double *facesAreas,
           double *facesMaxDistance,
           double *facesCenter,
           double *controlPoints,
           double *p1,
           double *p2,
           double *p3,
           double *e1,
           double *e2,
           double *e3,
           double *freestream,
           double *sigma,
           int nSpanLeftWing,
           int nWakeLeftWing,
           int *leftWingGrid,
           double *leftWingVertices,
           int *leftWingFaces,
           int nSpanRightWing,
           int nWakeRightWing,
           int *rightWingGrid,
           double *rightWingVertices,
           int *rightWingFaces,
           int nSpanTail,
           int nWakeTail,
           int *tailGrid,
           double *tailVertices,
           int *tailFaces,
           double *doublet,
           double *velx, double *vely, double *velz, double *velNorm,
           double *cp,
           double *mach,
           double *delta,
           double *A, double *B,
           double *Psi,
           double *Ctau1, double *Ctau2,
           double *tau_x, double *tau_y, double *tau_z,
           double density, double viscosity, double sound_speed,
           double *transpiration,
           double *sigma_v,
           double *doublet_v,
           double *cp_v,
           double *velx_v, double *vely_v, double *velz_v, double *velNorm_v,
           double *transpiration_v,
           double *delta_v,
           double *A_v, double *B_v,
           double *Psi_v,
           double *Ctau1_v, double *Ctau2_v,
           double *tau_x_v, double *tau_y_v, double *tau_z_v)
{

    printf("Aerodynamic solver\n");

    printf("  - Potential flow\n");

    /* Potential flow parameters */
    double *matrix, *array;
    double *matrixVelx, *matrixVely, *matrixVelz, *arrayVel;
    double freestreamNorm;

    /* Initialize */
    matrix = (double *)malloc(nf * nf * sizeof(double));
    array = (double *)malloc(nf * sizeof(double));
    matrixVelx = (double *)malloc(nf * nf * sizeof(double));
    matrixVely = (double *)malloc(nf * nf * sizeof(double));
    matrixVelz = (double *)malloc(nf * nf * sizeof(double));
    arrayVel = (double *)malloc(nf * 3 * sizeof(double));
    freestreamNorm = sqrt(freestream[0] * freestream[0] + freestream[1] * freestream[1] + freestream[2] * freestream[2]);

    /* Potential flow system */
    printf("    > Creating linear system\n");
    createLinearSystem(nf, facesAreas, facesMaxDistance, facesCenter, controlPoints, p1, p2, p3, e1, e2, e3, freestream, sigma, nSpanLeftWing, nWakeLeftWing, leftWingGrid, leftWingVertices, leftWingFaces, nSpanRightWing, nWakeRightWing, rightWingGrid, rightWingVertices, rightWingFaces, nSpanTail, nWakeTail, tailGrid, tailVertices, tailFaces, matrix, array, matrixVelx, matrixVely, matrixVelz, arrayVel);

    /* Solve linear system with zero transpiration */
    printf("    > Solving linear system\n");
    calculateDoubletDistribution(nf, matrix, array, transpiration, doublet);

    /* Calculate potential surface parameters */
    calculateSurfaceParameters(nf, matrixVelx, matrixVely, matrixVelz, arrayVel, doublet, freestreamNorm, velx, vely, velz, velNorm, cp, mach, sound_speed);

    /* Boundary layer correction */

    /* Boundary layer parameters */
    double *tau_wall_x, *tau_wall_y;
    struct VerticeConnection *vertices_connection;

    /* Initialize */
    tau_wall_x = (double *)calloc(nf, sizeof(double));
    tau_wall_y = (double *)calloc(nf, sizeof(double));
    vertices_connection = (struct VerticeConnection *)malloc(nv * sizeof(struct VerticeConnection));

    calculateVerticesConnection(nv, nf, vertices, faces, vertices_connection);

    if (type == 1)
    {
        printf("  - Boundary layer correction\n");
        solveBoundaryLayer(nf, nv, vertices_connection, vertices, faces, facesCenter, facesAreas, e1, e2, e3, p1, p2, p3, transpiration, delta, A, B, Psi, Ctau1, Ctau2, tau_wall_x, tau_wall_y, velNorm, velx, vely, velz, mach, density, viscosity, cp, sound_speed, matrix, array, matrixVelx, matrixVely, matrixVelz, arrayVel, doublet, freestreamNorm);
        printf("\n");
    }

    /* Parameters */
    int i, j;
    struct Point e3_point, vel_point, s1, s2;
    double aux;

    // Faces
    for (i = 0; i < nf; i++)
    {

        e3_point.x = e3[3 * i];
        e3_point.y = e3[3 * i + 1];
        e3_point.z = e3[3 * i + 2];

        vel_point.x = velx[i] - transpiration[i] * e3_point.x;
        vel_point.y = vely[i] - transpiration[i] * e3_point.y;
        vel_point.z = velz[i] - transpiration[i] * e3_point.z;

        aux = norm(vel_point);

        s1.x = vel_point.x / aux;
        s1.y = vel_point.y / aux;
        s1.z = vel_point.z / aux;

        s2 = cross(s1, e3_point);

        tau_x[i] = s1.x * tau_wall_x[i] + s2.x * tau_wall_y[i];
        tau_y[i] = s1.y * tau_wall_x[i] + s2.y * tau_wall_y[i];
        tau_z[i] = s1.z * tau_wall_x[i] + s2.z * tau_wall_y[i];
    }

    // Vertices
    for (i = 0; i < nv; i++)
    {

        sigma_v[i] = 0.0;
        doublet_v[i] = 0.0;
        cp_v[i] = 0.0;
        velx_v[i] = 0.0;
        vely_v[i] = 0.0;
        velz_v[i] = 0.0;
        velNorm_v[i] = 0.0;
        transpiration_v[i] = 0.0;
        delta_v[i] = 0.0;
        A_v[i] = 0.0;
        B_v[i] = 0.0;
        Psi_v[i] = 0.0;
        Ctau1_v[i] = 0.0;
        Ctau2_v[i] = 0.0;
        tau_x_v[i] = 0.0;
        tau_y_v[i] = 0.0;
        tau_z_v[i] = 0.0;

        for (j = 0; j < vertices_connection[i].n; j++)
        {
            sigma_v[i] = sigma_v[i] + sigma[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            doublet_v[i] = doublet_v[i] + doublet[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            cp_v[i] = cp_v[i] + cp[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            velx_v[i] = velx_v[i] + velx[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            vely_v[i] = vely_v[i] + vely[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            velz_v[i] = velz_v[i] + velz[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            velNorm_v[i] = velNorm_v[i] + velNorm[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            transpiration_v[i] = transpiration_v[i] + transpiration[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            delta_v[i] = delta_v[i] + delta[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            A_v[i] = A_v[i] + A[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            B_v[i] = B_v[i] + B[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            Psi_v[i] = Psi_v[i] + Psi[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            Ctau1_v[i] = Ctau1_v[i] + Ctau1[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            Ctau2_v[i] = Ctau2_v[i] + Ctau2[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            tau_x_v[i] = tau_x_v[i] + tau_x[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            tau_y_v[i] = tau_y_v[i] + tau_y[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            tau_z_v[i] = tau_z_v[i] + tau_z[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
        }
    }

    free(matrix);
    free(array);
    free(matrixVelx);
    free(matrixVely);
    free(matrixVelz);
    free(arrayVel);
}