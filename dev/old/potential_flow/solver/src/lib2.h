struct input_struct {

	// Number of faces
	int nf;

    // Vertices
    double **vertices;

    // Faces
    int **faces;
    double *facesAreas;
    double *facesMaxDistance;

    // Face centers
    double **facesCenter;

    // Control points
	double **controlPoints;

    // Base vectors
	double **e1;
    double **e2;
    double **e3;

    // Freestream
    double *freestream;

    // Sigma
    double *sigma;
    
    // // Wake
    // int nWakeLeftWing;
    // int nSpanLeftWing;
    // double **leftWingWakeVertices;
    // int **leftWingWakeGrid;
    // int **leftWingWakeFaces;
    // double **leftWingTrailingVectors;

    // int nWakeRightWing;
    // int nSpanRightWing;
    // double **rightWingWakeVertices;
    // int **rightWingWakeGrid;
    // int **rightWingWakeFaces;
    // double **rightWingTrailingVectors;

    // int nWakeTail;
    // int nSpanTail;
    // double **tailWakeVertices;
    // int **tailWakeGrid;
    // int **tailWakeFaces;
    // double **tailTrailingVectors;
};

struct output_struct {

	double **matrix;
    double *array;
    double ***velMatrix;
    double **velArray;
    
};

struct Point3D {

	double x;
    double y;
    double z;
    
};

struct Point2D {

	double x;
    double y;
    
};

void create(int nf,
            int nv,
            double *vertices,
            int *faces,
            double *facesAreas,
            double *facesMaxDistance,
            double *facesCenter,
            double *controlPoints,
            double *e1, double *e2, double *e3,
            double *freestream,
            double *sigma,
            double *matrix,
            double *array,
            double *velMatrixX,
            double *velMatrixY,
            double *velMatrixZ,
            double *velArray);
double division(double a, double b);
double norm(struct Point3D p);
struct Point3D cross(struct Point3D p1, struct Point3D p2);
void source(struct Point3D p, struct Point2D p1, struct Point2D p2, struct Point2D p3, struct Point3D e1, struct Point3D e2, struct Point3D e3, double area, double maxDistance, double* vel);
void line(struct Point3D p, struct Point3D p1, struct Point3D p2, double* vel);
void doublet(struct Point3D p, struct Point3D p1, struct Point3D p2, struct Point3D p3, struct Point3D e1, struct Point3D e2, struct Point3D e3, double area, double maxDistance, double* vel);