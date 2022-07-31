struct input_struct {

	// Number of faces
	int n;

    // Vertices
    double **vertices;

    // Faces
    int **faces;

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
};

struct output_struct {

    // Sigma and doublet
	double **matrix;
    double *array;
    
};

struct output_vel_struct {

	double **vel;
    double *velNorm;
    double *pressure;
    
};

// Function prototype
double division(double a, double b);
void create(struct input_struct *input, struct output_struct *output);
void parameters(struct input_struct *input, struct output_vel_struct *output, double *sigma);