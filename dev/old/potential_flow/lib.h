struct input_struct {

	// Number of faces
	int n;

    // Face centers
    double **facesCenters;

    // Control points
	double **controlPoints;

    // Base vectors
	double **e1;
    double **e2;
    double **e3;

    // Vertices positions
	double **p1;
    double **p2;
    double **p3;

    // Freestream
    double *freestream;
    double density;
    double pressure;
};

struct output_struct {

    // Sigma and doublet
	double *sigma;
    double *doublet;

	// Velocity
	double *velocityNorm;
    double **velocityField;

	// Pressure
	double *pressure;
};

// Function prototype
double solver(struct input_struct *input, struct output_struct *output);