from ctypes import CDLL, c_int
from numpy import empty, ndarray, copy, ctypeslib, float64

class Input:

    def __init__(self, facesAreas: ndarray,
                       facesMaxDistance: ndarray,
                       facesCenter: ndarray,
                       controlPoints: ndarray,
                       p1: ndarray,
                       p2: ndarray,
                       p3: ndarray,
                       e1: ndarray,
                       e2: ndarray,
                       e3: ndarray,
                       freestream: ndarray,
                       sigma: ndarray) -> None:
        self.facesAreas = copy(facesAreas)
        self.facesMaxDistance = copy(facesMaxDistance)
        self.facesCenter = copy(facesCenter)
        self.controlPoints = copy(controlPoints)
        self.p1 = copy(p1)
        self.p2 = copy(p2)
        self.p3 = copy(p3)
        self.e1 = copy(e1)
        self.e2 = copy(e2)
        self.e3 = copy(e3)
        self.freestream = copy(freestream)
        self.sigma = copy(sigma)

class Output:

    def __init__(self, n: int) -> None:
        self.doublet = empty(n, dtype=float64)
        self.cp = empty(n, dtype=float64)
        self.velNorm = empty(n, dtype=float64)
        self.velx = empty(n, dtype=float64)
        self.vely = empty(n, dtype=float64)
        self.velz = empty(n, dtype=float64)

ND_POINTER_1 = ctypeslib.ndpointer(dtype=float64, ndim=1, flags="C")

def solve(facesAreas: ndarray,
          facesMaxDistance: ndarray,
          facesCenter: ndarray,
          controlPoints: ndarray,
          p1: ndarray,
          p2: ndarray,
          p3: ndarray,
          e1: ndarray,
          e2: ndarray,
          e3: ndarray,
          freestream: ndarray,
          sigma: ndarray) -> Output:

    # Load library
    lib = CDLL('./solver_lib.so')

    # Set input and output
    lib.solve.argtypes = [
        c_int,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
    ]

    lib.solve.restype = None
    
    # Change input shape
    n = facesAreas.shape[0]

    input = Input(facesAreas=facesAreas.reshape(facesAreas.size).astype(float64),
                  facesMaxDistance=facesMaxDistance.reshape(facesMaxDistance.size).astype(float64),
                  facesCenter=facesCenter.reshape(facesCenter.size).astype(float64),
                  controlPoints=controlPoints.reshape(controlPoints.size).astype(float64),
                  p1=p1.reshape(p1.size).astype(float64),
                  p2=p2.reshape(p2.size).astype(float64),
                  p3=p3.reshape(p3.size).astype(float64),
                  e1=e1.reshape(e1.size).astype(float64),
                  e2=e2.reshape(e2.size).astype(float64),
                  e3=e3.reshape(e3.size).astype(float64),
                  freestream=freestream,
                  sigma=sigma)

    # Create output
    output = Output(n)

    # Call library function
    lib.solve(n,
              input.facesAreas,
              input.facesMaxDistance,
              input.facesCenter,
              input.controlPoints,
              input.p1,
              input.p2,
              input.p3,
              input.e1,
              input.e2,
              input.e3,
              input.freestream,
              input.sigma,
              output.doublet,
              output.cp,
              output.velNorm,
              output.velx,
              output.vely,
              output.velz)
    
    return output