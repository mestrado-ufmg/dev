import numpy as np

ZERO_ERROR = 1e-8

def division(a: float, b: float) -> float:
    if -ZERO_ERROR < b < ZERO_ERROR:
        sign = 1 if b > 0 else -1
        return sign * a / ZERO_ERROR
    else:
        return a / b

def source_velocity(sigma: float = None,
                    p1: np.ndarray = None,
                    p2: np.ndarray = None,
                    p3: np.ndarray = None,
                    p: np.ndarray = None,
                    e1: np.ndarray = None,
                    e2: np.ndarray = None,
                    e3: np.ndarray = None) -> np.ndarray:
    """
    Calculate the velocity indulced by a triangular panel source
    
    Parameters
    ----------

    - center: center of the panel
    - p1, p2 and p3: vertices of the panel
    - p: point where the velocity is indulced
    - e1, e2 and e3: base vectors

    Output
    ------
    - Velocity
    """

    # Params
    r1 = ((p[0] - p1[0]) ** 2 + (p[1] - p1[1]) ** 2 + p[2] ** 2) ** 0.5
    r2 = ((p[0] - p2[0]) ** 2 + (p[1] - p2[1]) ** 2 + p[2] ** 2) ** 0.5
    r3 = ((p[0] - p3[0]) ** 2 + (p[1] - p3[1]) ** 2 + p[2] ** 2) ** 0.5
    l1 = (p[0] - p1[0]) ** 2 + p[2] ** 2
    l2 = (p[0] - p2[0]) ** 2 + p[2] ** 2
    l3 = (p[0] - p3[0]) ** 2 + p[2] ** 2
    h1 = (p[0] - p1[0]) * (p[1] - p1[1])
    h2 = (p[0] - p2[0]) * (p[1] - p2[1])
    h3 = (p[0] - p3[0]) * (p[1] - p3[1])

    # Point 1 and 2
    d12 = ((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2) ** 0.5
    m12 = division((p2[1] - p1[1]), (p2[0] - p1[0]))

    # Point 2 and 3
    d23 = ((p3[0] - p2[0]) ** 2 + (p3[1] - p2[1]) ** 2) ** 0.5
    m23 = division((p3[1] - p2[1]), (p3[0] - p2[0]))

    # Point 3 and 1
    d31 = ((p1[0] - p3[0]) ** 2 + (p1[1] - p3[1]) ** 2) ** 0.5
    m31 = division((p1[1] - p3[1]), (p1[0] - p3[0]))

    # Velocities
    ln12 = np.log( division((r1 + r2 - d12), (r1 + r2 + d12)) )
    ln23 = np.log( division((r2 + r3 - d23), (r2 + r3 + d23)) )
    ln31 = np.log( division((r3 + r1 - d31), (r3 + r1 + d31)) )

    mult = (sigma / (4 * np.pi))

    u = -mult * ( division((p2[1] - p1[1]), d12) * ln12 + division((p3[1] - p2[1]), d23) * ln23 + division((p1[1] - p3[1]), d31) * ln31 )
    v = mult * ( division((p2[0] - p1[0]), d12) * ln12 + division((p3[0] - p2[0]), d23) * ln23 + division((p1[0] - p3[0]), d31) * ln31 )
    w = -mult * ( np.arctan(division(m12 * l1 - h1, p[2] * r1)) - np.arctan(division(m12 * l2 - h2, p[2] * r2)) + np.arctan(division(m23 * l2 - h2, p[2] * r2)) - np.arctan(division(m23 * l3 - h3, p[2] * r3)) + np.arctan(division(m31 * l3 - h3, p[2] * r3)) - np.arctan(division(m31 * l1 - h1, p[2] * r1)) )
    
    return u * e1 + v * e2 + w * e3