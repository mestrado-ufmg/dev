from numpy import array, copy, deg2rad
from scipy.spatial.transform import Rotation as R

from pybird.models.types import Point

class BaseSystem:

    def __init__(self) -> None:
        self.x0 = array([1, 0, 0])
        self.y0 = array([0, 1, 0])
        self.z0 = array([0, 0, 1])
        return

class WingSystem:
    """Contains the base vectors of each section"""

    def __init__(self, thetaRootY: float,
                       theta1: float,
                       theta2: float,
                       theta3: float,
                       theta4: float,
                       theta5: float,
                       theta6: float,
                       theta7: float,
                       x0: Point = None,
                       y0: Point = None,
                       z0: Point = None) -> None:
        
        self._thetaRootY = thetaRootY
        self._theta1 = theta1
        self._theta2 = theta2
        self._theta3 = theta3
        self._theta4 = theta4
        self._theta5 = theta5
        self._theta6 = theta6
        self._theta7 = theta7
        
        self.x0 = x0 if x0 is not None else array([1., 0., 0.])
        self.y0 = y0 if y0 is not None else array([0., 1., 0.])
        self.z0 = z0 if z0 is not None else array([0., 0., 1.])

        self._system1()
        self._system2()
        self._system3()

        return
    
    def _system1(self) -> None:
        
        # Left wing
        x1, y1, z1 = copy(self.x0), copy(self.y0), copy(self.z0)

        # Rotate about -y1
        r = R.from_rotvec(-deg2rad(self._thetaRootY + self._theta1) * y1)
        x1 = r.apply(x1)
        z1 = r.apply(z1)

        # Rotate about x1
        r = R.from_rotvec(deg2rad(self._theta2) * x1)
        y1 = r.apply(y1)
        z1 = r.apply(z1)

        # Rotate about z1
        r = R.from_rotvec(deg2rad(self._theta3) * z1)
        x1 = r.apply(x1)
        y1 = r.apply(y1)

        # Save
        self.x1 = x1
        self.y1 = y1
        self.z1 = z1

        return
    
    def _system2(self) -> None:

        # Left wing
        x2Base, y2Base, z2Base = copy(self.x1), copy(self.y1), copy(self.z1)

        # Rotate about z2
        r = R.from_rotvec(-deg2rad(self._theta4) * z2Base)
        x2Base = r.apply(x2Base)
        y2Base = r.apply(y2Base)

        x2Tip, y2Tip, z2Tip = copy(x2Base), copy(y2Base), copy(z2Base)

        # Rotate about y2
        r = R.from_rotvec(-deg2rad(self._theta5) * y2Tip)
        x2Tip = r.apply(x2Tip)
        z2Tip = r.apply(z2Tip)

        # Save
        self.x2Base = x2Base
        self.y2Base = y2Base
        self.z2Base = z2Base
        self.x2Tip = x2Tip
        self.y2Tip = y2Tip
        self.z2Tip = z2Tip

        return
    
    def _system3(self) -> None:

        # Left wing
        x3, y3, z3 = copy(self.x2Tip), copy(self.y2Tip), copy(self.z2Tip)

        # Rotate about x3
        r = R.from_rotvec(-deg2rad(self._theta6) * x3)
        y3 = r.apply(y3)
        z3 = r.apply(z3)

        # Rotate about z3
        r = R.from_rotvec(deg2rad(self._theta7) * z3)
        x3 = r.apply(x3)
        y3 = r.apply(y3)

        # Save
        self.x3 = x3
        self.y3 = y3
        self.z3 = z3

        return

class TailSystem:
    """Contains the base vectors of each section"""

    def __init__(self, theta8: float,
                       theta9: float,
                       theta10: float,
                       x0: Point = None,
                       y0: Point = None,
                       z0: Point = None) -> None:
        
        self._theta8 = theta8
        self._theta9 = theta9
        self._theta10 = theta10
        
        self.x0 = x0 if x0 is not None else -array([1., 0., 0.])
        self.y0 = y0 if y0 is not None else -array([0., 1., 0.])
        self.z0 = z0 if z0 is not None else array([0., 0., 1.])

        self._system()

        return
    
    def _system(self) -> None:
        
        # Left wing
        x4, y4, z4 = copy(self.x0), copy(self.y0), copy(self.z0)

        # Rotate about x4
        r = R.from_rotvec(deg2rad(self._theta8) * x4)
        y4 = r.apply(y4)
        z4 = r.apply(z4)

        # Rotate about y4
        r = R.from_rotvec(deg2rad(self._theta9) * y4)
        x4 = r.apply(x4)
        z4 = r.apply(z4)

        # Rotate about z4
        r = R.from_rotvec(deg2rad(self._theta10) * z4)
        x4 = r.apply(x4)
        y4 = r.apply(y4)

        # Save
        self.x4 = x4
        self.y4 = y4
        self.z4 = z4

        return