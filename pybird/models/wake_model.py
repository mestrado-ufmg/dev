from dataclasses import dataclass

from numpy import ndarray


@dataclass
class WakeModel:
    vertices: ndarray
    grid: ndarray
    faces: ndarray