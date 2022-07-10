from dataclasses import dataclass
from typing import Any

from pybird.models.convert import from_float, from_str, from_tail_shape, to_float, to_tail_shape
from pybird.models.enums import TailShape

@dataclass
class TailModel:
    theta8: float = None
    theta9: float = None
    theta10: float = None
    h20: float = None
    h21: float = None
    h22: float = None
    delta67: float = None
    foil: str = None
    shape: TailShape = None

    @staticmethod
    def from_dict(obj: Any) -> 'TailModel':
        assert isinstance(obj, dict)
        theta8 = from_float(obj.get("theta8"))
        theta9 = from_float(obj.get("theta9"))
        theta10 = from_float(obj.get("theta10"))
        h20 = from_float(obj.get("h20"))
        h21 = from_float(obj.get("h21"))
        h22 = from_float(obj.get("h22"))
        delta67 = from_float(obj.get("delta67"))
        foil = from_str(obj.get("foil"))
        shape = to_tail_shape(obj.get("shape"))
        return TailModel(theta8, theta9, theta10, h20, h21, h22, delta67, foil, shape)

    def to_dict(self) -> dict:
        result: dict = {}
        result["theta8"] = to_float(self.theta8)
        result["theta9"] = to_float(self.theta9)
        result["theta10"] = to_float(self.theta10)
        result["h20"] = to_float(self.h20)
        result["h21"] = to_float(self.h21)
        result["h22"] = to_float(self.h22)
        result["delta67"] = to_float(self.delta67)
        result["foil"] = from_str(self.foil)
        result["shape"] = from_tail_shape(TailShape, self.shape)
        return result