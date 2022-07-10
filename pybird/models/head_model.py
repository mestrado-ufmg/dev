from dataclasses import dataclass
from typing import Any

from pybird.models.convert import from_float, to_float

@dataclass
class HeadModel:
    h17: float = None
    h18: float = None
    h19: float = None
    delta65: float = None
    delta66: float = None

    @staticmethod
    def from_dict(obj: Any) -> 'HeadModel':
        assert isinstance(obj, dict)
        h17 = from_float(obj.get("h17"))
        h18 = from_float(obj.get("h18"))
        h19 = from_float(obj.get("h19"))
        delta65 = from_float(obj.get("delta65"))
        delta66 = from_float(obj.get("delta66"))
        return HeadModel(h17, h18, h19, delta65, delta66)

    def to_dict(self) -> dict:
        result: dict = {}
        result["h17"] = to_float(self.h17)
        result["h18"] = to_float(self.h18)
        result["h19"] = to_float(self.h19)
        result["delta65"] = to_float(self.delta65)
        result["delta66"] = to_float(self.delta66)
        return result