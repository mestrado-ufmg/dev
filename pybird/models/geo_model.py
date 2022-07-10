from dataclasses import dataclass
from typing import Any

from pybird.models.convert import to_class
from pybird.models.wing_model import WingModel
from pybird.models.body_model import BodyModel
from pybird.models.head_model import HeadModel
from pybird.models.tail_model import TailModel

@dataclass
class GeoModel:
    wing: WingModel
    body: BodyModel
    head: HeadModel
    tail: TailModel

    @staticmethod
    def from_dict(obj: Any) -> 'GeoModel':
        assert isinstance(obj, dict)
        wing = WingModel.from_dict(obj.get("wing"))
        body = BodyModel.from_dict(obj.get("body"))
        head = HeadModel.from_dict(obj.get("head"))
        tail = TailModel.from_dict(obj.get("tail"))
        return GeoModel(wing, body, head, tail)

    def to_dict(self) -> dict:
        result: dict = {}
        result["wing"] = to_class(WingModel, self.wing)
        result["body"] = to_class(BodyModel, self.body)
        result["head"] = to_class(HeadModel, self.head)
        result["tail"] = to_class(TailModel, self.tail)
        return result