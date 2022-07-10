from dataclasses import dataclass
from json import dump, load
from typing import Any

from pybird.models.geo_model import GeoModel
from pybird.models.convert import from_str, to_class

@dataclass
class CaseModel:
    name: str
    description: str
    geo: GeoModel

    @staticmethod
    def from_dict(obj: Any) -> 'CaseModel':
        assert isinstance(obj, dict)
        name = from_str(obj.get("name"))
        description = from_str(obj.get("description"))
        geo = GeoModel.from_dict(obj.get("geo"))
        return CaseModel(name, description, geo)

    def to_dict(self) -> dict:
        result: dict = {}
        result["name"] = from_str(self.name)
        result["description"] = from_str(self.description)
        result["geo"] = to_class(GeoModel, self.geo)
        return result
    
    @staticmethod
    def from_file(file: str) -> 'CaseModel':
        f = open(file)
        data = load(f)
        f.close()
        return CaseModel.from_dict(data)
    
    def to_file(self, file: str) -> None:
        out_file = open(file, "w") 
        dump(self.to_dict(), out_file, indent=2)
        out_file.close()
        return
    

