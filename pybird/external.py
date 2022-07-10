from pybird.models.case_model import CaseModel
from pybird.models.geo_model import GeoModel
from pybird.models.wing_model import WingModel
from pybird.models.body_model import BodyModel
from pybird.models.head_model import HeadModel
from pybird.models.tail_model import TailModel
from pybird.geo.geo import Geo
from pybird.mesh.mesh import Mesh
from pybird.view.view import View

class model:

    def __init__(self, name: str, description: str=None) -> None:
        self.data = CaseModel(name, description if description is not None else '', GeoModel(WingModel(), BodyModel(), HeadModel(), TailModel()))
        self.geo = Geo(self.data.geo)
        self.mesh = Mesh(self.geo)
        self.view = View(self.mesh)
        return

    def save(self, file: str) -> None:
        self.data.to_file(file)
        return
    
    def load(self, file: str) -> None:
        self.data = CaseModel.from_file(file)
        self.geo = Geo(self.data.geo)
        self.mesh = Mesh(self.geo)
        self.view = View(self.mesh)
        return