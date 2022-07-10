import sys

sys.path.append('./')

from models.case_model import CaseModel

if __name__ == '__main__':
    case = CaseModel.from_file('./data/data.case')

    print(case.geo.wing.l0)