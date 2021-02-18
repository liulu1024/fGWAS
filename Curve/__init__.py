class BaseCurve(object):
    curveType = 'Base'
    description = 'BaseCurve'

    def __init__(self, curveType, description):
        self.curveType = curveType
        self.description = description

    def show(self,obj):
        print("     Class :", obj.__class__, "\n")
        print("Curve Type :", obj.curveType, "\n")

        info = obj.get_param_info(None)
        print("Parameters :", info.name, "\n")
        print("   Formula :", info.formula, "\n")

class ParamInfo(object):
    def __init__(self,count=None,name=None,formula=None):
        self.count=count
        self.name=name
        self.formula=formula
