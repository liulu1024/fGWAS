class BaseCurve(object):
    curveType = 'Base'
    description = 'BaseCurve'

    def __init__(self, curveType, description):
        self.curveType = curveType
        self.description = description

    def show(self):
        print("     Class :", self.__class__, "\n")
        print("Curve Type :", self.curveType, "\n")

        info = self.get_param_info()
        print("Parameters :", info['names'], "\n")
        print("   Formula :", info['formula'], "\n")


