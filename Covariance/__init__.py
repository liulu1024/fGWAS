class BaseCovariance():
    covarType = 'Base'
    description = 'BaseCovariance'

    def __init__(self, covarType, description):
        self.covarType  = covarType
        self.description = description

    def show(obj):
        print("     Class :", obj.__class__, "\n")
        print("Covar Type :", obj.covarType, "\n")

        info = obj.get_param_info()
        print("Parameters :", info['names'], "\n")
        print("   Formula :", info['formula'], "\n")
class ParamInfo(object):
    def __init__(self,count=None,name=None,formula=None):
        self.count=count
        self.name=name
        self.formula=formula