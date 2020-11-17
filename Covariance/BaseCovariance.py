class BaseCovariance():
    type = 'Base'
    description = 'BaseCovariance'

    def __init__(self, tpe, description):
        self.type = type
        self.description = description

    def show(obj):
        print("     Class :", obj.__class__, "\n")
        print("Curve Type :", obj.type, "\n")

        info = obj.get_param_info()
        print("Parameters :", info['names'], "\n")
        print("   Formula :", info['formula'], "\n")


def get_matrix(obj, *args):
    obj.get_matrix(*args)


def get_gradient(obj, *args):
    obj.get_gradients(*args)


def get_param_info(obj, *args):
    obj.get_param_info(*args)


def check_param(obj, *args):
    obj.check_param(*args)


def get_simu_param(obj, *args):
    obj.get_simu_param(*args)


def est_init_param(obj, *args):
    obj.est_init_param(*args)
