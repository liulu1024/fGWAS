import abc


class AbstractCurve(metaclass=abc.ABCMeta):
    curveType = 'Base'
    description = 'BaseCurve'

    def show(self):
        pass

    @abc.abstractmethod
    def get_curve(self, param, times, options):
        pass


    @abc.abstractmethod
    def get_param_info(self):
        pass

    @abc.abstractmethod
    def show(self):
        pass

    @abc.abstractmethod
    def get_gradient(self):
        pass

    @abc.abstractmethod
    def check_param(self):
        pass

    @abc.abstractmethod
    def get_simu_param(self):
        pass

    @abc.abstractmethod
    def est_init_param(self):
        pass
