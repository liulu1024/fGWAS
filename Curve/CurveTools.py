import numpy as np
import torch


class CurveTools(object):
    def get_curve(obj, *args):
        return obj.get_curve(*args)

    def get_param_info(obj, *args):
        return obj.get_param_info(*args)

    def get_gradient(obj, *args):
        return obj.get_gradient(*args)

    def check_param(obj, *args):
        return obj.check_param(*args)

    def get_simu_param(obj, *args):
        return obj.get_simu_param(*args)

    def est_init_param(obj, *args):
        return obj.est_init_param(*args)

    # Common Tool
    '''
    pheT: 1:time.points np.array 便于标量向量运算
    pheY: numpy
    '''

    def get_mean_vector(self, pheY, pheT):
        # mean t
        t_count = len(np.unique(pheT))
        if t_count >= 20:
            t_min = min(pheT)
            t_max = max(pheT)
        pheT = np.round((pheT - t_min) / (t_max - t_min) * 20) / 20 * (t_max - t_min) + t_min
        t_all = np.sort(np.unique(pheT))
        # mean y
        y_all = []
        pheT = torch.tensor(pheT)
        for t in t_all:
            y_all = y_all.append(np.mean(pheY[pheT == t]))
        return {'t': t_all, 'y': y_all}

    '''
    all curve 
    '''
    def get_all_curve(self):
        pass

