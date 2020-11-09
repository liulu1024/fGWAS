import Curve.Logistic
import Curve.CurveTools as tool
class testCurve(object):

    test_log = Curve.Log.log()

    print( tool.CurveTools.get_curve(test_log,10))
    print(tool.CurveTools.get_param_info(test_log,"test"))

