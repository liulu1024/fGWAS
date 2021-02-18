import Simulate as simu
import Scan as scan

obj1=simu.fg_simulate("Logistic","ARMA1",20,10,range(0,3),phe_missing=0.05,snp_missing=0.05,sig_pos=5,plink_format=False,file_prefix="test",par_X=[2,3.2])
obj1_fast  =scan.fg_snpscan(obj1.obj_gen, obj1.obj_phe, method="fast", curve_type="Logistic", covariance_type="AR1",ncores=20,verbose=True)
