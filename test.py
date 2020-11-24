import Simulation.simulate

# base.get_all_curve()
# base.get_all_curve()
# base.get_all_curve()
Simulation.simulate.fg_simulate("Logistic", "ARMA1", 200, 100, range(0, 8), phe_missing=0.05, snp_missing=0.05,
                                sig_pos=51, plink_format=False, file_prefix='test', par_X=(2, 3.2))
