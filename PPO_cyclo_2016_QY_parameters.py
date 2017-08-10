import copy
import itertools
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import statsmodels.api as sm


import PTI.Corrections as PTICorr
from PTI.ReadDataFiles import PTIData
import PTI.QuantumYield as PTIQY

# <editor-fold desc="Importing">
# The file paths for the blank cyclohexane measurements
cyclo_paths = ["Henry/Emission/PPOcyclo/Jul7/cyclo2pt5g.txt"]

# The file paths for the flurophore measurements
PPO_cyclo_paths = ["Henry/Emission/PPOcyclo/Jul7/pt04mMPPOcyclo2pt5g.txt",
                   "Henry/Emission/PPOcyclo/Jul7/pt43mMPPOcyclo2pt5g.txt",
                   "Henry/Emission/PPOcyclo/Jul7/4pt3mMPPOcyclo2pt5g.txt"]


def convert_paths_to_PTIData_objs(list_of_paths):
    list_of_PTIData = list()
    for path in list_of_paths:
        list_of_PTIData.append(PTIData(path))
    return list_of_PTIData

cyclo = convert_paths_to_PTIData_objs(cyclo_paths)
PPO_cyclo = convert_paths_to_PTIData_objs(PPO_cyclo_paths)
# </editor-fold>

DEFAULT_EX_MONOCHROMATOR_SHIFT = 2.5
DEFAULT_EM_MONOCHROMATOR_SHIFT = 2.0
DEFAULT_CORRECTION_REGION_START = 365

baseline_diffs = list()
number_absorbed = list()
number_emitted = list()
qys = list()

concentrations = ['0.04 mM', '0.43 mM','4.3 mM']
f2 = open('QY Uncertainty Data/PPO_cyclo/effects_of_baseline.txt', 'w+')
f2.write("Concentration,Use of Intercept Error,Use of Slope Error,Num Absorbed,Num emitted,QY\n")
f2.close()

def QY_analysis(correction_region_start = DEFAULT_CORRECTION_REGION_START,
                ex_LUT_split = 'none', em_LUT_split = 'none',
                shift_LUT=False,
                ex_shift=DEFAULT_EX_MONOCHROMATOR_SHIFT, em_shift=DEFAULT_EM_MONOCHROMATOR_SHIFT,
                ex_LUT_interpolation = 'cubic', em_LUT_interpolation = 'cubic',
                const_diode = False,
                use_baseline_se = ('none', 'none')):

    corrected_cyclo = [PTICorr.correct_raw_to_cor(cyclo[i], baseline_fit_ranges = [[300, 305], [320, 600]],
                                                 ex_LUT_split=ex_LUT_split, em_LUT_split=em_LUT_split,
                                                 shift_LUT=shift_LUT,
                                                 ex_shift=ex_shift, em_shift=em_shift,
                                                 ex_LUT_interpolation = ex_LUT_interpolation,
                                                 em_LUT_interpolation=em_LUT_interpolation,
                                                 const_diode = const_diode,
                                                 use_baseline_se=use_baseline_se)
                      for i in range(len(cyclo))]

    corrected_PPO_cyclo = [PTICorr.correct_raw_to_cor(PPO_cyclo[i],baseline_fit_ranges = [[300, 305], [450, 600]],
                                                     ex_LUT_split=ex_LUT_split, em_LUT_split=em_LUT_split,
                                                     shift_LUT=shift_LUT,
                                                     ex_shift=ex_shift, em_shift=em_shift,
                                                     ex_LUT_interpolation = ex_LUT_interpolation,
                                                     em_LUT_interpolation=em_LUT_interpolation,
                                                     const_diode = const_diode,
                                                     use_baseline_se=use_baseline_se)
                      for i in range(len(PPO_cyclo))]
    QYs = list()
    correction_ratios = list()
    for blank, fluor in zip(3*corrected_cyclo, corrected_PPO_cyclo):
        # Define the emission region used for the quantum yield calculation
        ex_wavelength = blank.ex_range[0]
        ex_delta = 10
        ex_int_range = [ex_wavelength - ex_delta, ex_wavelength + ex_delta]

        em_int_range = [330, 450]
        correction_int_range = [correction_region_start, 450]

        # QY
        num_absorbed = PTIQY.integrate_between(blank, fluor, ex_int_range)
        num_emitted = PTIQY.integrate_between(fluor, blank, em_int_range)
        correction_area = PTIQY.integrate_between(fluor, blank, correction_int_range)

        QYs.append(num_emitted / num_absorbed)
        correction_ratios.append(correction_area / num_emitted)

        # # plt.plot(blank.wavelengths, blank.baseline-fluor.baseline, 'k')
        # baseline_diffs.append(np.mean(blank.baseline-fluor.baseline))
        # number_absorbed.append(num_absorbed)
        # number_emitted.append(num_emitted)
        # qys.append(num_emitted / num_absorbed)
        # concentration = concentrations[corrected_PPO_cyclo.index(fluor)]
        #
        # f2 = open('QY Uncertainty Data/PPO_cyclo/effects_of_baseline.txt', 'a+')
        # f2.write("%s,%s,%s,%.3f,%.3f,%.3f\n" % (concentration, use_baseline_se[0], use_baseline_se[1],
        #                                       num_absorbed, num_emitted, num_emitted/num_absorbed))
        # f2.close()
        #
        #
        #
        # plt.figure(figsize=(15,10))
        # plt.plot(blank.wavelengths, blank.baseline, 'k')
        # plt.plot(fluor.wavelengths, fluor.baseline, 'r')
        # plt.title('(intercept error, slope error) = ' + str(use_baseline_se) + '\nEmitted: %0.2f' %num_emitted + '   ' + 'Absorbed: %0.2f' %num_absorbed + '\nQY = %0.2f' %(num_emitted/num_absorbed))
        # plt.grid()
        # plt.savefig('QY Uncertainty Data/PPO_cyclo/%s -'%concentration + str(use_baseline_se) + 'BASELINES.png')
        # plt.close('all')
        #
        # plt.figure(figsize=(15, 10))
        # plt.plot(blank.wavelengths, blank.cor_data, 'k')
        # plt.plot(fluor.wavelengths, fluor.cor_data, 'r')
        # plt.axvline(x=340)
        # plt.title('(intercept error, slope error) = '+ str(use_baseline_se) + '\nEmitted: %0.2f' %num_emitted + '   ' + 'Absorbed: %0.2f' %num_absorbed + '\nQY = %0.2f' %(num_emitted/num_absorbed))
        # ex_limits_bool = (blank.wavelengths >= ex_int_range[0]) & (blank.wavelengths <= ex_int_range[1])
        # em_limits_bool = (blank.wavelengths >= em_int_range[0]) & (blank.wavelengths <= em_int_range[1])
        # plt.fill_between(blank.wavelengths, blank.cor_data, fluor.cor_data, ex_limits_bool, alpha=0.4)
        # plt.fill_between(blank.wavelengths, blank.cor_data, fluor.cor_data, em_limits_bool, alpha=0.4)
        # plt.grid()
        # plt.savefig('/home/nmcneal/Documents/WbLS/QY Uncertainty Data/PPO_cyclo/%s -'%concentration + str(use_baseline_se) + 'SPECTRA.png')
        # plt.close('all')

    correction_ratios = np.array(correction_ratios)
    first_ratio = correction_ratios[0]
    ratio_accepted_error = 0.1
    similar_ratios = correction_ratios[np.where(abs(correction_ratios - first_ratio) < ratio_accepted_error)]

    corrected_QYs = [QYs[i] * correction_ratios[i] / np.mean(similar_ratios)
                                 for i in range(len(QYs))]
    return corrected_QYs, correction_ratios


QYs = list()

LUT_interpolation_options = [(a,b) for a in ['linear', 'slinear', 'quadratic', 'cubic'] for b in ['linear', 'slinear', 'quadratic', 'cubic']]
LUT_splitting_options = [(a,b) for a in ['none', 'even', 'odd'] for b in ['none', 'even', 'odd']]
const_diode_options = [False, True]
baseline_se_options = [(a,b) for a in ['none','plus','minus'] for b in ['none','plus','minus']]
LUT_shifting_options = [False, True]
correction_region_initial_wavelengths = range(360, 370+2, 2)

def run_baseline_options():
    f = open("QY Uncertainty Data/PPO_cyclo/baseline_options.txt",'w+')
    f.write("Intercept SE,Slope SE,0.04 mM,0.43 mM,4.3 mM\n")
    for option in baseline_se_options:
        f.write(','.join(option))
        for item in QY_analysis(use_baseline_se=option)[0]:
            f.write(',' + str(item))
        f.write('\n')


    f.close()


def run_const_diode_options():
    f = open("QY Uncertainty Data/PPO_cyclo/const_diode.txt",'w+')
    f.write("Constant Diode?,0.04 mM,0.43 mM,4.3 mM\n")
    for option in const_diode_options:
        f.write(str(option))
        for item in QY_analysis(const_diode=option)[0]:
            f.write(',' + str(item))
        f.write('\n')
    f.close()


def run_LUT_interpolation_options():
    f = open("QY Uncertainty Data/PPO_cyclo/LUT_interpolation.txt",'w+')
    f.write("Ex LUT Interpolation,Em LUT Interpolation,0.04 mM,0.43 mM,4.3 mM\n")
    for option in LUT_interpolation_options:
        f.write(str(option[0])+',')
        f.write(str(option[1]))
        for item in QY_analysis(ex_LUT_interpolation=option[0],
                                em_LUT_interpolation=option[1])[0]:
            f.write(',' + str(item))
        f.write('\n')
    f.close()


def run_LUT_splitting_options():
    f = open("QY Uncertainty Data/PPO_cyclo/LUT_splitting.txt",'w+')
    f.write("Ex LUT Splitting,Em LUT Splitting,0.04 mM,0.43 mM,4.3 mM\n")
    for option in LUT_splitting_options:
        f.write(str(option[0])+',')
        f.write(str(option[1]))
        for item in QY_analysis(ex_LUT_split=option[0],
                                em_LUT_split=option[1])[0]:
            f.write(',' + str(item))
        f.write('\n')
    f.close()


def run_LUT_shifting_options():
    f = open("QY Uncertainty Data/PPO_cyclo/LUT_shifting.txt",'w+')
    f.write("Offsetting LUT?,0.04 mM,0.43 mM,4.3 mM\n")
    for option in LUT_shifting_options:
        f.write(str(option))
        for item in QY_analysis(shift_LUT=option)[0]:
            f.write(',' + str(item))
        f.write('\n')
    f.close()
    print "Finished running look-up table shifting options"


def run_correction_region_options():
    f = open("QY Uncertainty Data/PPO_cyclo/correction_region_start.txt",'w+')
    f.write("Beginning of Correction Region,0.04 mM Ratio,0.43 mM Ratio,4.3 mM Ratio,"+
            "0.04 mM,0.43 mM,4.3 mM\n")
    for option in correction_region_initial_wavelengths:
        f.write(str(option))
        qys,ratios = QY_analysis(correction_region_start=option)
        for item in np.append(ratios, qys):
            f.write(',' + str(item))
        f.write('\n')
    f.close()
    print "Finished running correction region options"


def run_all_options():
    f = open("QY Uncertainty Data/PPO_cyclo/all_options.txt",'w+')
    f.write("Shift LUT?,Intercept SE,Slope SE,Ex LUT Interpolation,Em LUT Interpolation," +
            "Ex LUT Split,Em LUT Split,Constant Diode,Start of Correction Region,0.04 mM,0.43 mM,4.3 mM\n")
    for start in correction_region_initial_wavelengths:
        for shift_LUT in LUT_shifting_options:
            for use_baseline_se in baseline_se_options:
                for ex_LUT_interpolation, em_LUT_interpolation in LUT_interpolation_options:
                    for ex_LUT_split, em_LUT_split in LUT_splitting_options:
                        for const_diode in const_diode_options:
                            f.write(str(shift_LUT) + ',')
                            f.write(use_baseline_se[0] + ',' + use_baseline_se[0] + ',')
                            f.write(ex_LUT_interpolation + ',' + em_LUT_interpolation + ',')
                            f.write(ex_LUT_split + ',' + em_LUT_split + ',')
                            f.write(str(const_diode) + ',')
                            f.write(str(start))

                            for item in QY_analysis(correction_region_start=start,
                                                    ex_LUT_interpolation= ex_LUT_interpolation,
                                                    em_LUT_interpolation=em_LUT_interpolation,
                                                    shift_LUT=shift_LUT,
                                                    ex_LUT_split=ex_LUT_split,
                                                    em_LUT_split=em_LUT_split,
                                                    const_diode=const_diode)[0]:
                                f.write(',' + str(item))
                            f.write('\n')
    f.close()

print QY_analysis()[0]

run_baseline_options()
run_const_diode_options()
run_LUT_interpolation_options()
run_LUT_splitting_options()
run_LUT_shifting_options()
run_correction_region_options()
run_all_options()

fig = plt.figure(figsize=(16,20))
axl = plt.subplot2grid((2,2), (0,0))
axr = plt.subplot2grid((2,2), (0,1), sharey=axl)
axb = plt.subplot2grid((2,2), (1,0), colspan=2)
fig.suptitle("Effect of Difference between Blank Baseline and Fluor Baseline\n(27 points - 9 error variations for each concentration", fontsize=20)


axl.scatter(baseline_diffs[::3], number_absorbed[::3], c='r', label='0.04 mM')
axl.scatter(baseline_diffs[1::3], number_absorbed[1::3], c='g', label='0.43 mM')
axl.scatter(baseline_diffs[2::3], number_absorbed[2::3], c='b', label='4.3 mM')
axl.set_title("Number of Photons Absorbed")
axl.legend()

axr.scatter(baseline_diffs[::3], number_emitted[::3], c='r', label='0.04 mM')
axr.scatter(baseline_diffs[1::3], number_emitted[1::3], c='g', label='0.43 mM')
axr.scatter(baseline_diffs[2::3], number_emitted[2::3], c='b', label='4.3 mM')
axr.set_title("Number of Photons Emitted")
axr.legend()

axb.scatter(baseline_diffs[::3], qys[::3], c='r', label='0.04 mM')
axb.scatter(baseline_diffs[1::3], qys[1::3], c='g', label='0.43 mM')
axb.scatter(baseline_diffs[2::3], qys[2::3], c='b', label='4.3 mM')
axb.set_xlabel("Average Difference between Blank Baseline and Fluor Baselines")
axb.set_title("Calculated Quantum Yield")
axb.legend()

plt.show()
