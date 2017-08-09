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
# The file paths for the blank LAB measurements
# The file paths for the blank ethanol measurements
EtOH_paths = ["Henry/Sphere/PPO_ETOH/EmissionScan_ETOH_ex310_2sec_160830.txt",
              "Henry/Sphere/PPO_ETOH/EmissionScan_ETOH_ex320_2sec_160830.txt",
              "Henry/Sphere/PPO_ETOH/EmissionScan_ETOH_ex330_2sec_160830.txt",
              "Henry/Sphere/PPO_ETOH/EmissionScan_ETOH_ex340_2sec_160830.txt"]

# The file paths for the 0.31 mg/L PPO in ethanol measurements
PPO_3x14_paths = ["Henry/Sphere/PPO_ETOH/EmissionScan_3x14gperL_PPOinETOH_ex310_2sec_160831.txt",
                  "Henry/Sphere/PPO_ETOH/EmissionScan_3x14gperL_PPOinETOH_ex320_2sec_160831.txt",
                  "Henry/Sphere/PPO_ETOH/EmissionScan_3x14gperL_PPOinETOH_ex330_2sec_160831.txt",
                  "Henry/Sphere/PPO_ETOH/EmissionScan_3x14gperL_PPOinETOH_ex340_2sec_160831.txt"]


def convert_paths_to_PTIData_objs(list_of_paths):
    list_of_PTIData = list()
    for path in list_of_paths:
        list_of_PTIData.append(PTIData(path))
    return list_of_PTIData

ETOH = convert_paths_to_PTIData_objs(EtOH_paths)
PPO_3x14 = convert_paths_to_PTIData_objs(PPO_3x14_paths)
# </editor-fold>

DEFAULT_EX_MONOCHROMATOR_SHIFT = 2.5
DEFAULT_EM_MONOCHROMATOR_SHIFT = 2.0
DEFAULT_CORRECTION_REGION_START = 360


def QY_analysis(correction_region_start = DEFAULT_CORRECTION_REGION_START,
                ex_LUT_split = 'none', em_LUT_split = 'none',
                shift_LUT=False,
                ex_shift=DEFAULT_EX_MONOCHROMATOR_SHIFT, em_shift=DEFAULT_EM_MONOCHROMATOR_SHIFT,
                ex_LUT_interpolation = 'cubic', em_LUT_interpolation = 'cubic',
                const_diode = False,
                use_baseline_se = ('none', 'none')):

    ETOH_ex_wavelengths = [310, 320, 330, 340]
    corrected_ETOH = [PTICorr.correct_raw_to_cor(ETOH[i], baseline_fit_ranges = [[300, ETOH_ex_wavelengths[i] - 5], [450, 600]],
                                                 ex_LUT_split=ex_LUT_split, em_LUT_split=em_LUT_split,
                                                 shift_LUT=shift_LUT,
                                                 ex_shift=ex_shift, em_shift=em_shift,
                                                 ex_LUT_interpolation = ex_LUT_interpolation,
                                                 em_LUT_interpolation=em_LUT_interpolation,
                                                 const_diode = const_diode,
                                                 use_baseline_se=use_baseline_se)
                      for i in range(len(ETOH))]

    corrected_PPO_3x14 = [PTICorr.correct_raw_to_cor(PPO_3x14[i],baseline_fit_ranges = [[300, ETOH_ex_wavelengths[i] - 5], [450, 600]],
                                                     ex_LUT_split=ex_LUT_split, em_LUT_split=em_LUT_split,
                                                     shift_LUT=shift_LUT,
                                                     ex_shift=ex_shift, em_shift=em_shift,
                                                     ex_LUT_interpolation = ex_LUT_interpolation,
                                                     em_LUT_interpolation=em_LUT_interpolation,
                                                     const_diode = const_diode,
                                                     use_baseline_se=use_baseline_se)
                      for i in range(len(PPO_3x14))]
    QYs = list()
    correction_ratios = list()
    for blank, fluor in zip(corrected_ETOH, corrected_PPO_3x14):
        # Define the emission region used for the quantum yield calculation
        ex_wavelength = blank.ex_range[0]
        ex_delta = 10
        ex_int_range = [ex_wavelength - ex_delta, ex_wavelength + ex_delta]

        em_int_range = [320, 650]
        correction_int_range = [correction_region_start, 650]

        # QY
        num_absorbed = PTIQY.integrate_between(blank, fluor, ex_int_range)
        num_emitted = PTIQY.integrate_between(fluor, blank, em_int_range)
        correction_area = PTIQY.integrate_between(fluor, blank, correction_int_range)

        QYs.append(num_emitted / num_absorbed)
        correction_ratios.append(correction_area / num_emitted)

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
correction_region_initial_wavelengths = range(360, 365+2, 2)


def run_baseline_options():
    f = open("QY Uncertainty Data/PPO_3x14/baseline_options.txt",'w+')
    f.write("Intercept SE,Slope SE,310 nm,320 nm,330 nm,340 nm\n")
    for option in baseline_se_options:
        f.write(','.join(option))
        for item in QY_analysis(use_baseline_se=option)[0]:
            f.write(',' + str(item))
        f.write('\n')


    f.close()


def run_const_diode_options():
    f = open("QY Uncertainty Data/PPO_3x14/const_diode.txt",'w+')
    f.write("Constant Diode?,310 nm,320 nm,330 nm,340 nm\n")
    for option in const_diode_options:
        f.write(str(option))
        for item in QY_analysis(const_diode=option)[0]:
            f.write(',' + str(item))
        f.write('\n')
    f.close()


def run_LUT_interpolation_options():
    f = open("QY Uncertainty Data/PPO_3x14/LUT_interpolation.txt",'w+')
    f.write("Ex LUT Interpolation,Em LUT Interpolation,310 nm,320 nm,330 nm,340 nm\n")
    for option in LUT_interpolation_options:
        f.write(str(option[0])+',')
        f.write(str(option[1]))
        for item in QY_analysis(ex_LUT_interpolation=option[0],
                                em_LUT_interpolation=option[1])[0]:
            f.write(',' + str(item))
        f.write('\n')
    f.close()


def run_LUT_splitting_options():
    f = open("QY Uncertainty Data/PPO_3x14/LUT_splitting.txt",'w+')
    f.write("Ex LUT Splitting,Em LUT Splitting,310 nm,320 nm,330 nm,340 nm\n")
    for option in LUT_splitting_options:
        f.write(str(option[0])+',')
        f.write(str(option[1]))
        for item in QY_analysis(ex_LUT_split=option[0],
                                em_LUT_split=option[1])[0]:
            f.write(',' + str(item))
        f.write('\n')
    f.close()


def run_LUT_shifting_options():
    f = open("QY Uncertainty Data/PPO_3x14/LUT_shifting.txt",'w+')
    f.write("Offsetting LUT?,310 nm,320 nm,330 nm,340 nm\n")
    for option in LUT_shifting_options:
        f.write(str(option))
        for item in QY_analysis(shift_LUT=option)[0]:
            f.write(',' + str(item))
        f.write('\n')
    f.close()
    print "Finished running look-up table shifting options"


def run_correction_region_options():
    f = open("QY Uncertainty Data/PPO_3x14/correction_region_start.txt",'w+')
    f.write("Beginning of Correction Region,310 nm Ratio,320 nm Ratio,330 nm Ratio,340 nm Ratio,"+
            "310 nm,320 nm,330 nm,340 nm\n")
    for option in correction_region_initial_wavelengths:
        f.write(str(option))
        qys,ratios = QY_analysis(correction_region_start=option)
        for item in np.append(ratios, qys):
            f.write(',' + str(item))
        f.write('\n')
    f.close()
    print "Finished running correction region options"


def run_all_options():
    f = open("QY Uncertainty Data/PPO_3x14/all_options.txt",'w+')
    f.write("Shift LUT?,Intercept SE,Slope SE,Ex LUT Interpolation,Em LUT Interpolation," +
            "Ex LUT Split,Em LUT Split,Constant Diode,Start of Correction Region,310 nm,320 nm,330 nm,340 nm\n")
    for start in correction_region_initial_wavelengths_long_step:
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
