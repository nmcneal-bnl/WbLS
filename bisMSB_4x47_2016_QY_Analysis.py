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
LAB_paths = ["Henry/Sphere/bisMSB_LAB/EmissionScan_LAB_ex350_2sec_160823.txt",
             "Henry/Sphere/bisMSB_LAB/EmissionScan_LAB_ex360_2sec_160823.txt",
             "Henry/Sphere/bisMSB_LAB/EmissionScan_LAB_ex370_2sec_160823.txt",
             "Henry/Sphere/bisMSB_LAB/EmissionScan_LAB_ex380_2sec_160823.txt"]

# The file paths for the 4.47 mg/L bisMSB in LAB measurements
bisMSB_4x47_paths = ["Henry/Sphere/bisMSB_LAB/EmissionScan_bisMSBinLAB_4.47mgL_ex350_2sec_160824.txt",
                     "Henry/Sphere/bisMSB_LAB/EmissionScan_bisMSBinLAB_4.47mgL_ex360_2sec_160824.txt",
                     "Henry/Sphere/bisMSB_LAB/EmissionScan_bisMSBinLAB_4.47mgL_ex370_2sec_160824.txt",
                     "Henry/Sphere/bisMSB_LAB/EmissionScan_bisMSBinLAB_4.47mgL_ex380_2sec_160824.txt"]

def convert_paths_to_PTIData_objs(list_of_paths):
    list_of_PTIData = list()
    for path in list_of_paths:
        list_of_PTIData.append(PTIData(path))
    return list_of_PTIData

LAB = convert_paths_to_PTIData_objs(LAB_paths)
bisMSB_4x47 = convert_paths_to_PTIData_objs(bisMSB_4x47_paths)
# </editor-fold>


def QY_analysis(ex_LUT_split = 'none', em_LUT_split = 'none',
                ex_LUT_interpolation = 'cubic', em_LUT_interpolation = 'cubic',
                const_diode = False,
                use_baseline_se = ('none', 'none')):

    corrected_LAB = [PTICorr.correct_raw_to_cor(data, baseline_fit_ranges=[[300, 325], [550, 650]],
                                                ex_LUT_split=ex_LUT_split, em_LUT_split=em_LUT_split,
                                                ex_LUT_interpolation = ex_LUT_interpolation,
                                                em_LUT_interpolation=em_LUT_interpolation,
                                                const_diode = const_diode,
                                                use_baseline_se=use_baseline_se)
                     for data in LAB]
    corrected_bisMSB_4x47 = [PTICorr.correct_raw_to_cor(data, baseline_fit_ranges=[[300, 325], [550, 650]],
                                                        ex_LUT_split=ex_LUT_split, em_LUT_split=em_LUT_split,
                                                        ex_LUT_interpolation = ex_LUT_interpolation,
                                                        em_LUT_interpolation= em_LUT_interpolation,
                                                        const_diode = const_diode,
                                                        use_baseline_se=use_baseline_se)
                             for data in bisMSB_4x47]

    QYs = list()
    correction_ratios = list()
    for blank, fluor in zip(corrected_LAB, corrected_bisMSB_4x47):
        # Define the emission region used for the quantum yield calculation
        ex_wavelength = blank.ex_range[0]
        ex_delta = 5
        ex_int_range = [ex_wavelength - ex_delta, ex_wavelength + ex_delta]

        em_int_range = [365, 650]
        correction_int_range = [400, 650]

        # QY
        num_absorbed = PTIQY.integrate_between(blank, fluor, ex_int_range)
        num_emitted = PTIQY.integrate_between(fluor, blank, em_int_range)
        correction_area = PTIQY.integrate_between(fluor, blank, correction_int_range)

        QYs.append(num_emitted / num_absorbed)
        correction_ratios.append(correction_area / num_emitted)

    true_correction_ratio = np.mean(correction_ratios[:2])
    corrected_QYs = [QYs[i] * correction_ratios[i] / true_correction_ratio
                                 for i in range(len(QYs))]

    return corrected_QYs


QYs = list()

LUT_interpolation_options = [(a,b) for a in ['linear', 'slinear', 'quadratic', 'cubic'] for b in ['linear', 'slinear', 'quadratic', 'cubic']]
LUT_splitting_options = [(a,b) for a in ['none', 'even', 'odd'] for b in ['none', 'even', 'odd']]
const_diode_options = [False, True]
baseline_se_options = [(a,b) for a in ['none','plus','minus'] for b in ['none','plus','minus']]


def run_baseline_options():
    f = open("QY Uncertainty Data/bisMSB_4x47/baseline_options.txt", 'w+')
    f.write("Intercept SE, Slope SE, 350 nm, 360 nm, 370 nm, 380 nm\n")
    for option in baseline_se_options:
        f.write(','.join(option))
        for item in QY_analysis(use_baseline_se=option):
            f.write(',' + str(item))
        f.write('\n')


    f.close()


def run_const_diode_options():
    f = open("QY Uncertainty Data/bisMSB_4x47/const_diode.txt", 'w+')
    f.write("Constant Diode?, 350 nm, 360 nm, 370 nm, 380 nm\n")
    for option in const_diode_options:
        f.write(str(option))
        for item in QY_analysis(const_diode=option):
            f.write(',' + str(item))
        f.write('\n')
    f.close()


def run_LUT_interpolation_options():
    f = open("QY Uncertainty Data/bisMSB_4x47/LUT_interpolation.txt", 'w+')
    f.write("Ex LUT Interpolation, Em LUT Interpolation, 350 nm, 360 nm, 370 nm, 380 nm\n")
    for option in LUT_interpolation_options:
        f.write(str(option[0])+',')
        f.write(str(option[1]))
        for item in QY_analysis(ex_LUT_interpolation=option[0],
                                em_LUT_interpolation=option[1]):
            f.write(',' + str(item))
        f.write('\n')
    f.close()


def run_LUT_splitting_options():
    f = open("QY Uncertainty Data/bisMSB_4x47/LUT_splitting.txt", 'w+')
    f.write("Ex LUT Splitting, Em LUT Splitting, 350 nm, 360 nm, 370 nm, 380 nm\n")
    for option in LUT_splitting_options:
        f.write(str(option[0])+',')
        f.write(str(option[1]))
        for item in QY_analysis(ex_LUT_split=option[0],
                                em_LUT_split=option[1]):
            f.write(',' + str(item))
        f.write('\n')
    f.close()


def run_all_options():
    f = open("QY Uncertainty Data/bisMSB_4x47/all_options.txt", 'w+')
    f.write("Intercept SE, Slope SE, Ex LUT Interpolation, Em LUT Interpolation," +
            "Ex LUT Split, Em LUT Split, Constant Diode, 350 nm, 360 nm, 370 nm, 380 nm\n")
    for use_baseline_se in baseline_se_options:
        for ex_LUT_interpolation, em_LUT_interpolation in LUT_interpolation_options:
            for ex_LUT_split, em_LUT_split in LUT_splitting_options:
                for const_diode in const_diode_options:
                    f.write(use_baseline_se[0] + ',' + use_baseline_se[0] + ',')
                    f.write(ex_LUT_interpolation + ',' + em_LUT_interpolation + ',')
                    f.write(ex_LUT_split + ',' + em_LUT_split + ',')
                    f.write(str(const_diode))

                    for item in QY_analysis(use_baseline_se=use_baseline_se,
                                            ex_LUT_interpolation= ex_LUT_interpolation,
                                            em_LUT_interpolation=em_LUT_interpolation,
                                            ex_LUT_split=ex_LUT_split,
                                            em_LUT_split=em_LUT_split,
                                            const_diode=const_diode):
                        f.write(',' + str(item))
                    f.write('\n')
    f.close()

# print QY_analysis()

run_baseline_options()
# run_const_diode_options()
# run_LUT_interpolation_options()
# run_LUT_splitting_options()
run_all_options()