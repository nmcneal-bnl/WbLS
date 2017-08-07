itimport copy
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

# The file paths for the blank ethanol measurements
EtOH_paths = ["Henry/Sphere/PPO_ETOH/EmissionScan_ETOH_ex310_2sec_160830.txt",
              "Henry/Sphere/PPO_ETOH/EmissionScan_ETOH_ex320_2sec_160830.txt",
              "Henry/Sphere/PPO_ETOH/EmissionScan_ETOH_ex330_2sec_160830.txt",
              "Henry/Sphere/PPO_ETOH/EmissionScan_ETOH_ex340_2sec_160830.txt"]

# The file paths for the 0.31 mg/L PPO in ethanol measurements
PPO_0x31_paths = ["Henry/Sphere/PPO_ETOH/EmissionScan_0x31gperL_PPOinETOH_ex310_2sec_160831.txt",
                  "Henry/Sphere/PPO_ETOH/EmissionScan_0x31gperL_PPOinETOH_ex320_2sec_160831.txt",
                  "Henry/Sphere/PPO_ETOH/EmissionScan_0x31gperL_PPOinETOH_ex330_2sec_160831.txt",
                  "Henry/Sphere/PPO_ETOH/EmissionScan_0x31gperL_PPOinETOH_ex340_2sec_160831.txt"]

# The file paths for the 3.14 mg/L PPO in ethanol measurements
PPO_3x14_paths = ["Henry/Sphere/PPO_ETOH/EmissionScan_3x14gperL_PPOinETOH_ex310_2sec_160831.txt",
                  "Henry/Sphere/PPO_ETOH/EmissionScan_3x14gperL_PPOinETOH_ex320_2sec_160831.txt",
                  "Henry/Sphere/PPO_ETOH/EmissionScan_3x14gperL_PPOinETOH_ex330_2sec_160831.txt",
                  "Henry/Sphere/PPO_ETOH/EmissionScan_3x14gperL_PPOinETOH_ex340_2sec_160831.txt"]

'''Requires: A list of file paths for PTI data files
   Effect:   Takes each path and uses it to read into a PTIData instance
   Returns:  A list of PTIData instances. Each element comes from the corresponding
             path in the list argument.'''
def convert_paths_to_PTIData_objs(list_of_paths):
    list_of_PTIData = list()
    for path in list_of_paths:
        list_of_PTIData.append(PTIData(path))
    return list_of_PTIData

LAB = convert_paths_to_PTIData_objs(LAB_paths)
bisMSB_4x47 = convert_paths_to_PTIData_objs(bisMSB_4x47_paths)

ETOH = convert_paths_to_PTIData_objs(EtOH_paths)
PPO_0x31 = convert_paths_to_PTIData_objs(PPO_0x31_paths)
PPO_3x14 = convert_paths_to_PTIData_objs(PPO_3x14_paths)
# </editor-fold>

def QY_analysis(ex_LUT_split, em_LUT_split,
                LUT_interpolation, const_diode):
    corrected_LAB = [PTICorr.correct_raw_to_cor(data, baseline_fit_ranges=[[300, 325], [550, 650]],
                                                ex_LUT_split=ex_LUT_split, em_LUT_split=em_LUT_split,
                                                LUT_interpolation = LUT_interpolation,
                                                const_diode = const_diode)
                     for data in LAB]
    corrected_bisMSB_4x47 = [PTICorr.correct_raw_to_cor(data, baseline_fit_ranges=[[300, 325], [550, 650]])
                             for data in bisMSB_4x47]

    ETOH_ex_wavelengths = [310, 320, 330, 340]
    corrected_ETOH = [
        PTICorr.correct_raw_to_cor(ETOH[i], baseline_fit_ranges=[[300, ETOH_ex_wavelengths[i] - 5], [500, 650]])
        for i in range(len(ETOH))]
    corrected_PPO_0x31 = [
        PTICorr.correct_raw_to_cor(PPO_0x31[i], baseline_fit_ranges=[[300, ETOH_ex_wavelengths[i] - 5], [500, 650]])
        for i in range(len(ETOH))]
    corrected_PPO_3x14 = [
        PTICorr.correct_raw_to_cor(PPO_3x14[i], baseline_fit_ranges=[[300, ETOH_ex_wavelengths[i] - 5], [500, 650]])
        for i in range(len(ETOH))]

    bisMSB_4x47_QYs = list()
    bisMSB_4x47_corrections = list()
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

        bisMSB_4x47_QYs.append(num_emitted / num_absorbed)
        bisMSB_4x47_corrections.append(correction_area / num_emitted)

    bisMSB_4x47_QY_correction = np.mean(bisMSB_4x47_corrections[:2])
    bisMSB_4x47_corrected_QYs = [bisMSB_4x47_QYs[i] * bisMSB_4x47_corrections[i] / bisMSB_4x47_QY_correction
                                 for i in range(len(bisMSB_4x47_QYs))]

    return bisMSB_4x47_corrected_QYs

QYs = list()
for ex_LUT_split in ['none', 'even', 'odd']:
    for em_LUT_split in ['none', 'even', 'odd']:
        for LUT_interpolation in ['linear', 'slinear', 'quadratic', 'cubic']:
            for const_diode in [False, True]:
                QYs.append(QY_analysis(ex_LUT_split, em_LUT_split, LUT_interpolation, const_diode))

QYs= np.array(QYs)
mean_QYs = np.mean(QYs, axis = 0)
se_means = np.std(QYs, axis = 0) / np.sqrt(np.shape(QYs)[0])

print "Quantum Yield Means:", mean_QYs
print "Quantum Yield SE of the Mean:", se_means
print "95% Confidence intervals:\n", np.column_stack((mean_QYs - 1.96*se_means, mean_QYs + 1.96*se_means))



