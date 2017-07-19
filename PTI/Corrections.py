'''Takes a PTIData instance and the excorr and emcorr LUTs as numpy arrays.
   Removes the diode, excitation LUT, and emission LUT corrections from the corrected data.
   Modifies the instances 'cor_data' member variable array.'''

import copy
from ReadDataFiles import PTIData
from scipy.interpolate import interp1d
import numpy as np

def load_excorr_file(PTIData_instance, interp_method = 'cubic'):
    excorr = np.genfromtxt('PTI/correction_data/excorr.txt',
                           skip_header = 6,
                           skip_footer = 1,
                           usecols = 1)
    
    step = PTIData_instance.step_size
    excorr_range = np.arange(250, 750+1)
    xvals = np.arange(250, 750 + step , step)
    excorr = interp1d(excorr_range, excorr, interp_method)(xvals)
    
    min_data_wavelength = PTIData_instance.wavelengths[0]
    max_data_wavelength = PTIData_instance.wavelengths[-1]
    
    needed_wavelengths = np.where((xvals >= min_data_wavelength) & 
                                  (xvals <= max_data_wavelength))
    
    excorr = excorr[needed_wavelengths]
    
    if min_data_wavelength < 250:
        extra_x = np.arange(min_data_wavelength, 250, step)
        left_interp = excorr[0]*np.ones(extra_x.size)
        excorr = np.append(left_interp, excorr)
    if max_data_wavelength > 750:
        extra_x = np.arange(750, max_data_wavelength, step)
        right_interp = excorr[-1]*np.ones(extra_x.size)
        excorr = np.append(excorr, right_interp)
    
    return excorr

def load_emcorr_file(PTIData_instance, interp_method = 'cubic',
                     in_FS = False):
    if in_FS:
        fname = 'PTI/correction_data/emcorri.txt'
        emcorr_min_wave = 250
        emcorr_max_wave = 850
    if not in_FS:
        fname = 'PTI/correction_data/emcorr-sphere-quanta.txt'
        emcorr_min_wave = 300
        emcorr_max_wave = 848
    
    emcorr = np.genfromtxt(fname,
                           skip_header = 6,
                           skip_footer = 1,
                           usecols = 1)
    
    step = PTIData_instance.step_size
    emcorr_range = np.arange(emcorr_min_wave, emcorr_max_wave+2,2)
    xvals = np.arange(emcorr_min_wave, emcorr_max_wave + step , step )
    
    emcorr = interp1d(emcorr_range, emcorr, interp_method)(xvals)

    min_data_wavelength = PTIData_instance.wavelengths[0]
    max_data_wavelength = PTIData_instance.wavelengths[-1]
    
    needed_wavelengths = np.where((xvals >= min_data_wavelength) & 
                                  (xvals <= max_data_wavelength))
    
    emcorr = emcorr[needed_wavelengths]
    
    if min_data_wavelength < emcorr_min_wave:
        extra_x = np.arange(min_data_wavelength, emcorr_min_wave, step)
        left_interp = emcorr[0]*np.ones(extra_x.size)
        emcorr = np.append(left_interp, emcorr)
    if max_data_wavelength > emcorr_max_wave:
        extra_x = np.arange(emcorr_max_wave, max_data_wavelength, step)
        right_interp = emcorr[-1]*np.ones(extra_x.size)
        emcorr = np.append(emcorr, right_interp)
    
    return emcorr
    
def get_corrections(PTIData_instance, 
                          interp_method = 'cubic',
                          in_FS = False,
                          diode = True,
                          excorr = True,
                          emcorr = True):
    
    # Create a copy of the data instance
    corrections = np.ones(PTIData_instance.wavelengths.size)
    
    # Perform the diode correction made from the ExCorr RCQC signal
    if diode:
        corrections  /= PTIData_instance.diode
    
    # Perform the LUT corrections
    if excorr:
        corrections /= load_excorr_file(PTIData_instance)
    
    if emcorr:
        corrections *= load_emcorr_file(PTIData_instance, in_FS)
    
    return corrections


def decorrect_data(PTIData_instance, 
                           interp_method = 'cubic',
                           in_FS = False,
                           diode = True,
                           excorr = True,
                           emcorr = True):

    # Create a copy of the data instance
    data = copy.deepcopy(PTIData_instance)
    
    corrections = get_corrections(PTIData_instance, interp_method, in_FS,
                                  diode, excorr, emcorr)
    
    data.raw_data = data.cor_data / corrections

    return data

def correct_data(PTIData_instance, 
                           interp_method = 'cubic',
                           in_FS = False,
                           diode = True,
                           excorr = True,
                           emcorr = True):

    # Create a copy of the data instance
    data = copy.deepcopy(PTIData_instance)
    
    corrections = get_corrections(PTIData_instance, interp_method, in_FS,
                                  diode, excorr, emcorr)
    
    data.cor_data = data.raw_data * corrections

    return data




    
