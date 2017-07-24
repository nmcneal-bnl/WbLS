import copy
from scipy.interpolate import interp1d
import numpy
import statsmodels.api as sm

from ReadDataFiles import PTIData


''' Uses a least-squared routine to fit  data to a linear function. 
    The fit is performed over the given wavelength ranges.
    Returns a tuple containing: ('''
def linear_baseline_params(PTIData, list_of_ranges):

    X = numpy.array([])
    Y = numpy.array([])
    for arange in list_of_ranges:
        start = arange[0]
        end = arange[1]

        X = numpy.append(X, numpy.arange(start, end+PTIData.step_size, PTIData.step_size))

        select_by_wavelength = numpy.where((PTIData.wavelengths >= start) &
                                           (PTIData.wavelengths <= end))
        Y = numpy.append(Y, PTIData.raw_data[select_by_wavelength])

    X = sm.add_constant(X)
    model = sm.OLS(Y, X)
    results = model.fit()
   
    return (results.params, results.bse)


''' Uses a least-squared routine to fit  data to a polynomial of a given
    degree. The fit is performed over the given wavelength ranges.
    Returns an array with the polynomial evaluated over the full spectrum range
    with the same step size.'''
def linear_baseline(PTIData, list_of_ranges,
                    use_incpt_se = 'none', use_slope_se = 'none'):

    # Calculate the parameters for the fit and their standard errors. Order: (intercept, slope)
    fit_params, errors = linear_baseline_params(PTIData, list_of_ranges)
    fit_params = list(fit_params)
    errors = list(errors)

    # Use the standard errors to modify either parameter
    incpt ={'none': fit_params[0],
            'minus': fit_params[0] - errors[0],
            'plus': fit_params[0] + errors[0]}

    slope ={'none': fit_params[1],
            'minus': fit_params[1] - errors[1],
            'plus': fit_params[1] + errors[1]}

    baseline = slope[use_slope_se] * PTIData.wavelengths + incpt[use_incpt_se]

    # Return the baseline array as well as the fit parameters and errors.
    return baseline, fit_params, errors

''''''
def load_excorr_file(PTIData_instance, interp_method = 'cubic', split = 'none'):
    
    excorr = numpy.genfromtxt('PTI/correction_data/excorr.txt',
                           skip_header = 6,
                           skip_footer = 1,
                           usecols = 1)
    LUT_step = 1
    LUT_start = 250
    LUT_end = 750

    if split.lower() == 'none':
        pass
    elif split.lower() == 'even':
        excorr =  excorr[::2]
        LUT_step *= 2
    elif split.lower() == 'odd':
        excorr =  excorr[1::2]
        LUT_step *= 2
        LUT_start += 1
        LUT_end -= 1
    else:  
        print "ERROR: Not a valid method for splitting LUT"
        return None

    step = PTIData_instance.step_size
    excorr_range = numpy.arange(LUT_start, LUT_end + LUT_step, LUT_step)
    xvals = numpy.arange(LUT_start, LUT_end + step , step)
    excorr = interp1d(excorr_range, excorr, interp_method)(xvals)
    
    min_data_wavelength = PTIData_instance.wavelengths[0]
    max_data_wavelength = PTIData_instance.wavelengths[-1]
    
    needed_wavelengths = numpy.where((xvals >= min_data_wavelength) & 
                                     (xvals <= max_data_wavelength))
    
    excorr = excorr[needed_wavelengths]
    
    if min_data_wavelength < LUT_start:
        extra_x = numpy.arange(min_data_wavelength, LUT_start, step)
        left_interp = excorr[0]*numpy.ones(extra_x.size)
        excorr = numpy.append(left_interp, excorr)
    if max_data_wavelength > LUT_end:
        extra_x = numpy.arange(LUT_end, max_data_wavelength, step)
        right_interp = excorr[-1]*numpy.ones(extra_x.size)
        excorr = numpy.append(excorr, right_interp)
    
    return excorr

''''''
def load_emcorr_file(PTIData_instance, interp_method = 'cubic', FS = False):
    if FS:
        fname = 'PTI/correction_data/emcorri.txt'
        emcorr_min_wave = 250
        emcorr_max_wave = 850
    if not FS:
        fname = 'PTI/correction_data/emcorr-sphere-quanta.txt'
        emcorr_min_wave = 300
    
        emcorr_max_wave = 848
    emcorr = numpy.genfromtxt(fname,
                              skip_header = 6,
                              skip_footer = 1,
                              usecols = 1)
    
    step = PTIData_instance.step_size
    emcorr_range = numpy.arange(emcorr_min_wave, emcorr_max_wave+2,2)
    xvals = numpy.arange(emcorr_min_wave, emcorr_max_wave + step , step )
    
    emcorr = interp1d(emcorr_range, emcorr, interp_method)(xvals)

    min_data_wavelength = PTIData_instance.wavelengths[0]
    max_data_wavelength = PTIData_instance.wavelengths[-1]
    
    needed_wavelengths = numpy.where((xvals >= min_data_wavelength) & 
                                     (xvals <= max_data_wavelength))
    
    emcorr = emcorr[needed_wavelengths]
    
    if min_data_wavelength < emcorr_min_wave:
        extra_x = numpy.arange(min_data_wavelength, emcorr_min_wave, step)
        left_interp = emcorr[0]*numpy.ones(extra_x.size)
        emcorr = numpy.append(left_interp, emcorr)
    if max_data_wavelength > emcorr_max_wave:
        extra_x = numpy.arange(emcorr_max_wave, max_data_wavelength, step)
        right_interp = emcorr[-1]*numpy.ones(extra_x.size)
        emcorr = numpy.append(emcorr, right_interp)
    
    return emcorr
    

''''''
def get_corrections(PTIData_instance, 
                    ex_interp_method = 'cubic', em_interp_method = 'cubic', FS = False,
                    ex_split = 'none', em_split = 'none',
                    diode = True, excorr = True, emcorr = True,
                    const_diode=False):
    
    # Create a copy of the data instance
    corrections = numpy.ones(PTIData_instance.wavelengths.size)
    
    # Perform the diode correction made from the ExCorr RCQC signal
    if diode:
        if const_diode:
            corrections  /= numpy.mean(PTIData_instance.diode)
        else:
            corrections /= PTIData_instance.diode
    
    # Perform the LUT corrections
    if excorr:
        corrections /= load_excorr_file(PTIData_instance, ex_interp_method, ex_split)
    
    if emcorr:
        corrections *= load_emcorr_file(PTIData_instance, em_interp_method, FS)
    
    return corrections


''''''
def decorrect_cor_to_raw(PTIData = None,
                         ex_LUT_interpolation = 'cubic', em_LUT_interpolation = 'cubic', FS = False,
                         undo_diode = True, undo_ex_LUT = True, undo_em_LUT = True):
    
    data = copy.deepcopy(PTIData)

    corrections = get_corrections(data, 
                                  ex_LUT_interpolation, em_LUT_interpolation, FS,
                                  undo_diode, undo_ex_LUT, undo_em_LUT)

    data.raw_data = data.cor_data / corrections
    
    return data


''''''
def correct_raw_to_cor(PTIData = None, use_decorrected_as_raw = False,
                     baseline_fit_ranges = None, baseline_polynomial_degree = 1,
                     use_baseline_se = ('none', 'none'),
                     ex_LUT_interpolation = 'cubic', em_LUT_interpolation = 'cubic', FS = False,
                     ex_LUT_split =  'none', em_LUT_split = 'none',
                     undo_diode = True, undo_ex_LUT = True, undo_em_LUT = True,
                     apply_diode = True, apply_ex_LUT = True, apply_em_LUT = True,
                     const_diode=False):
    
    data = copy.deepcopy(PTIData)
    
    if use_decorrected_as_raw:
        # Changes the .raw_data member variable and sets it to the decorrected spectrum
        data = decorrect_cor_to_raw(data,
                                    ex_LUT_interpolation, em_LUT_interpolation, FS,
                                    undo_diode, undo_ex_LUT, undo_em_LUT)
        
    # Baseline subtract from the raw data
    baseline, params, errors = linear_baseline(PTIData = data,
                                               list_of_ranges = baseline_fit_ranges,
                                               use_incpt_se = use_baseline_se[0],
                                               use_slope_se=use_baseline_se[1])
    data.baseline = baseline
    data.baseline_incpt = params[0]
    data.baseline_slope = params[1]
    data.baseline_incpt_se = errors[0]
    data.baseline_slope_se = errors[1]

    
    data.raw_data = data.raw_data - baseline

    corrections = get_corrections(data,
                                  ex_LUT_interpolation, em_LUT_interpolation, FS,
                                  ex_LUT_split, em_LUT_split,
                                  apply_diode, apply_ex_LUT, apply_em_LUT,
                                  const_diode)

    data.cor_data = data.raw_data * corrections
    
    
    return data

    
