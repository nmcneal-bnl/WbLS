import copy
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import numpy
import matplotlib.pyplot as plt


def linear_func(x, b, m):
    return m * x + b


def linear_baseline_params(PTIData, list_of_ranges):
    """ Uses a least-squared routine to fit  data to a linear function.
        The fit is performed over the given wavelength ranges.
        Returns a tuple containing: ("""

    X = numpy.array([])
    Y = numpy.array([])
    for arange in list_of_ranges:
        start = arange[0]
        end = arange[1]

        X = numpy.append(X, numpy.arange(start, end+PTIData.step_size, PTIData.step_size))

        select_by_wavelength = numpy.where((PTIData.wavelengths >= start) &
                                           (PTIData.wavelengths <= end))
        Y = numpy.append(Y, PTIData.raw_data[select_by_wavelength])

    params, cov_matrix = curve_fit(linear_func, X, Y)

    return(params, cov_matrix)


def linear_baseline(PTIData, list_of_ranges,
                    use_incpt_se = 'none', use_slope_se = 'none'):
    """Uses a least-squared routine to fit  data to a polynomial of a given
        degree. The fit is performed over the given wavelength ranges.
        Returns an array with the polynomial evaluated over the full spectrum range
        with the same step size."""

    # Calculate the parameters for the fit and their standard errors. Order: (intercept, slope)
    fit_params, cov_matrix = linear_baseline_params(PTIData, list_of_ranges)
    fit_params = list(fit_params)
    # print cov_matrix
    errors = [numpy.sqrt(cov_matrix[0][0]), numpy.sqrt(cov_matrix[1][1])]

    # Use the standard errors to modify either parameter
    incpt ={'none': fit_params[0],
            'minus': fit_params[0] - errors[0],
            'plus': fit_params[0] + errors[0]}

    slope ={'none': fit_params[1],
            'minus': fit_params[1] - errors[1],
            'plus': fit_params[1] + errors[1]}

    baseline = linear_func(x=PTIData.wavelengths, m=slope[use_slope_se],b=incpt[use_incpt_se])

    # Return the baseline array as well as the fit parameters and errors.
    return baseline, fit_params, errors


def gaussian_func(x, a, b, c, d):
    return a * numpy.exp((-(x - b) ** 2) / (2 * c)) + d


def gaussian_fit(x_data, y_data, guess=(1, 1, 1, 0)):
    params, cov_matrix = curve_fit(gaussian_func, x_data, y_data, p0=guess)
    return params, cov_matrix


def get_true_excitation_wavelength_from_secondary_peak(PTIData, dx_around_peak = 5):
    global_peak = numpy.max(PTIData.raw_data)
    peak_wavelength = PTIData.wavelengths[numpy.where(PTIData.raw_data == global_peak)][0]

    peak_range = numpy.where((PTIData.wavelengths >= peak_wavelength - dx_around_peak) &
                             (PTIData.wavelengths <= peak_wavelength + dx_around_peak))
    peak_range2 = numpy.where((PTIData.wavelengths >= 2*peak_wavelength - dx_around_peak) &
                              (PTIData.wavelengths <= 2*peak_wavelength + dx_around_peak))

    x_data = PTIData.wavelengths[peak_range]
    y_data = PTIData.raw_data[peak_range]
    x_data2 = PTIData.wavelengths[peak_range2]
    y_data2 = PTIData.raw_data[peak_range2]

    secondary_peak = numpy.max(y_data2)

    guess = (global_peak, peak_wavelength, 2*dx_around_peak/2.35482, 0)
    guess2 = (secondary_peak, 2*peak_wavelength, 2*dx_around_peak/2.35482, 0)

    params, cov_matrix = gaussian_fit(x_data, y_data, guess)
    params2, cov_matrix2 = gaussian_fit(x_data2, y_data2, guess2)

    # x = numpy.arange(peak_wavelength - dx_around_peak, peak_wavelength + dx_around_peak, 0.1)
    # x2 = numpy.arange(2*peak_wavelength - dx_around_peak, 2*peak_wavelength + dx_around_peak, 0.1)
    # y = gaussian_func(x, params[0], params[1], params[2], params[3])
    # y2 = gaussian_func(x2, params2[0], params2[1], params2[2], params2[3])
    #
    # plt.plot(PTIData.wavelengths, PTIData.raw_data)
    # plt.plot(x, y)
    # plt.plot(x2, y2)
    # plt.show()
    return params2[1] - params[1]


def get_excitation_monochromator_offset(PTIData, dx_around_peak = 5):
    theoretical_excitation= PTIData.ex_range[0]

    actual_excitation = get_true_excitation_wavelength_from_secondary_peak(PTIData=PTIData,
                                                                           dx_around_peak=dx_around_peak)
    # print actual_excitation, theoretical_excitation,
    offset = actual_excitation - theoretical_excitation

    return offset


def get_emission_monochromator_shift(PTIData, dx_around_peak = 5):
    global_peak = numpy.max(PTIData.raw_data)
    peak_wavelength = PTIData.wavelengths[numpy.where(PTIData.raw_data == global_peak)][0]
    peak_range = numpy.where((PTIData.wavelengths >= peak_wavelength - dx_around_peak) &
                             (PTIData.wavelengths <= peak_wavelength + dx_around_peak))
    guess = (global_peak, peak_wavelength, 2 * dx_around_peak / 2.35482, 0)

    gaussian_params, _ = gaussian_fit(x_data=PTIData.wavelengths[peak_range],
                                      y_data=PTIData.raw_data[peak_range],
                                      guess=guess)

    actual_excitation = get_true_excitation_wavelength_from_secondary_peak(PTIData=PTIData,
                                                                           dx_around_peak=dx_around_peak)
    # print gaussian_params[1]
    offset = actual_excitation - gaussian_params[1]

    return offset


def load_excorr_file(PTIData_instance, interp_method = 'cubic', split = 'none', shift = 0):
    
    excorr = numpy.genfromtxt('PTI/correction_data/excorr.txt',
                           skip_header = 6,
                           skip_footer = 1,
                           usecols = 1)
    LUT_start = 250
    LUT_end = 750
    LUT_step = 1

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
    excorr_range = numpy.arange(LUT_start, LUT_end + LUT_step, LUT_step)
    ex_wavelength = PTIData_instance.ex_range[0]

    left_fill_value = excorr[0]
    right_fill_value = excorr[-1]
    fill_value = (left_fill_value, right_fill_value)

    excorr = interp1d(x=excorr_range,
                      y=excorr,
                      kind=interp_method,
                      bounds_error=False,
                      fill_value=fill_value)(ex_wavelength + shift)

    return excorr


def load_emcorr_file(PTIData_instance, interp_method = 'cubic', FS = False, split = 'none', shift = 0):

    if FS:
        fname = 'PTI/correction_data/emcorri.txt'
        LUT_start = 250
        LUT_end = 850

    if not FS:
        fname = 'PTI/correction_data/emcorr-sphere-quanta.txt'
        LUT_start = 300
        LUT_end = 848
    LUT_step = 2

    emcorr = numpy.genfromtxt(fname,
                              skip_header = 6,
                              skip_footer = 1,
                              usecols = 1)

    if split.lower() == 'none':
        pass
    elif split.lower() == 'even':
        emcorr = emcorr[::2]
        LUT_step *= 2
    elif split.lower() == 'odd':
        emcorr = emcorr[1::2]
        LUT_step *= 2
        LUT_start += 2
        LUT_end -= 2
    else:
        print "ERROR: Not a valid method for splitting LUT"
        return None

    step = PTIData_instance.step_size
    emcorr_wavelengths = numpy.arange(LUT_start, LUT_end + LUT_step, LUT_step)

    xvals = numpy.arange(LUT_start, LUT_end + step, step)
    emcorr = interp1d(x=emcorr_wavelengths,
                      y=emcorr,
                      kind=interp_method,
                      fill_value='extrapolate')(xvals + shift)

    min_data_wavelength = PTIData_instance.wavelengths[0]
    max_data_wavelength = PTIData_instance.wavelengths[-1]
    
    needed_wavelengths = numpy.where((xvals >= min_data_wavelength) & 
                                     (xvals <= max_data_wavelength))

    emcorr = emcorr[needed_wavelengths]
    
    if min_data_wavelength < LUT_start:
        extra_x = numpy.arange(min_data_wavelength, LUT_start, step)
        left_interp = emcorr[0]*numpy.ones(extra_x.size)
        emcorr = numpy.append(left_interp, emcorr)
    if max_data_wavelength > LUT_end:
        extra_x = numpy.arange(LUT_end, max_data_wavelength, step)
        right_interp = emcorr[-1]*numpy.ones(extra_x.size)
        emcorr = numpy.append(emcorr, right_interp)
    
    return emcorr
    

def get_corrections(PTIData_instance,
                    ex_interp_method = 'cubic', em_interp_method = 'cubic', FS = False,
                    ex_split = 'none', em_split = 'none',
                    ex_shift = 0, em_shift = 0,
                    diode = True, excorr = True, emcorr = True,
                    const_diode=False):
    
    # Create a copy of the data instance
    corrections = numpy.ones(PTIData_instance.wavelengths.size)
    
    # Perform the diode correction made from the ExCorr RCQC signal
    if diode:
        if const_diode:
            corrections /= numpy.mean(PTIData_instance.diode)
        else:
            corrections /= PTIData_instance.diode
    
    # Perform the LUT corrections
    if excorr:
        corrections /= load_excorr_file(PTIData_instance=PTIData_instance,
                                        interp_method=ex_interp_method,
                                        split=ex_split,
                                        shift=ex_shift)
    
    if emcorr:
        corrections *= load_emcorr_file(PTIData_instance=PTIData_instance,
                                        interp_method=em_interp_method,
                                        FS=FS,
                                        split=em_split,
                                        shift=em_shift)
    
    return corrections


def decorrect_cor_to_raw(PTIData = None,
                         ex_LUT_interpolation = 'cubic', em_LUT_interpolation = 'cubic', FS = False,
                         undo_diode = True, undo_ex_LUT = True, undo_em_LUT = True):
    
    data = copy.deepcopy(PTIData)

    corrections = get_corrections(PTIData_instance=data,
                                  ex_interp_method=ex_LUT_interpolation, em_interp_method=em_LUT_interpolation,
                                  FS=FS,
                                  diode=undo_diode, excorr=undo_ex_LUT, emcorr=undo_em_LUT)

    data.raw_data = data.cor_data / corrections
    
    return data


def correct_raw_to_cor(PTIData = None, use_decorrected_as_raw = False,
                       baseline_fit_ranges = None, baseline_polynomial_degree = 1,
                       use_baseline_se = ('none', 'none'), gaussian_fit_dx_around_peak = 10,
                       ex_LUT_split='none', em_LUT_split='none',
                       ex_LUT_interpolation = 'cubic', em_LUT_interpolation = 'cubic', FS = False,
                       shift_LUT = False, ex_shift=0, em_shift=0,
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

    if not shift_LUT:
        ex_shift = 0
        em_shift = 0
    else:
        if 2 * PTIData.ex_range[0] < PTIData.em_range[1]:
            ex_shift = get_excitation_monochromator_offset(PTIData, dx_around_peak = 5)
            em_shift = get_emission_monochromator_shift(PTIData, dx_around_peak = 5)
        else:
            pass
    data.ex_monochromator_offset = ex_shift
    data.em_monochromator_offset = em_shift

    corrections = get_corrections(PTIData_instance=data,
                                  ex_interp_method=ex_LUT_interpolation, em_interp_method=em_LUT_interpolation, FS=FS,
                                  ex_split=ex_LUT_split, em_split=em_LUT_split,
                                  ex_shift=ex_shift, em_shift=em_shift,
                                  diode=apply_diode, excorr=apply_ex_LUT, emcorr=apply_em_LUT,
                                  const_diode=const_diode)



    data.cor_data = data.raw_data * corrections
    data.raw_data += baseline
    
    return data
