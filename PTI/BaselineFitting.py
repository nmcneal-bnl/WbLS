#!/usr/bin/env python2

import numpy
from scipy.optimize import curve_fit

''' Uses a least-squared routine to fit the spectrum to a purely linear function
    The function may be a constant or be a first degree polynomial'''
def polynomial_baseline_params(PTIData, poly_degree, list_of_ranges):

    x = numpy.array([])
    y = numpy.array([])
    for arange in list_of_ranges:
        start = arange[0]
        end = arange[1]
        x = numpy.append(x, numpy.arange(start, end+PTIData.step_size, PTIData.step_size))
        
        select_by_wavelength = numpy.where((PTIData.wavelengths >= start) &
                                           (PTIData.wavelengths <= end))

        y = numpy.append(y, PTIData.raw_data[select_by_wavelength])

    fit_params =  numpy.polyfit(x, y, deg=poly_degree)
    
    return fit_params

def polynomial_baseline(PTIData, poly_degree, list_of_ranges):

    fit_params = polynomial_baseline_params(PTIData, poly_degree, list_of_ranges)

    baseline = numpy.zeros(PTIData.wavelengths.size)
    
    for i in range(poly_degree + 1):
        baseline += fit_params[i] * PTIData.wavelengths**(poly_degree - i)

    return baseline

