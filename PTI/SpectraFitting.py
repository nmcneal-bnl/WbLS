#!/usr/bin/env python2

import numpy
from scipy.optimize import curve_fit

''' Uses a least-squared routine to fit the spectrum to a purely linear function
    The function may be a constant or be a first degree polynomial'''
def linear_baseline(PTIData, poly_degree, list_of_ranges=[]):

    x = numpy.array([])
    y = numpy.array([])
    for range in list_of_ranges:
        start = range[0]
        end = range[1]
        x = numpy.append(x, numpy.arange(start, end+PTIData.step_size, PTIData.step_size))
        
        select_by_wavelength = numpy.where((PTIData.wavelengths >= start) &
                                           (PTIData.wavelengths <= end))

        y = numpy.append(y, PTIData.cor_data[select_by_wavelength])
        fit_params =  numpy.polyfit(x, y, deg=poly_degree)
        
    slope,intercept = (0,0)
    if poly_degree == 0:
        slope = 0
        intercept = fit_params[0]  
    else:
        slope, intercept = fit_params

    return slope, intercept

