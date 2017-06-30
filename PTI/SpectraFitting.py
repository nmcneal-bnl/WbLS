#!/usr/bin/env python2

import numpy
from scipy.optimize import curve_fit

''' Uses a least-squared routine to fit the spectrum to a purely linear function
    The function may be a constant or be a first degree polynomial'''
def linear_fit(PTIData, poly_degree, list_of_ranges=[]):

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

'''Calculates a baseline from the linear fit and returns the baseline-subtracted data'''
def get_linear_subtracted(PTIData, poly_deg, list_of_ranges=[]):
    slope, intercept = (0,0)
    slope, intercept = linear_fit(PTIData, poly_deg, list_of_ranges)

    return PTIData.cor_data - (slope*PTIData.wavelengths + intercept)

'''Function to take in an array and five parameters. Returns the Gaussian-Linear
   values as an array'''
def gaussian_func(x, params):
    s_x = params[0] * numpy.exp(-params[1]*numpy.square(x-params[2]))
    return s_x


'''Fits a spectrum to '''
def gaussian_fit(PTIData, wavelength_ranges, param_guess):
    def gaussian_linear(x,a,b,c):
        s_x = a * numpy.exp(-b*numpy.square(x-c))
        return s_x

    x_data = numpy.array([])
    y_data = numpy.array([])
    for arange in wavelength_ranges:
        where = numpy.where((PTIData.wavelengths >= arange[0]) &
                            (PTIData.wavelengths <= arange[1]))
        x_data = numpy.append(x_data, PTIData.wavelengths[where])
        y_data = numpy.append(y_data, PTIData.cor_data[where])
        
    params = curve_fit(gaussian_linear, x_data, y_data, param_guess)
    return params[0]


'''Function to take in an array and five parameters. Returns the Gaussian-Linear
   values as an array'''
def gaussian_linear_func(x, params):
    s_x = params[0] * numpy.exp(-params[1]*numpy.square(x-params[2])) + \
          params[3]*x + params[4]
    return s_x


'''Fits a spectrum to '''
def gaussian_linear_fit(PTIData, wavelength_ranges, param_guess):
    def gaussian_linear(x,a,b,c,d,f):
        s_x = a * numpy.exp(-b*numpy.square(x-c))  + d*x + f
        return s_x

    x_data = numpy.array([])
    y_data = numpy.array([])
    for arange in wavelength_ranges:
        where = numpy.where((PTIData.wavelengths >= arange[0]) &
                            (PTIData.wavelengths <= arange[1]))
        x_data = numpy.append(x_data, PTIData.wavelengths[where])
        y_data = numpy.append(y_data, PTIData.cor_data[where])
        
    params = curve_fit(gaussian_linear, x_data, y_data, param_guess)
    return params[0]
