import numpy
from scipy.integrate import simps

def integrate_between(PTIData1, PTIData2, int_range):
    difference = PTIData1.cor_data - PTIData2.cor_data

    # Integration
    limits = numpy.where((PTIData1.wavelengths >= int_range[0]) & 
                         (PTIData1.wavelengths <= int_range[1]))

    def_int = simps(y = difference[limits],
                    dx = PTIData1.step_size)
                                       
    return def_int

def calc_QY_PTI(pti_blank_data, pti_fluor_data, 
                ex_integration_range, em_integration_range):
    
    # Get difference
    difference = pti_blank_data.cor_data -  pti_fluor_data.cor_data    

    # Integration
    ex_limits = numpy.where((pti_blank_data.wavelengths >= ex_integration_range[0]) & 
                            (pti_blank_data.wavelengths <= ex_integration_range[1]))

    em_limits = numpy.where((pti_blank_data.wavelengths >= em_integration_range[0]) & 
                            (pti_blank_data.wavelengths <= em_integration_range[1]))

    num_emitted  = simps(y = -difference[em_limits],
                        dx = pti_blank_data.step_size)
                                       
    num_absorbed = simps(y = difference[ex_limits],
                        dx = pti_blank_data.step_size)
                                       
    return num_emitted/num_absorbed

def calc_QY_fit(wavelengths, fit_blank_data, fit_fluor_data, 
                ex_integration_range, em_integration_range):
    
    # Get diffrence
    difference = fit_blank_data - fit_fluor_data
    
    # Integration
    step_size = wavelengths[1] - wavelengths[0]

    ex_limits = numpy.where((wavelengths >= ex_integration_range[0]) & 
                            (wavelengths <= ex_integration_range[1]))
    
    em_limits = numpy.where((wavelengths >= em_integration_range[0]) & 
                            (wavelengths <= em_integration_range[1]))

    num_emitted =  simps(y = -difference[em_limits],
                        dx = step_size)

    num_absorbed = simps(y = difference[ex_limits],
                        dx = step_size)
    return num_emitted/num_absorbed
