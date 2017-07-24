import copy
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps

def integrate_between(blank, fluor, int_range):
    difference = blank.cor_data - fluor.cor_data

    limits = np.where((blank.wavelengths >= int_range[0]) &
                      (blank.wavelengths <= int_range[1]))

    int = simps(y = difference[limits], dx = blank.step_size)

    return int

def calculate_quantum_yield(blank, fluor, ex_int_range, em_int_range):

    num_absorbed = integrate_between(blank, fluor, ex_int_range)

    num_emitted = integrate_between(fluor, blank, em_int_range)

    return num_emitted/num_absorbed

def calculate_QY(blank, fluor, ex_int_range, em_int_range):
    return calculate_quantum_yield(blank, fluor, ex_int_range, em_int_range)

