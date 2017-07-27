import numpy as np
import matplotlib.pyplot as plt

def monochromator(wavelength_setting, set_of_wavelengths, offset = 0):
    if wavelength_setting in set_of_wavelengths or wavelength_setting/2 in set_of_wavelengths:
        tmp = np.array([wavelength_setting + offset])
        return tmp.astype(int)
    else:
        return np.array([0])

xe_arc_lamp = range(250000, 650000)  # Picometers
ex_offset = 1000 * 0
em_offset = 1000 * 0

ex_monochromator_output = monochromator(300000, xe_arc_lamp, ex_offset)
for x in range(300000, 650000, 500):
    em_monochromator_output = monochromator(x, ex_monochromator_output, em_offset)
    plt.axvline(x=0.22058956)
    plt.plot(x, em_monochromator_output[0], 'bo')
plt.show()


