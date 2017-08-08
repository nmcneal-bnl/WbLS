import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

files = ['baseline_options.txt', 'const_diode.txt', 'LUT_interpolation.txt', 'LUT_splitting.txt', 'LUT_shifting.txt',
         'correction_region_start.txt', 'all_options.txt']

# The number of bins to use for the histograms for the distributions caused by each parameter
nums_bins = [4, 2, 4, 3, 2, 5, 50]

# The names of each of the parameters
parameters=["Baseline error", "Constant Diode Correction", "LUT Interpolation Method",
            "LUT Even-odd Splitting", "LUT Offset Shifting", "Start of Correction Region", "All Variations"]

# Create the files to hold the statistics for each wavelength
ex_wavelengths = [str(wavelength) + ' nm' for wavelength in range(310, 350, 10)]
for j in range(len(ex_wavelengths)):
    data_file = open("Summaries/" + ex_wavelengths[j] + "_parameter_summary_statistics.txt", 'w+')
    data_file.write("Parameter,Min,Max,Mean,Std,Half Range\n")
    data_file.close()

for i in range(len(files)):
    '''Calculating the statistics for a given parameter'''
    # File I/O
    f = files[i]
    data = pd.read_csv(f, header=0)

    # Run the calculations and store the results for each wavelength in a list
    means = [np.mean(data[col]) for col in ex_wavelengths]
    stdevs = [np.std(data[col]) for col in ex_wavelengths]
    maxs = np.array([np.max(data[col]) for col in ex_wavelengths])
    mins = np.array([np.min(data[col]) for col in ex_wavelengths])
    half_ranges = 0.5*(maxs - mins)

    if f == "all_options.txt":
        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_subplot(111)
        ax.errorbar([310, 320, 330, 340], means, yerr=half_ranges, fmt='o')
        ax.set_xscale("log")
        plt.savefig("QYs.png")

    '''Write the parameter results into each of the four files, one for each excitation wavelength'''
    for j in range(len(ex_wavelengths)):
        data_file = open("Summaries/" +ex_wavelengths[j] + "_parameter_summary_statistics.txt", 'a')
        data_file.write(parameters[i] + ',' + "%f,%f,%f,%f,%f\n" %(mins[j], maxs[j], means[j], stdevs[j], half_ranges[j]))
        data_file.close()

    '''Plot the data in a histogram'''
    # Set what error to use
    error = half_ranges

    # Create a histogram for each concentration
    num_bins = nums_bins[i]

    fig = plt.figure(figsize=(20, 10))
    fig.suptitle(parameters[i]+"\nQuantum Yield for 3.14 g/L PPO in EtOH", fontsize=20)
    for k in range(4):
        ax = fig.add_subplot(1,4,k+1)
        n, bins,_ = ax.hist(data[ex_wavelengths[k]], bins=num_bins, color="#6dd5ff",label='_nolegend_')
        step = bins[1]-bins[0]
        ax.set_title(ex_wavelengths[k] +"\n"+r'$\mu$=%0.3f' % means[k] +
                     "\n"+r'$\epsilon$=%0.3f (%0.3f%%)' %(error[k], 100*error[k]/means[k]))
        ax.axvline(x=means[k], color='r', ls='--', label=r'$\mu$')

        error_y = 0.8 * np.max(n)
        error_percent = 0.05

        ax.vlines(x=means[k] + error[k], color='r', linestyles='-',
                  ymin=(1 - error_percent) * error_y, ymax=(1 + error_percent) * error_y)
        ax.vlines(x=means[k] - error[k], color='r', linestyles='-',
                  ymin=(1 - error_percent) * error_y, ymax=(1 + error_percent) * error_y,
                  label = r'$\epsilon$')

        ax.annotate("", xytext=(means[k], error_y),
                    xy=(means[k] + error[k], error_y),
                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle="arc3",
                                    color='r')
                    )
        ax.annotate("", xytext=(means[k], error_y),
                    xy=(means[k] - error[k], error_y),
                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle="arc3",
                                    color='r')
                    )
        ax.legend()

    plt.subplots_adjust(top=.8)
    plt.savefig("Distributions/%s.svg" % parameters[i], format='svg', dpi=300)
    plt.close('all')

# Make the table of relative uncertainties
data = [pd.read_csv("Summaries/" + wavelength + "_parameter_summary_statistics.txt") for wavelength in ex_wavelengths]
table_file = open("relative_error_contributions.txt", 'w+')
table_file.write("Parameter, 310 nm, 320 nm, 330  nm, 340 nm\n")
for i in range(7):
    table_file.write("%s," % data[0]['Parameter'][i])
    for j in range(0,4):
        table_file.write("%f," % (100*data[j]['Half Range'][i] / data[j]['Mean'][i]))
    table_file.write('\n')

table_file.close()

# Plot the change in the QY as we iterate through the option
data = pd.read_csv("all_options.txt")
fig = plt.figure(figsize=(15,5))
fig.suptitle("Quantum Yield for 0.31 g/L PPO in EtOH", fontsize=20)

for k in range(4):
    ax = fig.add_subplot(2,4,k+1)
    ax2 = fig.add_subplot(2, 4, k+5)
    ax.plot(data[ex_wavelengths[k]], '+')

    ax.set_title(ex_wavelengths[k])
    ax2.set_title("Difference")


    ax2.plot(np.diff(data[ex_wavelengths[k]]), 'r+')
plt.xlabel("Variation Number")
plt.subplots_adjust(top=.9)
plt.savefig("Distributions/variations.png", format='png', dpi=300)
plt.close('all')



