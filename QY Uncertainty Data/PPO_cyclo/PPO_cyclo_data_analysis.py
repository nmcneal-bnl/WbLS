import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("all_options.txt", header=0)
concentrations = ['0.04 mM','0.43 mM','4.3 mM']

# Gather the statistics
means = [np.mean(data[col]) for col in concentrations]
stdevs = [np.std(data[col]) for col in concentrations]
maxs = np.array([np.max(data[col]) for col in concentrations])
mins = np.array([np.min(data[col]) for col in concentrations])
ranges = 0.5*(maxs - mins)

#Choose what to use for the error
error = ranges

# Create a histogram for each concentration
num_bins = data.shape[0] / 500

fig = plt.figure(figsize=(16,5))
for i in range(3):
    ax = fig.add_subplot(1,3,i+1)
    n, bins,_ = ax.hist(data[concentrations[i]], bins=num_bins, color="#6dd5ff")
    step = bins[1]-bins[0]
    ax.set_title(concentrations[i])

    ax.axvline(x=means[i], color='r', ls='--')
    ax.text(means[i]- 5*step, 0.9*np.max(n), r'$\mu$', color='r')

    error_y = 0.8 * np.max(n)
    error_percent = 0.1

    ax.vlines(x=means[i] + error[i], color='r', linestyles='--',
              ymin=(1 - error_percent) * error_y, ymax=(1 + error_percent) * error_y)
    ax.vlines(x=means[i] - error[i], color='r', linestyles='--',
              ymin=(1 - error_percent) * error_y, ymax=(1 + error_percent) * error_y)

    ax.annotate("", xytext=(means[i], error_y),
                xy=(means[i] + error[i], error_y),
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3",
                                color='r')
                )
    ax.annotate("", xytext=(means[i], error_y),
                xy=(means[i] - error[i], error_y),
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3",
                                color='r')
                )

plt.show()





