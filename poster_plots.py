import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy
mpl.rc('font',family='Times New Roman')

import PTI.Corrections as PTICorr
from PTI.ReadDataFiles import PTIData

test = "Henry/Sphere/PPO_ETOH/EmissionScan_3x14gperL_PPOinETOH_ex340_2sec_160831.txt"

'''Showing how the baseline and corrections affect data'''
data1 = PTIData(test)
data1_baseline,_,_ = PTICorr.linear_baseline(data1, list_of_ranges=[[300,325], [450,650]])

fig1 = plt.figure(figsize=(15,10))
ax1 = plt.subplot2grid((2,2), (0,0))
ax2 = plt.subplot2grid((2,2), (0,1))
ax3 = plt.subplot2grid((2,2), (1,0), colspan=2)

fig1.suptitle("Raw and Corrected Emission Spectra of PPO in Ethanol Excited at 340 nm\n",
              weight='bold', y=0.95, fontsize=20)
ax1.set_title("Initial Raw Emission Spectrum", fontsize = 12)
ax2.set_title("Total Correction Factor", fontsize = 12)
ax3.set_title("Final Corrected Emission Spectrum", fontsize = 12)


ax1.plot(data1.wavelengths, data1.raw_data/1000,"#d60027")
ax1.set_ylabel("Emission Intensity [a.u.]", fontsize = 12)
ax1.axhline(y=0, xmin=0, xmax=1, linewidth=1, color = 'k')
ax1.legend(["Raw Data"], fontsize = 12)
ax1.xaxis.set_ticks(range(300, 510, 100))
ax1.yaxis.set_ticks(range(0, 7, 3))

corrections = PTICorr.get_corrections(data1)
ax2.plot(data1.wavelengths, corrections/1000,"#268a47")
ax2.set_ylabel("Magnitude [a.u.]", fontsize = 12)
ax2.xaxis.set_ticks(range(300, 510, 100))
ax2.yaxis.set_ticks(numpy.arange(0, 0.07, .03))


data1 = PTICorr.correct_raw_to_cor(PTIData=data1, baseline_fit_ranges=[[300,335], [450,650]])
ax3.axhline(y=0, xmin=0, xmax=1, linewidth=1, color = 'k')
ax3.plot(data1.wavelengths, data1.cor_data/1000, "#3857ff")
ax3.set_ylabel("Emission Intensity [a.u.]", fontsize = 12)
ax3.xaxis.set_ticks(range(300, 510, 100))
ax3.yaxis.set_ticks(range(0, 150, 50))


fig1.text(0.5, 0.04, "Wavelength [nm]", weight='semibold', ha='center', fontsize = 15)

for ax in [ax1,ax2,ax3]:
    ax.set_xlim(300,500)
    ax.grid()


plt.savefig('test.svg', format='svg', dpi=1200)

'''Showing an example of the quantum yield calculation'''
'''
fig2 = plt.figure(figsize = (16,10))
ax1 = plt.subplot2grid((2,2), (0,0), colspan = 2)
ax2 = plt.subplot2grid((2,2), (1,0))
ax3 = plt.subplot2grid((2,2), (1,1))

fig2.suptitle("Corrected Emission Spectra of LAB and bisMSB in LAB",
              weight='semibold', y=0.95, fontsize=20)
#ax1.set_title("Comparison of Fully Corrected Spectra")
ax2.set_title("Absorption Wavelength Range")
ax3.set_title("Emission Wavelength Range")

blank = PTIData("Henry/Sphere/bisMSB_LAB/EmissionScan_LAB_ex380_2sec_160823.txt")
fluor = PTIData("Henry/Sphere/bisMSB_LAB/EmissionScan_bisMSBinLAB_4.47mgL_ex380_2sec_160824.txt")

blank = PTICorr.correct_raw_to_cor(blank, baseline_fit_ranges=[[300, 325], [500, 650]])
fluor = PTICorr.correct_raw_to_cor(fluor, baseline_fit_ranges=[[300, 325], [500, 650]])

blank.cor_data /= 1000
fluor.cor_data /= 1000

ex_int_range = [375, 385]
em_int_range = [365, 650]
ex_limits_bool = (blank.wavelengths >= ex_int_range[0]) & (blank.wavelengths <= ex_int_range[1])
em_limits_bool = (blank.wavelengths >= em_int_range[0]) & (blank.wavelengths <= em_int_range[1])

ax1.plot(blank.wavelengths, blank.cor_data,"#636266")
ax1.plot(fluor.wavelengths, fluor.cor_data,"#dc143c")
ax1.legend(["Only LAB", "bisMSB in LAB"],  fontsize = 12)
ax1.set_xlim(350, 550)
ax1.xaxis.set_ticks(range(350, 560, 50))
ax1.yaxis.set_ticks(range(0, 350, 100))


ax2.plot(blank.wavelengths, blank.cor_data,"#636266",label='_nolegend_')
ax2.plot(fluor.wavelengths, fluor.cor_data,"#dc143c",label='_nolegend_')
ax2.fill_between(blank.wavelengths, blank.cor_data, fluor.cor_data,  ex_limits_bool, 
                 color="#97f9ca", alpha=0.8)
ax2.legend(["Photons Absorbed"], fontsize = 12)
ax2.set_xlim(370, 390)
ax2.xaxis.set_ticks(range(370, 400, 10))
ax2.yaxis.set_ticks(range(0, 350, 100))

ax3.plot(blank.wavelengths, blank.cor_data,"#636266",label='_nolegend_')
ax3.plot(fluor.wavelengths, fluor.cor_data,"#dc143c",label='_nolegend_')
ax3.fill_between(blank.wavelengths, blank.cor_data, fluor.cor_data,  em_limits_bool, 
                 color="#97f9ca", alpha=0.8)
ax3.legend(["Photons Emitted"], fontsize = 12)
ax3.set_xlim(370, 550)
ax3.set_ylim(-0.5, 4.5)
ax3.xaxis.set_ticks(range(350, 560, 100))
ax3.yaxis.set_ticks(range(0, 5, 2))


fig2.text(0.5, 0.04, "Wavelength [nm]", weight='semibold', ha='center', fontsize = 15)
fig2.text(0.08, 0.5, "Emission Intensity [a.u.]", va='center', 
          weight='semibold',rotation='vertical', fontsize = 15)



plt.show()
'''



