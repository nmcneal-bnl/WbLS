import matplotlib.pyplot as plt
import numpy as np

bisMSB_qy = [0.75, 0.76, 0.73, 0.76]
bisMSB_err = [0.02, 0.01, 0.02, 0.01]


ppo_etoh_0x31_qy =  [1.04, 0.68,0.74,0.68]
ppo_etoh_0x31_err = [0.15,0.03,0.02,0.01]

ppo_etoh_3x14_qy =  [0.66,0.60,0.64,0.46]
ppo_etoh_3x14_err = [0.09,0.06,0.04,0.03]

ppo_cx_qy =  [1.02,1.19,1.21]
ppo_cx_err = [.07,0.08,0.08]


fig=plt.figure(figsize=(20,8))
ax0 = fig.add_subplot(1,3,1)
ax1 = fig.add_subplot(1,3,2, sharey=ax0)
ax2 = fig.add_subplot(1,3,3, sharey=ax0)
fig.suptitle("Measured Quantum Yield Results for bisMSB and PPO", fontsize=20, weight='bold')

ax0.set_title("Measured Quantum Yields of 4.47 g/L bisMSB in LAB", fontsize=15)
ax0.plot([350, 360, 370, 380], bisMSB_qy, color='r', ls='none', marker='o', ms=3)
ax0.errorbar(x=[350, 360, 370, 380], y=bisMSB_qy, yerr=bisMSB_err, fmt="none", color='r', capsize=3)
ax0.set_ylabel("Measured QY", fontsize=12)
# ax0.axhline(y=0.926, color='k')
# ax0.axhline(y=0.926-0.053, color='k', ls='--')
# ax0.axhline(y=0.926+0.053, color='k', ls='--')

ax1.set_title("Measured Quantum Yields of PPO in Ethanol", fontsize=15)
ax1.plot([310, 320, 330, 340], ppo_etoh_0x31_qy, color='c', ls='none', marker='o', ms=3, label="0.31 g/L")
ax1.errorbar(x=[310, 320, 330, 340], y=ppo_etoh_0x31_qy, yerr=ppo_etoh_0x31_err, fmt="none", color='c', capsize=3)

ax1.plot([310, 320, 330, 340], ppo_etoh_3x14_qy, color='b', ls='none', marker='o', ms=3, label="3.14 g/L")
ax1.errorbar(x=[310, 320, 330, 340], y=ppo_etoh_3x14_qy, yerr=ppo_etoh_3x14_err, fmt="none", color='b', capsize=3)

ax2.set_title("Measured Quantum Yields of PPO in Cyclohexane", fontsize=15)
ax2.plot(313, ppo_cx_qy[0], color='r', ls='none', marker='o', ms=3, label = "0.85 mg/L")
ax2.errorbar(x=313, y=ppo_cx_qy[0], yerr=ppo_cx_err[0], fmt="none", color='r', capsize=3, label=None)
ax2.plot(313, ppo_cx_qy[1], color='g', ls='none', marker='o', ms=3, label = "95.14 mg/L")
ax2.errorbar(x=313, y=ppo_cx_qy[1], yerr=ppo_cx_err[1], fmt="none", color='g', capsize=3, label=None)
ax2.plot(313, ppo_cx_qy[2], color='b', ls='none', marker='o', ms=3, label = "0.95 g/L")
ax2.errorbar(x=313, y=ppo_cx_qy[2], yerr=ppo_cx_err[2], fmt="none", color='b', capsize=3, label=None)

for ax in fig.get_axes():
    ax.grid()
    ax.legend(fontsize=12)
    ax.set_yticks(np.arange(0.4, 1.4, .1))
    ax.set_xlabel("Wavelength (nm)", fontsize=12)
plt.show()

