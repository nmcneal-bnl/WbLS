from ReadDataFiles import PTIData

test = PTIData("/home/nmcneal/Documents/WbLS/Noah/Integrating Sphere Tests for PPO Contamination/EmScan_IS_3x14gperL_PPOinETOH_ex310_em300-650_2sec_20170630_1147.txt")
print len(test.excorr)
print len(test.cor_data)

