import matplotlib.pyplot as plt
from PTI.ReadDataFiles import PTIData

import glob, os

all_paths = list()
for root, dirs, files in os.walk("Henry"):
    for f in files:
        fullpath = os.path.join(root, f)
        if os.path.splitext(fullpath)[1] == '.txt':
            all_paths.append(fullpath)

for path in all_paths:
    if 'box' in path.lower():      
        try:

            data = PTIData(path)

            fig = data.plot()
            fig.get_axes()[-1].axvline(x=data.ex_range[0], color = 'r', ls ='--')
            fname = path.split('/')[-1]
            plt.savefig('All_Henry_BOX_Plots/' + fname + '.png')
            plt.close('all')
        except KeyboardInterrupt:
            break
        except:
            plt.close('all')

