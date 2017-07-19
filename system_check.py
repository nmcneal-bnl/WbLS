import matplotlib.pyplot as plt
from PTI.ReadDataFiles import PTIData

import glob, os

all_paths = list()
for root, dirs, files in os.walk("Noah/PTI System Check"):
    for f in files:
        fullpath = os.path.join(root, f)
        if os.path.splitext(fullpath)[1] == '.txt':
            all_paths.append(fullpath)

for path in all_paths:
    try:
        data = PTIData(path)

        fig = data.plot()

        fname = path.split('/')[-1]
        plt.savefig('Noah/PTI System Check/Plots/' + fname + '.png')
        plt.close('all')

    except KeyboardInterrupt:
        break
    except:
        plt.close('all')
