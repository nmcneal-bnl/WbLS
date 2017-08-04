import matplotlib.pyplot as plt
from PTI.ReadDataFiles import PTIData
import os 

all_paths = list()
for root, dirs, files in os.walk("Henry"):
    for f in files:
        fullpath = os.path.join(root, f)
        if os.path.splitext(fullpath)[1] == '.txt':
            all_paths.append(fullpath)

for path in all_paths:
    new_path = path.split('/')[1:]
    new_path = 'All_Henry_Plots/'+'/'.join(new_path)
    new_path = new_path[:-3] + 'png'
    
    if not os.path.exists(os.path.dirname(new_path)):
        try:
            os.makedirs(os.path.dirname(new_path))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    print new_path
    try:
        data = PTIData(path)
        fig = data.plot()
        fig.get_axes()[-1].axvline(x=data.ex_range[0], color = 'r', ls ='--')
        fig.get_axes()[0].axvline(x=data.ex_range[0], color = 'r', ls ='--')
        fname = path.split('/')[-1].split('.')[0]
        plt.savefig(new_path)
        plt.close('all')
    except KeyboardInterrupt:
        break
    except:
        pass


    
