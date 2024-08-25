import numpy as np
from numba import njit, prange
import h5py
import os
from tqdm import trange
import socket

@njit(parallel=True)
def window_baseline(fluoTrace, min_w, p=10):
    numFrames=len(fluoTrace)
    smoothBaseline=np.zeros_like(fluoTrace)
    for i in prange(numFrames):
        idx1 = max(0, i - min_w)
        idx2 = min(numFrames-1, i + min_w)
        smoothBaseline[i]=np.nanpercentile(fluoTrace[idx1:idx2],p)
    return smoothBaseline


import argparse
 
parser = argparse.ArgumentParser(description='baseline correction for fluorescent traces of merged cells.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("experiment_filename", help="experiment_filename")
parser.add_argument("server", help="server")
parser.add_argument("experimenter", help="experimenter")
parser.add_argument("gain", type=int, help="gain")
# parser.add_argument("analyzer", help="analyzer")
parser.add_argument("-a", "--analyzer", help="analyzer", nargs='?', const="chuyu", default="chuyu")
args = parser.parse_args()
config = vars(args)
print(config)



data_path = "/nfs/data%s/%s/data/%s"%(config['server'], config['experimenter'], config['experiment_filename'])
save_path = "/nfs/data%s/%s/data/%s"%(config['server'], config['analyzer'], config['experiment_filename'])

hostname = socket.gethostname()
if hostname[5:] == config['server']:
    data_path = "/data/%s/data/%s"%(config['experimenter'], config['experiment_filename'])
    save_path = "/data/%s/data/%s"%(config['analyzer'], config['experiment_filename'])
    
if not os.path.exists(save_path):
   os.makedirs(save_path)


data = h5py.File(os.path.join(save_path, "NMF_merge.h5"), "r+")
fluoTraces= np.array(data["A_all_merge"])


print("Baseline correction...")
min_w = 2 * 60 * 15


A_baseline =np.zeros_like(fluoTraces)
A_dFF =np.zeros_like(fluoTraces)
for i in range(fluoTraces.shape[0]):
    neural_activity = fluoTraces[i]
    baseline = window_baseline(neural_activity, min_w, 10)
    A_baseline[i] = baseline
    
    neural_activity_gain_corrected = neural_activity
    baseline_gain_corrected = baseline
    baseline_min = np.nanmin(baseline_gain_corrected)
    
    outlier_threshold = 20
    offset = 0
    if  baseline_min < outlier_threshold:
        offset = outlier_threshold - baseline_min
        
    dFF = (neural_activity_gain_corrected - baseline_gain_corrected)/ (baseline_gain_corrected+offset)
    A_dFF[i] = dFF


A_dF = fluoTraces - A_baseline


data.create_dataset('A_baseline', data=A_baseline)
data.create_dataset('A_dF', data=A_dF)
data.create_dataset('A_dFF', data=A_dFF)
data.close()
print("Baseline correction done!")