from labdata.schema import *
from labdata import chronic_paper as cp
import numpy as np
from tqdm import tqdm

def get_concatenated_spike_data(spike_detection_keys, t_start_sec, t_end_sec):
    session_offset_sec = t_end_sec - t_start_sec
    all_spikes = []
    session_breaks = []
    spks = {}
    all_shank, all_amp, all_depth, all_t = [],[],[],[]
    for o,key in enumerate(tqdm(spike_detection_keys)):
        peaks, peak_locations = cp.DredgeSpikeDetection().get(key)
        fs = np.float32((EphysRecording.ProbeSetting() & key).fetch1('sampling_rate'))

        t_seconds = np.float32(peaks['sample_index']) / fs
        amps = np.abs(peaks['amplitude'])
        depth_um = peak_locations['y']
        x = peak_locations['x']
        
        good_t = np.logical_and(t_seconds > t_start_sec, t_seconds < t_end_sec)
        shank = np.digitize(x, [175, 400, 650, 1200]) # group the spikes by the closest shank

        for s in np.unique(shank):
            inds = np.vstack([good_t, shank == s])
            inds = np.all(inds, axis=0)

            all_shank.append(shank[inds])
            all_amp.append(amps[inds])
            all_depth.append(depth_um[inds])
            all_t.append(t_seconds[inds] + session_offset_sec*o - t_start_sec)

        session_breaks.append(session_offset_sec*o)
    
    all_shank = np.concatenate(all_shank) 
    all_amp = np.concatenate(all_amp)
    all_depth = np.concatenate(all_depth)
    all_t = np.concatenate(all_t)
    session_breaks = np.array(session_breaks[1:])
    return all_shank, all_amp, all_depth, all_t, session_breaks

def decode_distance(chA,chB,device = None):
    import torch
    if device is None:
        if torch.backends.mps.is_available():
            device = 'mps'
        elif torch.cuda.is_available():
            device = 'cuda'
        else:
            device = 'cpu'
    chA = torch.from_numpy(chA.astype('float32')).to(device)
    chB = torch.from_numpy(chB.astype('float32')).to(device)
    distance = torch.zeros(chA.shape)
    from tqdm import tqdm
    for i in tqdm(range(len(chA)-1)):
        if np.mod(len(chB[(chB > chA[i]) & (chB < chA[i+1])]),2):
            distance[i+1] = +1;
        else:
            distance[i+1] = -1;
    
    return np.cumsum(distance.to('cpu').numpy().flatten())
