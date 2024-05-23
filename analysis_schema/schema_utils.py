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
        fs = (EphysRecording.ProbeSetting() & key).fetch1('sampling_rate')

        t_seconds = peaks['sample_index'] / fs
        amps = np.abs(peaks['amplitude'])
        depth_um = peak_locations['y']
        x = peak_locations['x']
        
        good_t = np.logical_and(t_seconds > t_start_sec, t_seconds < t_end_sec)
        shank = np.digitize(x, [175, 400, 650, 1200]) # group the spikes into the closest shank

        for s in np.unique(shank):
            #if s in (0,1):
            #    good_y =  np.logical_and(depth_um < 4790, depth_um > 3000)
            #elif s in (2,3):
            #    good_y = depth_um < 3500
            #good_y = depth_um < 3600
            good_y = np.ones_like(depth_um, dtype=bool)
            
            inds = np.vstack([good_y, good_t, shank == s])
            inds = np.all(inds, axis=0)

            #spks['shank'] = shank[inds]
            #spks['amps'] = amps[inds]
            #spks['depths_um'] = depth_um[inds]
            #spks['times_s'] = t_seconds[inds] + session_offset_sec*o - t_start_sec
            all_shank.append(shank[inds])
            all_amp.append(amps[inds])
            all_depth.append(depth_um[inds])
            all_t.append(t_seconds[inds] + session_offset_sec*o - t_start_sec)

        session_breaks.append(session_offset_sec*o)
    
    all_shank = np.concatenate(all_shank) 
    all_amp = np.concatenate(all_amp)
    all_depth = np.concatenate(all_depth)
    all_t = np.concatenate(all_t)
    #all_spikes = pd.DataFrame(all_spikes).apply(lambda col: col.explode())
    session_breaks = np.array(session_breaks[1:])
    #return all_spikes, session_breaks
    return all_shank, all_amp, all_depth, all_t, session_breaks