from labdata.schema import *
import numpy as np
import spikeinterface.full as si
from spikeinterface.sortingcomponents.peak_detection import detect_peaks
from spikeinterface.sortingcomponents.peak_localization import localize_peaks
fname = 'peaks.npy'
fname2 = 'peak_locations.npy'

def extract_spikes(subject, session_name, probe_num):
    files2get = (EphysRecording().ProbeFile() & dict(probe_num=probe_num,
                                                 subject_name=subject,
                                                 session_name=session_name))
    path = (File() & files2get).get()[0].parent # get files from s3 if not local
    print(path)
    found_files = list(path.glob(fname)) + list(path.glob(fname2))
    if len(found_files) == 2:
        print('Spike detection already ran for',subject,session_name)
        peaks = np.load(path / fname)
        peak_locations = np.load(path / fname2)
        return peaks, peak_locations
    else:
        print('Running spike detection for',subject,session_name)

        # preprocessing
        rec = si.read_cbin_ibl(path)
        # restrict to channels with activity in this recording
        #rec = rec.channel_slice(rec.channel_ids[:300])
        rec = rec.channel_slice(rec.channel_ids[:-1]) # remove sync
        rec = si.bandpass_filter(rec)
        rec = si.phase_shift(rec)
        rec = si.common_reference(rec)
        noise_levels = si.get_noise_levels(rec, return_scaled=False)

        peaks = detect_peaks(
            rec,
            method="locally_exclusive",
            detect_threshold=6,
            peak_sign="both",
            noise_levels=noise_levels,
            #n_jobs=-1,
            n_jobs=16,
        )

        peak_locations = localize_peaks(
            rec,
            peaks,
            method="monopolar_triangulation",
            #local_radius_um=75,
            #n_jobs=20,
            n_jobs=16,
        )
        np.save(path / fname, peaks, allow_pickle=False)
        np.save(path / fname2, peak_locations, allow_pickle=False)
        return peaks, peak_locations