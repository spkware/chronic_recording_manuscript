from labdata.schema import *
from .schema_utils import * 

paperschema = dj.schema('paper_2024_chronic_holder')

__all__ = ['paperschema','IncludedSubjects',
           'SelectedVideoSessions',
           'SortingChannelMAD',
           'DredgeSpikeDetection',
           'ChronicHolderType',
           'ChronicInsertion',
           'DredgeMotionEstimate',
           'DredgeParams',
           'ConcatenatedSpikes',
           'LocomotionBehaviorTreadmill']
@paperschema
class IncludedSubjects(dj.Manual):
    # Lists the mice included in the study
    definition = '''
    -> Subject
    ---
    mouse_id = NULL : int 
    descriptor = NULL : varchar(52)
    '''

@paperschema
class SelectedVideoSessions(dj.Manual):
    # Table to hold which video sessions were used
    definition = '''
    -> Session
    ---
    recording_type: enum('top', 'side')   # side will be used to map movement
    '''

    def populate_sync(self):
        for session_key in self:
            insert_video(session_key)

@paperschema
class SortingChannelMAD(dj.Computed):
    definition = """
    -> SpikeSorting
    average_mad  : double
    std_mad      : double
    ---
    mad          : longblob
    """

    def make(self,key):
        seg = (SpikeSorting.Segment()*EphysRecording.ProbeSetting()*ProbeConfiguration() & key).fetch1()
        from spks import mad
        m = mad(seg['segment']*seg['probe_gain']).astype(np.float32)
        self.insert1(dict(key,
                          mad =m,
                         average_mad = np.mean(m),
                         std_mad = np.std(m)))

@paperschema
class DredgeSpikeDetection(dj.Manual):
    # Table to hold which video sessions were used
    # FIXME: This table accidentally populates a second shank for NP1's sometimes. Use channelmaps and pick the closest shank in the future.
    definition = '''
    -> EphysRecording.ProbeSetting
    ---
     -> [nullable] AnalysisFile.proj(peak_locations = 'file_path',peak_locations_storage='storage')
    -> [nullable] AnalysisFile.proj(peaks = 'file_path',peaks_storage='storage')
    '''
    def get(self, key):
        paths = (self & key).fetch1('peaks','peak_locations')
        peak_locations_path, peaks_path = (AnalysisFile() & 'file_path IN ' + str(paths)).get()

        peaks = np.load(peaks_path)
        peak_locations = np.load(peak_locations_path)
        return peaks, peak_locations
    
    def get_spikes(self,shank_num = None):
        keys = []
        xlims = [(-1000, 175),
                 (175, 400),
                 (400, 650),
                 (650, 1200)]
        for key in self:
            print(key)
            paths = (DredgeSpikeDetection & key).fetch1('peaks','peak_locations')
            peak_locations_path, peaks_path = (AnalysisFile() & 'file_path IN ' + str(paths)).get()

            peaks = np.load(peaks_path)
            peak_locations = np.load(peak_locations_path)
            srate = np.float32((EphysRecording.ProbeSetting() & key).fetch1('sampling_rate'))
            t_seconds = peaks['sample_index'].astype(np.float32) / srate
            amps = np.abs(peaks['amplitude'])
            depth_um = peak_locations['y']
            x = peak_locations['x']
            if shank_num is None:
                good_x = np.arange(len(t_seconds))
            else:
                good_x = np.logical_and(x >= xlims[shank_num][0], x <= xlims[shank_num][1])
            key['t_seconds'] = t_seconds[good_x]
            key['amps'] = amps[good_x]
            key['depth_um'] = depth_um[good_x]
            keys.append(key)
        return pd.DataFrame(keys)
    
    def plot_raster(self, subject, session_name, probe_num, shank_num,**kwargs):
        from spks.viz import plot_drift_raster

        xlims = [(-1000, 175),
                 (175, 400),
                 (400, 650),
                 (650, 1200)]

        key = dict(subject_name=subject,
                       session_name=session_name,
                       probe_num=probe_num)
        spikes_query = (DredgeSpikeDetection & key)
        peaks_path = (AnalysisFile & spikes_query.proj(file_path='peaks')).get()[0]
        peak_locations_path = (AnalysisFile & spikes_query.proj(file_path='peak_locations')).get()[0]

        peaks = np.load(peaks_path)
        peak_locations = np.load(peak_locations_path)

        srate = np.float32((EphysRecording.ProbeSetting() & key).fetch1('sampling_rate'))
        t_seconds = peaks['sample_index'].astype(np.float32) / srate
        amps = np.abs(peaks['amplitude'])
        depth_um = peak_locations['y']
        x = peak_locations['x']

        good_x = np.logical_and(x >= xlims[shank_num][0], x <= xlims[shank_num][1])

        t_seconds = t_seconds[good_x]
        amps = amps[good_x]
        depth_um = depth_um[good_x]
            
        plot_drift_raster(t_seconds, depth_um, amps,**kwargs)
        

    def extract_spikes(self, subject, session_name, probe_num):
        """
        This function will download the binary file for a probe, perform spike extraction/localization,
        save those to .npy files, and add these files to the AnalysiFile and DredgeSpikeDetection tables.
        """
        import numpy as np
        import spikeinterface.full as si
        from spikeinterface.sortingcomponents.peak_detection import detect_peaks
        from spikeinterface.sortingcomponents.peak_localization import localize_peaks
        from spikeinterface.core import set_global_job_kwargs
        set_global_job_kwargs(mp_context = "fork")
        fname = 'peaks.npy'
        fname2 = 'peak_locations.npy'

        key = dict(subject_name=subject,
                       session_name=session_name,
                       probe_num=probe_num)
        has_file_on_database = len(self & key) == 1
        if has_file_on_database:
            print('Spike detection already ran for',subject,session_name)
            return
            #peaks = np.load(path / fname)
            #peak_locations = np.load(path / fname2)
            #return peaks, peak_locations
        else:
            print('Running spike detection for',subject,session_name)
            files2get = (EphysRecording().ProbeFile() & dict(probe_num=probe_num,
                                                         subject_name=subject,
                                                         session_name=session_name))
            path = (File() & files2get).get()[0].parent # get files from s3 if not local
            print(path)

            # preprocessing
            rec = si.read_cbin_ibl(path)
            # restrict to channels with activity in this recording
            #rec = rec.channel_slice(rec.channel_ids[:300])
            rec = rec.channel_slice(rec.channel_ids[:-1]) # remove sync
            rec = si.bandpass_filter(rec)
            rec = si.phase_shift(rec)
            rec = si.common_reference(rec)
            noise_levels = si.get_noise_levels(rec, return_scaled=False)
            if (path/"preprocessed").exists():
                import shutil
                shutil.rmtree(path/"preprocessed")
            rec_cache = rec.save(folder=path/"preprocessed", progress_bar=True, 
                                 n_jobs=-1)
            peaks = detect_peaks(
                rec_cache,
                method="locally_exclusive",
                detect_threshold=6,
                peak_sign="both",
                noise_levels=noise_levels,
                n_jobs=-1)
            
            peak_locations = localize_peaks(
                rec_cache,
                peaks,
                method="monopolar_triangulation",
                #local_radius_um=75,
                n_jobs=8)
            import shutil
            shutil.rmtree(path/"preprocessed")
            np.save(path / fname, peaks, allow_pickle=False)
            np.save(path / fname2, peak_locations, allow_pickle=False)

            # add the files to the database
            key = (EphysRecording.ProbeSetting & key).proj().fetch1()
            self.add_dataset(key, [path / fname, path / fname2])
            return [path / fname, path / fname2]
        
    def add_dataset(self,key,associated_filepaths):
        nk = (EphysRecording.ProbeSetting & key).proj().fetch1() # if fetch1 crashes when you get 2 then drop the assert 
        dataset = dict(subject_name = key['subject_name'],
                      session_name = key['session_name'],
                      dataset_name = f'dredge_spike_detection/probe_{key["probe_num"]}')
        filekeys = AnalysisFile().upload_files(associated_filepaths,dataset, force=False)
        # need to find which file is spike locations and which is spikes then addd them to nk
        for f in filekeys:
            if 'peak_locations' in f['file_path']:
                nk['peak_locations'] = f['file_path']
                nk['peak_locations_storage'] = f['storage']
            elif 'peaks' in f['file_path']:
                nk['peaks'] = f['file_path']
                nk['peaks_storage'] = f['storage']
            else:
                raise(ValueError(f'Uploaded the wrong files? {filekeys} you handle the delete'))
        DredgeSpikeDetection.insert1(nk)

@paperschema
class DredgeParams(dj.Manual):
    definition = '''
    params_id: int
    ---
    bin_s: float
    bin_um: float
    win_step_um: int
    win_scale_um: int
    max_disp_um: int
    '''

@paperschema
class DredgeMotionEstimate(dj.Manual):
    definition = '''
    -> DredgeSpikeDetection
    -> DredgeParams
    shank_num: tinyint
    ---
    displacement: longblob
    spatial_bin_centers_um = NULL : longblob
    time_bin_centers_s: longblob
    min_spike_depth_um: float           # spikes below this depth were excluded before running DREDGE
    max_spike_depth_um: float           # spikes above this depth were excluded before running DREDGE
 
    '''
    xlims = [(-1000, 175),
             (175, 400),
             (400, 650),
             (650, 1200)]

    def insert_one_session(self, key, min_spike_depth, max_spike_depth):
        from dredge.dredge_ap import register
        #key = key.copy()
        n_shanks = (EphysRecording.ProbeSetting * Probe & key).fetch1('probe_n_shanks')
        for xlim, shank, ymin, ymax in zip(self.xlims, range(n_shanks), min_spike_depth, max_spike_depth):
            spikes_query = (DredgeSpikeDetection & key)
            #peaks_path = (AnalysisFile & spikes_query.proj(file_path='peaks')).get()[0]
            #peak_locations_path = (AnalysisFile & spikes_query.proj(file_path='peak_locations')).get()[0]
            paths = spikes_query.fetch1('peaks','peak_locations') #+ AnalysisFile & spikes_query.proj(file_path='peak_locations')
            peak_locations_path, peaks_path = (AnalysisFile() & 'file_path IN ' + str(paths)).get()

            peaks = np.load(peaks_path)
            peak_locations = np.load(peak_locations_path)

            srate = (EphysRecording.ProbeSetting() & key).fetch1('sampling_rate')
            t_seconds = peaks['sample_index'] / srate
            amps = np.abs(peaks['amplitude'])
            depth_um = peak_locations['y']
            x = peak_locations['x']
        
            dredge_params = (DredgeParams & key).fetch1() 

            #SpikeDepthLims() & key
            good_y = np.vstack([depth_um >= ymin, depth_um <= ymax])
            good_x = np.vstack([x >= xlim[0], x <= xlim[1]])
            inds = np.vstack([good_x, good_y])
            inds = np.all(inds, axis=0)
            if np.sum(inds) == 0:
                continue

            t_seconds = t_seconds[inds]
            amps = amps[inds]
            depth_um = depth_um[inds]

            dredge_params.pop('params_id')
            dredge_params.pop('max_disp_um')
            motion_est, _ = register(amps, depth_um, t_seconds, pbar=False, **dredge_params)

            key['displacement'] = motion_est.displacement
            key['spatial_bin_centers_um'] = motion_est.spatial_bin_centers_um
            key['time_bin_centers_s'] = motion_est.time_bin_centers_s
            key['min_spike_depth_um'] = ymin
            key['max_spike_depth_um'] = ymax
            key['shank_num'] = shank

            self.insert1(key, skip_duplicates=True)
    
    def populate_subject(self,subject, probe_num, dredge_params_id, min_spike_depth=None, max_spike_depth=None, n_workers=1):
        from tqdm import tqdm
        from multiprocessing.pool import Pool
        from functools import partial

        keys = (DredgeSpikeDetection * DredgeParams & dict(probe_num=probe_num, subject_name=subject, params_id=dredge_params_id)).fetch('KEY')
        keys_to_insert = []
        for k in keys:
            n_shanks = (EphysRecording.ProbeSetting * Probe & k).fetch1('probe_n_shanks')
            if len(self & k) < n_shanks: # if all shanks are already in the table, skip
                keys_to_insert.append(k)
        #for key in tqdm(keys_to_insert):
        #    self.insert_one_session(key, min_spike_depth, max_spike_depth)
        #TODO: make this work so that each worker just populates one shank, rather than looping inside the worker
        worker = partial(self.insert_one_session, min_spike_depth=min_spike_depth, max_spike_depth=max_spike_depth)
        with Pool(n_workers) as p:
            _ = list(tqdm(p.imap(worker, keys_to_insert), total=len(keys_to_insert)))

@paperschema
class ConcatenatedSpikes(dj.Manual):
    definition = '''
    -> Subject
    -> ProbeConfiguration
    probe_num : tinyint
    shank_num : tinyint
    ---
    spike_times_s : longblob
    spike_depths_um : longblob
    spike_amps : longblob
    t_start : float
    t_end : float
    session_breaks : longblob
    '''
    class IncludedSessions(dj.Part):
    # FIXME: make this part table depend on DredgeSpikeDetection instead of EphysRecording, so that if a spike detection is deleted, the concatenated spikes entry is also deleted
        definition = '''
        -> master
        concatenation_order : smallint
        ---
        -> EphysRecording
        '''
    class DredgeResults(dj.Part):
        definition = '''
        -> master
        ---
        displacement: longblob
        spatial_bin_centers_um = NULL : longblob
        time_bin_centers_s: longblob
        spatial_bin_edges_um = NULL : longblob
        time_bin_edges_s : longblob
        min_spike_depth_um: float           # spikes below this depth were excluded before running DREDGE
        max_spike_depth_um: float           # spikes above this depth were excluded before running DREDGE
        '''
    
    def plot_raster(self, key, corrected=False, overlay_dredge=False, ax=None,**kwargs):
        from spks.viz import plot_drift_raster
        if ax is None:
            import matplotlib.pyplot as plt
            ax = plt.gca()
        spks = (self & key).fetch1()
        if corrected or overlay_dredge:
            dredge_rez = (self.DredgeResults & key).fetch1()
        if not corrected:
            plot_drift_raster(spks['spike_times_s'],
                              spks['spike_depths_um'],
                              spks['spike_amps'],**kwargs)
            plt.ylabel('Depth (um)')
        else:
            import dredge.motion_util as mu
            if dredge_rez['spatial_bin_centers_um'] is None:
                motion_estimate = mu.RigidMotionEstimate(displacement=dredge_rez['displacement'],
                                                         time_bin_edges_s=dredge_rez['time_bin_edges_s'],
                                                         time_bin_centers_s=dredge_rez['time_bin_centers_s'])
            else:
                motion_estimate = mu.NonrigidMotionEstimate(displacement=dredge_rez['displacement'],
                                                            time_bin_edges_s=dredge_rez['time_bin_edges_s'],
                                                            spatial_bin_edges_um=dredge_rez['spatial_bin_edges_um'],
                                                            time_bin_centers_s=dredge_rez['time_bin_centers_s'],
                                                            spatial_bin_centers_um=dredge_rez['spatial_bin_centers_um'])

            depths_corrected = motion_estimate.correct_s(spks['spike_times_s'].astype(np.float32), spks['spike_depths_um'].astype(np.float32), grid=False)
            plot_drift_raster(spks['spike_times_s'],
                              depths_corrected,
                              spks['spike_amps'],**kwargs)
            plt.ylabel('Corrected depth (um)')
        if overlay_dredge:
            t = dredge_rez['time_bin_centers_s']
            d = dredge_rez['displacement']
            bins = dredge_rez['spatial_bin_centers_um']
            if bins is None:
                middle = np.mean(ax.get_ylim())
                plt.plot(t, d + middle) # plot the displacement for rigid estimates
            else:
                plt.plot(t, bins + d.T) # plot the displacement for nonrigid estimates

        plt.vlines(spks['session_breaks'], *ax.get_ylim(), linewidth=.5, linestyles='--', colors='black', label='Session breaks')

    def create_entries(self, subject, session_names, probe_num, configuration_id=None, t_start=300, t_end=360, dredge_min_depths=None, dredge_max_depths=None, **dredge_params):
        from .schema_utils import get_concatenated_spike_data
        from dredge.dredge_ap import register

        session_names = [dict(session_name=s) for s in session_names]
        spike_detection_keys = (Session() * DredgeSpikeDetection() * EphysRecording.ProbeSetting() & dict(configuration_id=configuration_id,
                                                                                                             probe_num=probe_num,
                                                                                                             subject_name=subject) & session_names).proj('probe_id').fetch(order_by='session_datetime', as_dict=True)
                                                                                 

        insertiondict = dict(subject_name=subject,
                             probe_id=spike_detection_keys[0]['probe_id'],
                             probe_num=probe_num,
                             configuration_id=configuration_id,
                             t_start=t_start,
                             t_end=t_end,)
        if len(self & insertiondict) > 0:
            print(f'Already inserted concatenated spikes for subject {subject}, probe {probe_num}, configuration {configuration_id}')
            return

        shank, amp, depth, t, session_breaks = get_concatenated_spike_data(spike_detection_keys, t_start, t_end)

        insertiondict['session_breaks'] = session_breaks
        
        nshanks = len(np.unique(shank))
        dredge_min_depths = np.array(dredge_min_depths)
        dredge_max_depths = np.array(dredge_max_depths)
        print(f'Creating entries for {nshanks} shanks')
        assert len(dredge_min_depths) == nshanks == len(dredge_max_depths), 'Need to provide min and max spike depths for each shank'

        dredge_keys = dict()
        for i,s in enumerate(np.unique(shank)):
            # 1. Insert spikes for one shank
            insertiondict['shank_num'] = s
            insertiondict['spike_times_s'] = t[shank==s] # grab all spikes on shank
            insertiondict['spike_depths_um'] = depth[shank==s]
            insertiondict['spike_amps'] = amp[shank==s]
            self.insert1(insertiondict)

            # 2. Run Dredge for that shank and insert it
            # grab all spikes on shank but also restrict by depths
            stacked = np.vstack([shank==s, depth >= dredge_min_depths[s], depth <= dredge_max_depths[s]])
            inds = np.all(stacked, axis=0)
            motion_est, _ = register(amp[inds], depth[inds], t[inds], pbar=False, **dredge_params)

            dredge_keys['displacement'] = motion_est.displacement
            dredge_keys['spatial_bin_centers_um'] = motion_est.spatial_bin_centers_um
            dredge_keys['time_bin_centers_s'] = motion_est.time_bin_centers_s
            dredge_keys['spatial_bin_edges_um'] = motion_est.spatial_bin_edges_um
            dredge_keys['time_bin_edges_s'] = motion_est.time_bin_edges_s
            dredge_keys['min_spike_depth_um'] = dredge_min_depths[s]
            dredge_keys['max_spike_depth_um'] = dredge_max_depths[s]
            self.DredgeResults.insert1({**insertiondict, **dredge_keys}, ignore_extra_fields=True)

            # 3. Insert the included recordings
            self.IncludedSessions.insert([{**insertiondict,
                                           'session_name': k['session_name'],
                                           'dataset_name': k['dataset_name'],
                                           'concatenation_order' : i} for i,k in enumerate(spike_detection_keys)], ignore_extra_fields=True)

@paperschema
class ChronicHolderType(dj.Lookup):
    # Table to hold the chronic holder type
    definition = '''
    holder_id: int
    ---
    description : varchar(32)
    '''
    contents = [[0, 'NP1 head fixed'],
                [1, 'NP1 freely moving'],
                [2, 'NP24 head fixed'],
                [3, 'NP24 freely moving'],
                [4, 'NP24alpha head fixed'],
                [5, 'NP24alpha freely moving'],]

@paperschema
class ChronicInsertion(dj.Manual):
    definition = '''
    -> ProbeInsertion
    ---
    -> ChronicHolderType
    '''
    class TargetedRegion(dj.Part):
        definition = '''
        -> master
        -> Atlas.Region
        ---
        hemisphere: enum('left', 'right')
        '''

@paperschema
class LocomotionBehaviorTreadmill(dj.Manual):
    definition = '''
    -> EphysRecording
    ---
    ch_sma         : tinyint
    ch_encoder0    : tinyint
    ch_encoder1    : tinyint
    cm_conversion  : double
    time     = NULL      : longblob
    velocity = NULL      : longblob
    position = NULL      : longblob
    '''
    def insert_session(self,key, ch_sma = 7, ch_encoder0 = 4, ch_encoder1=5):
        # select data
        key = (Dataset() & (EphysRecording() & key)).proj().fetch1()
        ephys = (DatasetEvents.Digital() & key & 'stream_name = "imec0"').fetch(as_dict = True)
        srate = float((EphysRecording.ProbeSetting() & key & 'probe_num = 0').fetch1('sampling_rate'))
        if not len(ephys):
            print('Could not find sync for probe 0 - check session.')
            return None
        nidqfiles = (File() & (Dataset.DataFiles() & key & 'file_path LIKE "%.nidq.%"')).get()
        for f in nidqfiles:
            if str(f).endswith('.bin'):
                from spks.spikeglx_utils import load_spikeglx_binary
                dat,meta = load_spikeglx_binary(f)
                from spks.sync import unpackbits_gpu
                onsets,offsets = unpackbits_gpu(dat[:,-1])
        from spks.sync import interp1d
        ni_ap_interp = interp1d(onsets[ch_sma],ephys[0]['event_values'],fill_value='extrapolate')
        # extract the encoder traces this takes a while, could probably be made faster
        chA = ni_ap_interp(onsets[ch_encoder0])
        chB = ni_ap_interp(onsets[ch_encoder1])    
        position = decode_distance(chA,chB)
        enc_time = np.arange(0,len(dat)/srate,1/200.) # 200Hz resolution
        distance = interp1d(chA/srate, position, fill_value='extrapolate')(enc_time)
        # smooth with a gaussian
        from scipy.ndimage import gaussian_filter
        factor = 600./40 # roughly the encoder_pulses per cm 
        distance = gaussian_filter(distance,20)/factor
        velocity = np.diff(distance)/np.diff(enc_time[:2])
        self.insert1(dict(key,time = enc_time,
                          velocity = velocity,
                          position = distance,
                          cm_conversion = factor,
                          ch_encoder0 = ch_encoder0, 
                          ch_encoder1 = ch_encoder1,
                          ch_sma = ch_sma),ignore_extra_fields = True)
        return 
