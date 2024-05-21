from labdata.schema import *

paperschema = dj.schema('paper_2024_chronic_holder')

__all__ = ['paperschema','IncludedSubjects',
           'SelectedVideoSessions',
           'SortingChannelMAD',
           'DredgeSpikeDetection',
           'ChronicHolderType',
           'ChronicInsertion']
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
    definition = '''
    -> EphysRecording.ProbeSetting
    ---
     -> [nullable] AnalysisFile.proj(peak_locations = 'file_path',peak_locations_storage='storage')
    -> [nullable] AnalysisFile.proj(peaks = 'file_path',peaks_storage='storage')
    '''

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
                n_jobs=-1)
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
    min_depth: int
    max_depth: int
    bin_s: float
    bin_um: float
    win_step_um: int
    win_scale_um: int
    max_disp_um: int
    '''

@paperschema
class DredgeMotionEstimate(dj.Computed):
    definition = '''
    -> DredgeSpikeDetection
    -> DredgeParams
    ---
    shank_num: int
    displacement: longblob
    spatial_bin_centers_um: longblob
    time_bin_centers_s: longblob
    
    '''
    xlims = [(-1000, 175),
             (175, 400),
             (400, 650),
             (650, 1200)]

    def make(self, key):
        from dredge.dredge_ap import register
        n_shanks = (EphysRecording.ProbeSetting * Probe & key).fetch1('probe_n_shanks')
        print(key)
        print(n_shanks)
        for xlim, shank in zip(self.xlims, range(n_shanks)):
            spikes_query = (DredgeSpikeDetection & key)
            peaks_path = (AnalysisFile & spikes_query.proj(file_path='peaks')).get()[0]
            peak_locations_path = (AnalysisFile & spikes_query.proj(file_path='peak_locations')).get()[0]

            peaks = np.load(peaks_path)
            peak_locations = np.load(peak_locations_path)

            srate = (EphysRecording.ProbeSetting() & key).fetch1('sampling_rate')
            t_seconds = peaks['sample_index'] / srate
            amps = np.abs(peaks['amplitude'])
            depth_um = peak_locations['y']
            x = peak_locations['x']
        
            dredge_params = (DredgeParams & key).fetch1() 

            good_y = np.vstack([depth_um >= dredge_params.pop('min_depth'), depth_um <= dredge_params.pop('max_depth')])
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
            motion_est, _ = register(amps, depth_um, t_seconds, **dredge_params)

            key['displacement'] = motion_est.displacement
            key['spatial_bin_centers_um'] = motion_est.spatial_bin_centers_um
            key['time_bin_centers_s'] = motion_est.time_bin_centers_s
            key['shank_num'] = shank

            self.insert1(key)

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
