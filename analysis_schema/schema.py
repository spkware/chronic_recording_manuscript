from labdata.schema import *

paperschema = dj.schema('paper_2024_chronic_holder')

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
     -> [nullable] AnalysisFile.proj(spikes_locations = 'file_path',spikes_locations_storage='storage')
    -> [nullable] AnalysisFile.proj(spikes = 'file_path',spikes_storage='storage')
    '''
    def add_dataset(self,key,associated_filepaths):
        nk = (EphysRecording.ProbeSetting & key).proj().fetch1() # if fetch1 crashes when you get 2 then drop the assert 
        dataset = dict(subject_name = key['subject_name'],
                      session_name = key['session_name'],
                      dataset_name = 'dredge_spike_detection')
        filekeys = AnalysisFile().upload_files(associated_filepaths,dataset)
        # need to find which file is spike locations and which is spikes then addd them to nk
        for f in file_keys:
            if 'spikes_locations' in f['file_path']:
                nk['spikes_locations'] = f['file_path']
                nk['spikes_locations_storage'] = f['storage']
            elif 'spikes' in f['file_path']:
                nk['spikes'] = f['file_path']
                nk['spikes_storage'] = f['storage']
            else:
                raise(ValueError(f'Uploaded the wrong files? {filekeys} you handle the delete'))
        DredgeSpikeDetection.insert1(nk)
