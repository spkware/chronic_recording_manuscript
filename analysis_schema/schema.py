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
                [3, 'NP24 freely moving'],]

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