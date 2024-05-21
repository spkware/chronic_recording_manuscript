from labdata.utils import *
from labdata.compute.utils import BaseCompute

class SpikeDetectionCompute(BaseCompute):
    container = 'labdata-spks'
    cuda = False
    name = 'detect'
    url = 'http://github.com/spkware'
    def __init__(self,job_id, allow_s3 = None, delete_results = True, **kwargs):
        super(SpikeDetectionCompute,self).__init__(job_id, allow_s3 = allow_s3)
        self.file_filters = ['.ap.']
        # default parameters
        self.parameters = {}
        self._init_job()
        self.delete_results = delete_results
    
           
    def _secondary_parse(self,arguments,parameter_number = None):
        # Nothing to do here. Could add "probe"
        import argparse
        parser = argparse.ArgumentParser(
            description = 'Extract spikes for the chronic paper.',
            usage = 'detect -a <SUBJECT> -s <SESSION>')
        
        args = parser.parse_args(arguments[1:])
        return
    
    def find_datasets(self, subject_name = None, session_name = None):
        '''
        Searches for subjects and sessions in EphysRecording
        '''
        if subject_name is None and session_name is None:
            print("\n\nPlease specify a 'subject_name' and a 'session_name' to perform spike-sorting.\n\n")
        from labdata.schema import EphysRecording
        from .schema import DredgeSpikeDetection

        keys = []
        if not subject_name is None:
            if len(subject_name) > 1:
                raise ValueError(f'Please submit one subject at a time {subject_name}.')
            if not subject_name[0] == '':
                subject_name = subject_name[0]
        if not session_name is None:
            for s in session_name:
                if not s == '':
                    keys.append(dict(subject_name = subject_name,
                                     session_name = s))
        else:
            # find all sessions that can be spike sorted
            sessions = np.unique(((
                (EphysRecording() & f'subject_name = "{subject_name}"') -
                (DredgeSpikeDetection()))).fetch('session_name'))
            for ses in sessions:
                keys.append(dict(subject_name = subject_name,
                                 session_name = ses))
        datasets = []
        for k in keys:
            datasets += (EphysRecording()& k).proj('subject_name','session_name','dataset_name').fetch(as_dict = True)
        return datasets
        
    def _compute(self):
        from .schema import DredgeSpikeDetection
        from labdata.schema import EphysRecording
        datasets = pd.DataFrame((EphysRecording.ProbeFile() & self.dataset_key).fetch())

        for probe_num in np.unique(datasets.probe_num):
            self.set_job_status(job_log = f'Detecting spikes for probe {probe_num}')
            files = datasets[datasets.probe_num.values == probe_num]
            dset = []
            for i,f in files.iterrows():
                if 'ap.cbin' in f.file_path or 'ap.ch' in f.file_path:
                    dset.append(i)
                elif 'ap.meta' in f.file_path: # requires a metadata file (spikeglx)
                    dset.append(i)
            dset = files.loc[dset]
            if not len(dset):
                print(files)
                raise(ValueError(f'Could not find ap.cbin files for probe {probe_num}'))
            localfiles = self.get_files(dset, allowed_extensions = ['.ap.bin'])
            probepath = list(filter(lambda x: str(x).endswith('bin'),localfiles))
            paths = DredgeSpikeDetection().extract_spikes(subject = self.dataset_key['subject_name'],
                                                          session_name = self.dataset_key['session_name'],
                                                          probe_num = probe_num)
            if self.delete_results:
                # delete results_folder
                for p in paths:
                    os.unlink(p)
                # delete local files if they did not exist
                if not self.files_existed:
                    for f in localfiles:
                        os.unlink(f)
