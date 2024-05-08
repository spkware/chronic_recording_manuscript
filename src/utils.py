from spks import *
import json
import datajoint as dj
from labdatatools import *

with open('/home/joao/.djchurchland.prefs','r') as fd:
    preferences = json.load(fd)

for k in preferences['datajoint'].keys():
    dj.config[k] = preferences['datajoint'][k]

dj.config["database.use_tls"] = False  # Workaround to use python 3.10
    
dj.conn()

spkschema = dj.schema('2023_chronic_holder_manuscript') 

if not 'database.name' in preferences['datajoint'].keys():
    preferences['datajoint']['database.name'] = 'churchlandlab_datajoint'


def sessionfolder_to_datetime(folder,folderformat = '%Y%m%d_%H%M%S'):
    dd = datetime.datetime.strptime(folder,folderformat)
    return dd #dd.strftime('%Y-%m-%d %H:%M:%S')

def list_probe_paths(directory, probe_regex = '\s*.imec(\d).\s*' ,extension = 'ap.meta'):
    '''
    Return files and separates per probe
    '''
    bin_paths = list(Path(directory).expanduser().rglob('*.'+extension))

    bin_paths_no_extension = [p.with_suffix('') for p in bin_paths]
    unique_files = np.unique(bin_paths_no_extension)
    # probe number is attached to imec
    prbs = []
    for p in bin_paths:
        prbs.append(re.search(probe_regex,str(p)).group())
    # put each probe in a list
    all_probe_dirs = []
    for probe in np.sort(np.unique(prbs)):
        probe_dirs = [str(folder) for folder in bin_paths if probe in str(folder)]
        probe_dirs = natsorted(probe_dirs, key=str)
        all_probe_dirs.append(probe_dirs)
    return all_probe_dirs

def fix_int_on_keys(path):
    fd = h5.File(path,mode='a')
    for k in fd.keys():
        try:
            int(k)
        except:
            nk = str(int(float(k)))
            fd[nk] = fd[k]
            del fd[k]
    del fd

def get_units_for_session(session):
    session_path = session['session_path']
    
    sortingpaths = list_sorting_result_paths(Path(labdata_preferences['paths'][0])/session_path)
    rawpaths = list_probe_paths(Path(labdata_preferences['paths'][0])/session_path,
                                extension='ap.meta')
    assert len(rawpaths) == len(sortingpaths), f'raw paths not the same size as sorting {rawpaths} {sortingpaths}'
                            
    recs = []
    recunits = []
    probes = []
    for ifolder,(folder,metafile) in enumerate(zip(sortingpaths,rawpaths)):
        #TODO: ifolder needs to come from the regexp of folder whereever imec is
        from spks.spikeglx_utils import read_spikeglx_meta
        meta = read_spikeglx_meta(metafile[0])
        probe = dict(probe_serial = str(int(meta['imDatPrb_sn'])),
                     probe_type = 'NP1')
        if int(meta['imDatPrb_type']) in [24]:
            probe['probe_type'] = 'NP2alpha'
        recprobe = dict(session,probe_num = ifolder, sorting_path = folder[0],
                        probe_serial = probe['probe_serial'])
        if (Path(folder[0])/'cluster_waveforms.hdf').exists():
            # make sure it is int, not float index
            fix_int_on_keys(Path(folder[0])/'cluster_waveforms.hdf')
        clu = Clusters(folder = folder[0])
        # check if there are any fields missing and if so ask to recompute
        for field in ['cluster_id', 'depth', 'electrode', 'shank', 'position','isi_contamination', 'firing_rate',
                    'presence_ratio', 'amplitude_cutoff', 'trough_time', 'trough_amplitude',
                    'fw3m', 'trough_gradient', 'peak_gradient', 'peak_time',
                    'peak_amplitude', 'spike_duration', 'polarity', 'n_active_channels',
                    'active_channels']:
            if not field in clu.cluster_info.columns:
                clu.compute_statistics(recompute=True)
                break
        units = []
        for i,unit in clu.cluster_info.iterrows():
            units.append(dict(recprobe,
                            unit = unit.cluster_id,
                            shank = unit.shank,
                            electrode = unit.electrode,
                            depth = unit.depth,
                            #position = unit.position,
                            isi_contamination = unit.isi_contamination,
                            firing_rate = unit.firing_rate,
                            presence_ratio = unit.presence_ratio,
                            amplitude_cutoff = unit.amplitude_cutoff,
                            trough_time = unit.trough_time,
                            trough_amplitude = unit.trough_amplitude,
                            trough_gradient = unit.trough_gradient,
                            fw3m = unit.fw3m,
                            peak_amplitude = unit.peak_amplitude,
                            peak_time = unit.peak_time,
                            peak_gradient = unit.peak_gradient,
                            spike_duration = unit.spike_duration,
                            polarity = unit.polarity,
                            n_active_channels = unit.n_active_channels,
                            active_channels = unit.active_channels))
        del clu
        recs.append(recprobe)
        probes.append(probe)
        recunits.append(units)
    return probes,recs,recunits

@spkschema
class Subject(dj.Manual):
    ''' Experimental subject.'''
    definition = """
    subject_name  : varchar(20)          # unique mouse id
    ---
    subject_dob       : date                 # mouse date of birth
    subject_gender    : enum('M', 'F', 'U')  # sex of mouse - Male, Female, or Unknown
    """

@spkschema
class Subject(dj.Manual):
    ''' Experimental subject.'''
    definition = """
    subject_name  : varchar(20)          # unique mouse id
    ---
    subject_dob       : date                 # mouse date of birth
    subject_gender    : enum('M', 'F', 'U')  # sex of mouse - Male, Female, or Unknown
    """

@spkschema
class Probe(dj.Manual):
    ''' Experimental subject.'''
    definition = """
    probe_serial  = 'unknown' : varchar(20)          # unique probe serial
    ---
    probe_type    : enum('NP1','NP2','NP2alpha')
    """

@spkschema
class EphysSession(dj.Manual):
    definition = """
    -> Subject
    session_datetime              : datetime           # experiment date
    ---
    session_path                  : varchar(256)        # path to the
    """

@spkschema
class RecordingProbe(dj.Imported):
    definition = '''
    -> EphysSession
    probe_num               : tinyint
    ---
    -> Probe 
    sorting_path            : varchar(256)
    raw_path     = NULL           : varchar(256)
    '''
    class Units(dj.Part):
        definition = '''
        -> master
        unit                : int
        shank               : int
        ---
        electrode         = NULL   : int
        depth             = NULL   : float
        position          = NULL   : float
        isi_contamination = NULL   : float 
        firing_rate       = NULL   : float
        presence_ratio    = NULL   : float
        amplitude_cutoff  = NULL   : float
        fw3m              = NULL   : float
        trough_time       = NULL   : float
        trough_amplitude  = NULL   : float
        trough_gradient	  = NULL   : float
        peak_time         = NULL   : float
        peak_amplitude    = NULL   : float
        peak_gradient     = NULL   : float
        spike_duration    = NULL   : float
        polarity          = NULL   : int
        n_active_channels = NULL   : int
        active_channels   = NULL   : blob
        '''
    def make(self,key):
        # seskey = dict(subject = key['subject_name'],
        #              session = key['session_datetime'].strftime('%Y%m%d_%H%M%S'))
        nkey = pd.DataFrame((EphysSession() & key).fetch()).iloc[0]
        probes, recprobes,recunits = get_units_for_session(nkey)
        Probe().insert(probes,skip_duplicates=True)
        RecordingProbe().insert(recprobes,skip_duplicates=True,
                                ignore_extra_fields=True)
        for units in recunits:
            RecordingProbe.Units.insert(units,skip_duplicates=True,
                                        ignore_extra_fields=True)
@spkschema
class SessionStats(dj.Computed):
    #TODO: add Kilosort "good" here
    definition = '''
        -> RecordingProbe
        ---
        n_mua_units = 0 : int
        n_noise_units = 0 : int
        n_sua_units = 0 : int
        n_total_units = 0 :int 
        '''
    def make(self,key):
        data = pd.DataFrame((RecordingProbe.Units() & key).fetch())
        sua_idx = ((data.peak_amplitude>25) & 
                  (data.amplitude_cutoff<0.1) & 
                  #(data.n_active_channels<40) & # the threshold doesnt play well with noise.
                  (data.isi_contamination<0.1))
        
        ses = dict(key,n_total_units = len(data),
                            n_mua_units = np.sum(~sua_idx),
                            n_sua_units = np.sum(sua_idx),
                            n_noise = -1) # not used
        SessionStats().insert1(ses,ignore_extra_fields=True)
