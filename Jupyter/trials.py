import numpy as np
import h5py
import errno
import os
import salpa
import urllib.request
import io

def geturltrial(tri, pfx = 'ephysdata'):
    #url = f'file:///home/wagenaar/progs/em170428/data/ephysdata-{tri}.h5'
    password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    topurl = 'https://leechem.caltech.edu/170428'
    password_mgr.add_password(None, topurl, 'wagenaar', 'x0Vj6Qb8Aj2J')
    handler = urllib.request.HTTPBasicAuthHandler(password_mgr)
    opener = urllib.request.build_opener(handler)
    url = f'{topurl}/trialdata/{pfx}-{tri}.h5'
    with opener.open(url) as resp:
        L = resp.getheader('content-length')
        if L:
            L = int(L)
            BS = max(4096, L//20)
        else:
            BS = 512*1024
        buf = io.BytesIO()
        N = 0
        while True:
            part = resp.read(BS)
            if not part:
                break
            buf.write(part)
            N += len(part)
            if L:
                progress = f'{(100*N)//L}% of {L} bytes'
            else:
                progress = f'{N} bytes'
            print(f'Read {progress} for trial {tri}', end='\r')
        print(f'Trial {tri} download complete                     ')
        return h5py.File(buf, 'r')

class VSD:
    def __init__(self):
        '''Do not use this constructor. For use by TRIAL only.'''
        self.tt = None
        self.data = {}

    def timestamps(self):
        '''TIMESTAMPS - Timestamps for the VSD data
        times, units = TIMESTAMPS() returns the timestamps for the frames of
        the VSD traces.
        TIMES is returned in seconds. To make that abundantly clear, "s"
        is returned in UNITS.'''
        return (self.tt, 's')

    def keys(self):
        '''KEYS - Names of ROIs in the VSD data
        kk = KEYS() returns a list of the names of the ROIs in the dataset.
        These can be used to retrieve individual traces with TRACE.'''
        return self.data.keys()

    def trace(self, k, tau=100, poly=None):
        '''TRACE - Retrieve data from a single ROI
        dff, units = TRACE(key) retrieves the data from the named ROI.
        By default, the data are filtered through SALPA with a time constant
        of τ = 100 samples.
        dff, units = TRACE(key, tau=None) returns raw data.
        dff, units = TRACE(key, poly=N) subtracts a polynomial of degree N
        from the raw trace instead.
        DFF is the relative fluorescence change (dF/F) as a percentage. To
        make that abundantly clear, "%" is returned in UNITS.'''
        y = self.data[k]
        if poly is not None:
            x = np.arange(len(y))
            p = np.polyfit(x, y, poly)
            for k in range(poly+1):
                y -= p[k] * x**(poly-k)
        elif tau is not None:
            s = salpa.Salpa(tau=tau)
            y = s.apply(y)
        return (y, '%')

class EPhys:
    def __init__(self):
        '''Do not use this constructor. For use by TRIAL only.'''
        self.tt = None
        self.data = None
        self.units = 'mV'
        self.roi = None

    def timestamps(self):
        '''TIMESTAMPS - Timestamps for the electrophysiology data
        times, units = TIMESTAMPS() returns the timestamps for the frames of
        the electrophysiology traces.
        TIMES is returned in seconds. To make that abundantly clear, "s"
        is returned in UNITS.'''
        return (self.tt, 's')

    def trace(self):
        '''TRACE - Retrieve data from electrophysiology
        yy, units = TRACE() retrieves the data from the intracellular electrode.
        from the raw trace instead.
        For stimuli, YY is a current measured in nanoamperes; for recordings,
        a voltage measured in millivolts. To make that abundantly clear, "nA"
        or "mV" is returned in UNITS.'''
        y = self.data
        return (y, self.units)

    def correspondingROI(self):
        '''CORRESPONDINGROI - ROI corresponding to this recording.
        roi = CORRESPONDINGROI() returns the ROI name corresponding to this
        electrophysiology trace, if there is one. That ROI can then be
        looked up in the VSD data.'''
        return self.roi

usedchannels_ = {
    # Based on Yusuke's spreadsheet "170428_ephys_summary.xlsx"
    6: { 2: ('CV', 'p') },
    8: { 2: ('CV', 'p') },
    9: { 1: ('P_VL', 'im') },
    10: { 1: ('P_VL', 'im') },
    11: { 1: ('P_VL', 'im') },
    12: { 1: ('P_VR', None) },
    15: { 1: ('AE_R', 'bx') },
    17: { 1: ('AE_R', 'bx') },
    }

channelnames_ = {
    1: ('ME1_10Vm', 'ME1_I'),
    2: ('ME2_Vm', 'ME2_I')
    }

def makestr_(a):
    def unpack(x):
        if type(x)==np.ndarray:
            return x[0]
        else:
            return x
    return ''.join(chr(unpack(x)) for x in a['value'])

def lookupchannel_(f, name):
    kk = f['ephys_ch']['value'].keys()
    for k in kk:
        if name==makestr_(f['ephys_ch']['value'][k]):
            return int(k[1:])

class Trial:
    def __init__(self, trial):
        '''TRIAL - Load VSD and electrophysiology data
        x = TRIAL(n) loads data from the given trial. Usable numbers are:
          6 - Swim trial
          8 - Swim and crawl trial
          9 - Local bend (P_VL) trial
         10 - Local bend (P_VL) trial
         11 - Local bend (P_VL) trial
         12 - Local bend (P_VR) trial
         15 - Crawl trial
         17 - Crawl trial
        The returned object has methods VSD, STIMULI, and INTRACELLULAR
        to retrieve the voltage-dye data, stimuli, and intracellular
        recordings associated with the trial respectively.'''
        self.trial = trial
        #here = os.path.dirname(__file__)
        #ifn = f'{here}/../data/ephysdata-{trial}.h5'
        self.vsddata = VSD()
        self.ephysrec = {}
        self.ephysstim = {}
        with geturltrial(trial) as f:
            ids = f['vsd_id']['value']
            dat = np.array(f['vsd_F']['value'])
            for idx in ids.keys():
                if idx=='dims':
                    continue
                vsdid = makestr_(ids[idx])
                vsddata = dat[int(idx[1:]),:]
                vsddata = 100 * (vsddata / np.mean(vsddata) - 1)
                self.vsddata.data[vsdid] = vsddata
            tt = np.array(f['vsd_t']['value']['_1']['value'])
            self.vsddata.tt = tt.flatten()
            if trial in usedchannels_:
                ch = usedchannels_[trial]
                for elc,use in ch.items():
                    rec = EPhys()
                    krec = channelnames_[elc][0]
                    idx = lookupchannel_(f, krec)
                    rec.tt = np.array(f['ephys_t']['value']).flatten()
                    rec.data = np.array(f['ephys_v']['value'])[idx,:]
                    rec.units = makestr_(f['ephys_u']['value'][f'_{idx}'])
                    rec.roi = use[1]
                    self.ephysrec[use[0]] = rec

                    stm = EPhys()
                    kstim = channelnames_[elc][1]
                    idx = lookupchannel_(f, kstim)
                    stm.tt = np.array(f['ephys_t']['value']).flatten()
                    stm.data = np.array(f['ephys_v']['value'])[idx,:]
                    stm.units = makestr_(f['ephys_u']['value'][f'_{idx}'])
                    stm.roi = use[1]
                    self.ephysstim[use[0]] = stm

    def vsd(self):
        '''VSD - Access to voltage-sensitive dye data
        x = VSD() returns access to all the VSD traces in the trial.
        The returned object has methods TIMESTAMPS, KEYS, and TRACE
        to retrieve the contained data.'''
        return self.vsddata

    def stimuli(self):
        '''STIMULI - Access to stimulus data
        s = STIMULI() returns a dict mapping canonical neuron names
        to stimuli delivered to those neurons. The values in the dict
        are objects with methods TIMESTAMPS, TRACE, and CORRESPONDINGROI
        that can be used to retrieve the contained data.'''
        return self.ephysstim

    def intracellular(self):
        '''INTRACELLULAR - Access to intracellular recordings
        s = INTRACELLULAR() returns a dict mapping canonical neuron names
        to recordings made intracellularly from those neurons. The values
        in the dict are objects with methods TIMESTAMPS, TRACE, and
        CORRESPONDINGROI that can be used to retrieve the contained data.'''
        return self.ephysrec

class Coherence:
    def __init__(self, trial):
        '''COHERENCE - Load coherence data for trial.
        x = COHERENCE(n) loads coherence data for the given trial.
        Usable numbers are:
          6 - Swim trial
          8 - Swim and crawl trial
          9 - Local bend (P_VL) trial
         10 - Local bend (P_VL) trial
         11 - Local bend (P_VL) trial
         12 - Local bend (P_VR) trial
         15 - Crawl trial
         17 - Crawl trial
        The result is a class with (read-only) member variables:
          DATA - A dict containing:
            f - frequency at which results were computed
            psd - power spectral densities of all signals at that frequency
            coh - complex coherence
            mag - coherence magnitudes (ditto)
            phase - coherence phases (ditto)
            mag_lo and mag_hi - low and high end of confidence intervals 
                                for mag
            phase_lo and phase_hi - ditto for phase
            thr - threshold for significance
            cc - color map (signals with mag<thr will be gray)
          ROIIDS - A list with the names of the ROIs
          EXTRA - A dict containing:
            f - frequency vector
            psd - power spectral densities of all signals at all freqs
            refpsd - psds of reference at all frequencies
            tt - time vector
            sig - detrended debleached signals at those times
            ref - detrended reference at those times
            img - raw data for first frame within the time window, 
                  ventral camera
            xx and yy - transformed coordinates for that image
            imgb - raw data for first frame within the time window, 
                   dorsal camera
            xxb and yyb - transformed coordinates for that image
            rois - original rois'''

        self.trial = trial
        self.data = {}
        self.extra = {}
        with geturltrial(trial, 'coh') as coh:
            coh = coh['coh']['value']
            for fld in ['psd', 'mag', 'phase',
                        'mag_lo', 'mag_hi', 'phase_lo', 'phase_hi']:
                self.data[fld] = np.array(coh[fld]['value']).flatten()
            self.data['cc'] = np.array(coh['cc']['value'])
            self.data['coh'] = self.data['mag'] * np.exp(1j*self.data['phase'])
            for fld in ['f', 'thr', 'pthr']:
                self.data[fld] = np.array(coh[fld]['value']).flatten()[0]
            for fld in ['img', 'xx', 'yy', 'tt', 'sig', 'ref',
                        'f', 'refpsd', 'psd']:
                self.extra[fld] = np.array(coh['extra']['value'][fld]['value'])
            self.extra['rois'] = {}
            rois = coh['extra']['value']['rois']['value']
            for k in rois.keys():
                if k.startswith('_'):
                    r = int(k[1:])
                    self.extra['rois'][r] = np.array(rois[k]['value'])
            N = len(self.data['coh'])
            self.roiids = []
            for n in range(N):
                if n<26:
                    self.roiids.append(chr(ord('a') + n))
                else:
                    self.roiids.append(chr(ord('a') + (n//26-1))
                                       + chr(ord('a') + (n%26)))
