

from rtlsdr import RtlSdr
import numpy as np
from numpy.fft import fft

class Rtlsdr:
    def __init__(self):
        self.sdr = RtlSdr()
        self.buffer_limit = 4000*2048
        
    def obs_1s(self,cfreq,sample_rate,fftsize,gain):
        self.sdr.center_freq = cfreq
        self.sdr.sample_rate = sample_rate
        self.sdr.gain = gain

        freq_range = np.linspace( cfreq - sample_rate/2,cfreq + sample_rate/2,fftsize )
        num_1s = int( np.round(sample_rate/fftsize) )
        samples_1s = self.sdr.read_samples(num_1s*fftsize).reshape(num_1s,fftsize)
        samples_1sfft = np.abs( np.sum( np.sqrt( fft(samples_1s)**2 ),axis=0 )/num_1s )
        if samples_1sfft[0] > np.sum(samples_1sfft[1:5])/4:
            samples_1sfft[0] = np.sum(samples_1sfft[1:5])/4
        else:
            pass
        samples_1sfft = np.fft.fftshift(samples_1sfft)
        return freq_range,samples_1sfft

    def power_spectrum(self,obstime,cfreq,sample_rate,fftsize,gain):
        self.sdr.sample_rate = sample_rate
        self.sdr.center_freq = cfreq
        self.sdr.gain = gain
        num_1s = int( np.round(sample_rate/fftsize) )
        loop_time = 0
        freql,samples = self.obs_1s(cfreq,sample_rate,fftsize,gain)
        while loop_time < obstime:
            loop_time += 1
            freql,temp_array = self.obs_1s(cfreq,sample_rate,fftsize,gain)
            samples = np.vstack( (samples,temp_array) )

        fft_samples = np.sum(samples,axis=0)/int( obstime )
        return freql, fft_samples


    def disconnect(self):
        self.sdr.close()

