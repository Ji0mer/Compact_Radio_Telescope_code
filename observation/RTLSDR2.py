#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Not titled yet
# GNU Radio version: 3.10.1.1

from gnuradio import blocks
from gnuradio import fft
from gnuradio.fft import window
from gnuradio import gr
from gnuradio.filter import firdes
import sys
import signal
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
import osmosdr
import time

para = ArgumentParser(description="Set parameters for this observation.")

para.add_argument("--freq","-f",help="center frequency of this observation.",default=1420e6)

para.add_argument("--bandwidth","-b",help="sample rate/ bandwidth of this observation.",default=1.8e6)

para.add_argument("--gain","-g",help="gain of observation.",default=50)

para.add_argument("--obstime","-o",help="observation time of this observation.",default=300)

para.add_argument("-fftsize","-fs",help="fftsize of this observation.",default=4096)

para.add_argument("--filename1","-n1",help="saved file name.",default="SDR1-data.txt")

para.add_argument("--filename2","-n2",help="saved file name.",default="SDR2-data.txt")

para.add_argument("--path","-p",help="file path",default="./")

paras = para.parse_args()

file_name1 = paras.path + time.strftime("%Y-%m-%d-%H:%M:%S_", time.localtime()) + str(paras.filename1)

file_name2 = paras.path + time.strftime("%Y-%m-%d-%H:%M:%S_", time.localtime()) + str(paras.filename2)

print(file_name1,file_name2)

class RTLSDR(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Not titled yet", catch_exceptions=True)

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = float( paras.bandwidth )
        self.fftsize = fftsize = int( paras.fftsize )
        self.second = second = int( samp_rate/fftsize )
        self.obs_time = obs_time = int( paras.obstime )
        self.gain = gain = float( paras.gain )
        self.cfreq = cfreq = float( paras.freq )

        ##################################################
        # Blocks
        ##################################################
        self.rtlsdr_source_0_0 = osmosdr.source(
            args="numchan=" + str(1) + " " + "rtl=1,bias=1"
        )
        self.rtlsdr_source_0_0.set_time_unknown_pps(osmosdr.time_spec_t())
        self.rtlsdr_source_0_0.set_sample_rate(samp_rate)
        self.rtlsdr_source_0_0.set_center_freq(cfreq, 0)
        self.rtlsdr_source_0_0.set_freq_corr(0, 0)
        self.rtlsdr_source_0_0.set_dc_offset_mode(0, 0)
        self.rtlsdr_source_0_0.set_iq_balance_mode(0, 0)
        self.rtlsdr_source_0_0.set_gain_mode(False, 0)
        self.rtlsdr_source_0_0.set_gain(gain, 0)
        self.rtlsdr_source_0_0.set_if_gain(gain, 0)
        self.rtlsdr_source_0_0.set_bb_gain(gain, 0)
        self.rtlsdr_source_0_0.set_antenna('', 0)
        self.rtlsdr_source_0_0.set_bandwidth(0, 0)
        self.rtlsdr_source_0 = osmosdr.source(
            args="numchan=" + str(1) + " " + "rtl=0,bias=1"
        )
        self.rtlsdr_source_0.set_time_unknown_pps(osmosdr.time_spec_t())
        self.rtlsdr_source_0.set_sample_rate(samp_rate)
        self.rtlsdr_source_0.set_center_freq(cfreq, 0)
        self.rtlsdr_source_0.set_freq_corr(0, 0)
        self.rtlsdr_source_0.set_dc_offset_mode(0, 0)
        self.rtlsdr_source_0.set_iq_balance_mode(0, 0)
        self.rtlsdr_source_0.set_gain_mode(False, 0)
        self.rtlsdr_source_0.set_gain(gain, 0)
        self.rtlsdr_source_0.set_if_gain(gain, 0)
        self.rtlsdr_source_0.set_bb_gain(gain, 0)
        self.rtlsdr_source_0.set_antenna('', 0)
        self.rtlsdr_source_0.set_bandwidth(0, 0)
        self.fft_vxx_0_0 = fft.fft_vcc(fftsize, True, window.blackmanharris(fftsize), True, 1)
        self.fft_vxx_0 = fft.fft_vcc(fftsize, True, window.blackmanharris(fftsize), True, 1)
        self.blocks_stream_to_vector_0_0 = blocks.stream_to_vector(gr.sizeof_gr_complex*1, fftsize)
        self.blocks_stream_to_vector_0 = blocks.stream_to_vector(gr.sizeof_gr_complex*1, fftsize)
        self.blocks_multiply_const_xx_0_0_0 = blocks.multiply_const_ff(1/(fftsize*second), fftsize)
        self.blocks_multiply_const_xx_0_0 = blocks.multiply_const_ff(1/(fftsize*second), fftsize)
        self.blocks_integrate_xx_0_0 = blocks.integrate_ff(second, fftsize)
        self.blocks_integrate_xx_0 = blocks.integrate_ff(second, fftsize)
        self.blocks_head_0_0 = blocks.head(gr.sizeof_gr_complex*1, int(samp_rate*obs_time))
        self.blocks_head_0 = blocks.head(gr.sizeof_gr_complex*1, int(samp_rate*obs_time))
        self.blocks_file_sink_0_0 = blocks.file_sink(gr.sizeof_float*fftsize, file_name1, False)
        self.blocks_file_sink_0_0.set_unbuffered(False)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_float*fftsize, file_name2, False)
        self.blocks_file_sink_0.set_unbuffered(False)
        self.blocks_correctiq_0_0 = blocks.correctiq()
        self.blocks_correctiq_0 = blocks.correctiq()
        self.blocks_complex_to_mag_0_0 = blocks.complex_to_mag(fftsize)
        self.blocks_complex_to_mag_0 = blocks.complex_to_mag(fftsize)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_complex_to_mag_0, 0), (self.blocks_integrate_xx_0, 0))
        self.connect((self.blocks_complex_to_mag_0_0, 0), (self.blocks_integrate_xx_0_0, 0))
        self.connect((self.blocks_correctiq_0, 0), (self.blocks_stream_to_vector_0, 0))
        self.connect((self.blocks_correctiq_0_0, 0), (self.blocks_stream_to_vector_0_0, 0))
        self.connect((self.blocks_head_0, 0), (self.blocks_correctiq_0, 0))
        self.connect((self.blocks_head_0_0, 0), (self.blocks_correctiq_0_0, 0))
        self.connect((self.blocks_integrate_xx_0, 0), (self.blocks_multiply_const_xx_0_0, 0))
        self.connect((self.blocks_integrate_xx_0_0, 0), (self.blocks_multiply_const_xx_0_0_0, 0))
        self.connect((self.blocks_multiply_const_xx_0_0, 0), (self.blocks_file_sink_0, 0))
        self.connect((self.blocks_multiply_const_xx_0_0_0, 0), (self.blocks_file_sink_0_0, 0))
        self.connect((self.blocks_stream_to_vector_0, 0), (self.fft_vxx_0, 0))
        self.connect((self.blocks_stream_to_vector_0_0, 0), (self.fft_vxx_0_0, 0))
        self.connect((self.fft_vxx_0, 0), (self.blocks_complex_to_mag_0, 0))
        self.connect((self.fft_vxx_0_0, 0), (self.blocks_complex_to_mag_0_0, 0))
        self.connect((self.rtlsdr_source_0, 0), (self.blocks_head_0, 0))
        self.connect((self.rtlsdr_source_0_0, 0), (self.blocks_head_0_0, 0))


    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.set_second(int( self.samp_rate/self.fftsize ))
        self.blocks_head_0.set_length(int(self.samp_rate*self.obs_time))
        self.blocks_head_0_0.set_length(int(self.samp_rate*self.obs_time))
        self.rtlsdr_source_0.set_sample_rate(self.samp_rate)
        self.rtlsdr_source_0_0.set_sample_rate(self.samp_rate)

    def get_fftsize(self):
        return self.fftsize

    def set_fftsize(self, fftsize):
        self.fftsize = fftsize
        self.set_second(int( self.samp_rate/self.fftsize ))
        self.blocks_multiply_const_xx_0_0.set_k(1/(self.fftsize*self.second))
        self.blocks_multiply_const_xx_0_0_0.set_k(1/(self.fftsize*self.second))

    def get_second(self):
        return self.second

    def set_second(self, second):
        self.second = second
        self.blocks_multiply_const_xx_0_0.set_k(1/(self.fftsize*self.second))
        self.blocks_multiply_const_xx_0_0_0.set_k(1/(self.fftsize*self.second))

    def get_obs_time(self):
        return self.obs_time

    def set_obs_time(self, obs_time):
        self.obs_time = obs_time
        self.blocks_head_0.set_length(int(self.samp_rate*self.obs_time))
        self.blocks_head_0_0.set_length(int(self.samp_rate*self.obs_time))

    def get_gain(self):
        return self.gain

    def set_gain(self, gain):
        self.gain = gain
        self.rtlsdr_source_0.set_gain(self.gain, 0)
        self.rtlsdr_source_0.set_if_gain(self.gain, 0)
        self.rtlsdr_source_0.set_bb_gain(self.gain, 0)
        self.rtlsdr_source_0_0.set_gain(self.gain, 0)
        self.rtlsdr_source_0_0.set_if_gain(self.gain, 0)
        self.rtlsdr_source_0_0.set_bb_gain(self.gain, 0)

    def get_cfreq(self):
        return self.cfreq

    def set_cfreq(self, cfreq):
        self.cfreq = cfreq
        self.rtlsdr_source_0.set_center_freq(self.cfreq, 0)
        self.rtlsdr_source_0_0.set_center_freq(self.cfreq, 0)




def main(top_block_cls=RTLSDR, options=None):
    tb = top_block_cls()

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        sys.exit(0)

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    tb.start()

    tb.wait()


if __name__ == '__main__':
    main()
