#! /usr/bin/env python
"""
Feature extraction simulator for real-time visualization
"""

# Scientific Libraries
import math
import numpy as np
import scipy.signal as sigp

# I/O Handling
import scipy.io
import sys
import curses


# Utilities
from time import sleep
import time

class EEGElectrode():
    """
    A class for all EEG data associated with a single electrode entity i.e
    grid, strip, microwire
    """
    def __init__(self, data_file, winSize, winShift, trainTime, feat):
        # Assume the input data is a NxM mat file.
        # N is number of channels, M is the number of samples
        raw_mat_data = scipy.io.loadmat(data_file)
        self.data = raw_mat_data['data']
        self.data_prop = {\
                'NumChannels' : self.data.shape[0],\
                'Fs' : 2713,\
                'NumSamples' : self.data.shape[1]}
        self.sim_prop = {\
                'WindowSize' : winSize,\
                'WindowShift' : winShift,\
                'WindowTrain' : trainTime,\
                'Feature' : feat}


    def SlidingWindow(self):
            """Simulates the real-time arrival of data from a DAQ based on
            window size and window shift."""

            # Calculate window sizes and shifts in terms of samples
            winSize = self.sim_prop['WindowSize'] * self.data_prop['Fs']
            winShift = self.sim_prop['WindowShift'] * self.data_prop['Fs']
            winTrain = self.sim_prop['WindowTrain'] * self.data_prop['Fs']

            # Verify the inputs
            try:
                it = iter(self.data)
            except TypeError:
                raise Exception("**ERROR** signal must be iterable.")
            if  winShift> winSize:
                raise Exception("**ERROR** winShift must not be larger than\
                        winSize")

            # Pre-compute number of full chunks, discounting residual items
            numFullChunks = \
                    1 + ((self.data_prop['NumSamples'] - winSize)/winShift)

            # Initialize training variables
            # Running Code:
            # http://subluminal.wordpress.com/2008/07/31/running-standard-deviations/
            # TODO: Create a Dynamic Range Class
            runMean = 0
            runPwrSumMean = 0
            runStd = 0
            runN = 0
            dynamicRange = [0, 0]
            def DNRAdjust(featureVector, dynRange):
                if dynRange[1] < dynRange[0]:
                    raise Exception("**ERROR** Dynamic range should increase")
                for idx,featVal in enumerate(featureVector):
                    if featVal < dynRange[0]:
                        featureVector[idx] = 0
                    elif featVal > dynRange[1]:
                        featureVector[idx] = 1
                    else:
                        featureVector[idx] = \
                            (featVal - dynRange[0])/(dynRange[1] - dynRange[0])
                return featureVector


            # Generate windowed signal chunks
            # Iterate through each window
            for wIdx in range(0, int(numFullChunks*winShift), int(winShift)):
                feat = np.zeros(self.data_prop['NumChannels'])
                emptyFeat = feat.copy()

                # Calculate feature for specified signal window
                for i, j in enumerate(self.data):
                    feat[i] = self.sim_prop['Feature']\
                            (self.data[i, wIdx:wIdx+int(winSize)])

                #print DNRAdjust(feat, dynamicRange)

                # Check if we want to yield empty features or trained features
                if wIdx > winTrain:
                    trainedFeat = DNRAdjust(feat, dynamicRange)
                    #print "Before Adjust:", feat
                    #print "After Adjust:", trainedFeat
                    yield trainedFeat
                else:
                    # Dynamic Range Training using Running Statistics
                    runN += 1
                    runMean += (np.mean(feat) - runMean) / runN
                    runPwrSumMean += (np.mean(feat)**2 - runPwrSumMean)\
                            /runN
                    runStd = np.sqrt(\
                            (runPwrSumMean*runN - runN*runMean**2) /\
                            (runN-1))
                    dynamicRange = [runMean - 8*runStd, runMean + 8*runStd]

                    #print "Running Mean:", runMean
                    #print "Running Power Sum Average:", runPwrSumMean
                    #print "Running Standard Deviation:", runStd
                    #print "Dynamic Range:", dynamicRange

                    yield emptyFeat



if __name__ == '__main__':
    import matplotlib.pyplot as plt

    fname = '/home/akhambhati/Litt_Fourier/ankk/Neuralynx_RT Project/Mayo Data/NEO_SZ_01/NEO_SZ_01_Grid_slice_002.mat'

    # Test signal processing
    """
    20th order Cauer filter, 100 - 500 Hz bandpass, 
    65dB minimum lower/upper stopband attenuation
    0.5dB maximum pass band ripple
    25 Hz lower/upper transition width
    """
    filt_num, filt_den = sigp.iirdesign(\
            wp = [100./(2713./2), 500./(2713./2)],\
            ws= [75./(2713./2), 525./(2713./2)],\
            gstop= 60, gpass=0.5, ftype='ellip')

    HFOFeat = lambda winSignal: np.sqrt(np.mean(\
            sigp.filtfilt(filt_num, filt_den, winSignal)**2))

    # Test time-series plotting
    #fig = plt.figure()
    #l, = plt.plot(temporalGrid.data[0,:])
    #plt.show()


    tG = EEGElectrode(fname,\
            winSize=0.1,winShift=0.0,trainTime = 60,\
            feat=HFOFeat)

    tt = tG.SlidingWindow()

    while(1):
        sleep(0.05)
        print tt.next()
