#!/usr/bin/env python
"""
A tool for plotting and mapping representative electrode channels on cortical
surface renderings.
"""


# Scientific Libraries
import math
import numpy as np

# I/O Handling
from Volume.Formats import NIfTI
import sys
import curses

#Visualization Libraries
from vtk import *
try:
    from wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor
except ImportError:
    print 'Warning: If you are getting flickering get \
            wxVTKRenderWindowInteractor from: \
            "http://www.siafoo.net/snippet/312"'
    from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor
try:
    from vtkImageArray import vtkImageImportFromArray
except ImportError:
    print 'You need to get the updated vtkImageImportFromArray from:\
            "http://www.siafoo.net/snippet/313" '
    exit(1)
from cortex import Cortex
from electrodes import MappedElectrode
from electrodes import CTElectrode

# GUI Libraries
import wx


# Simulation Libraries
import feature_extraction as RTSim
import time


class ElectrodeMappingInteractor(wxVTKRenderWindowInteractor):
    def __init__(self, parent, brain_data,\
            simData, elec_map_fname, elec_ct_data=None):

        #----------------------------------------------------------------------
        # Renderer and Interactor
        wxVTKRenderWindowInteractor.__init__(self, parent, -1,\
                size=parent.GetSize())
        ren = vtk.vtkOpenGLRenderer()
        self.GetRenderWindow().AddRenderer(ren)

        #----------------------------------------------------------------------
        # Setup Surface Rendering
        myPatient = Cortex(brain_data)
        ren.AddViewProp(myPatient)

        #----------------------------------------------------------------------
        # Setup Simulation Module
        simStream = simData.SlidingWindow()

        #----------------------------------------------------------------------
        # Check and render segmented electrode CT surface
        print elec_ct_data
        if elec_ct_data is not None:
            myElectrodeCT = CTElectrode(elec_ct_data)
            #ren.AddViewProp(myElectrodeCT)
            #ren.AddViewProp(myElectrodeCT.grid)
            ren.AddViewProp(myElectrodeCT.triangulation)
            #ren.AddViewProp(myElectrodeCT.cbr)
            print "No electrode CT data specified"

        # Use the trackball camera for interaction
        style = vtk.vtkInteractorStyleTrackballCamera()
        self.SetInteractorStyle(style)

        #----------------------------------------------------------------------
        # Cortical Surface Picking
        posPicker = vtk.vtkCellPicker()
        myPatient.SurfacePicker()
        posPicker.AddLocator(myPatient.cortexLocator)

        # Electrode Actor Picker
        electPicker = vtk.vtkPropPicker()

        # Electrode class instance, there should be one instance for each
        # electrode that will be mapped onto a particular cortical surface
        # TODO: This needs to be made modular, addition of multiple electrodes
        newElectrode = MappedElectrode()

        # Add a electrode placement cursor to the window
        ren.AddViewProp(newElectrode.channelCursor)

        def MoveCursor(wxVTKRenderWindowInteractor, events=""):
            # Function for replacing mouse cursor with channel cursor (sphere)
            self.GetRenderWindow().HideCursor()
            x, y = wxVTKRenderWindowInteractor.GetEventPosition()
            posPicker.Pick(x, y, 0, ren)
            p = posPicker.GetPickPosition()
            # only actors in the channelActors collection are pickable
            electPicker.PickProp(x, y, ren, newElectrode.channelActors)
            # picker returns None when an actor not in channelActors is picked
            if electPicker.GetActor() is not None:
                newElectrode.UpdateChannelCursor(p[0], p[1], p[2],\
                    deleteCursor = 1)
            else:
                newElectrode.UpdateChannelCursor(p[0], p[1], p[2],\
                    deleteCursor = 0)
            wxVTKRenderWindowInteractor.Render()

        def middleClickMouse(wxVTKRenderWindowInteractor, events=""):
            # Function for adding new channel markers to the render/collection
            x, y = wxVTKRenderWindowInteractor.GetEventPosition()
            posPicker.Pick(x, y, 0, ren)
            p = posPicker.GetPickPosition()
            ren.AddViewProp(\
                    newElectrode.AddChannelRepresentationActor(\
                    p[0], p[1], p[2]))
            wxVTKRenderWindowInteractor.Render()

        def rightClickMouse(wxVTKRenderWindowInteractor, events=""):
            # Function for deleting channel markers from render/collection
            x, y = wxVTKRenderWindowInteractor.GetEventPosition()
            electPicker.PickProp(x, y, ren, newElectrode.channelActors)
            if electPicker.GetActor() is not None:
                newElectrode.RemoveChannelRepresentationActor(\
                        electPicker.GetActor())
                ren.RemoveActor(electPicker.GetActor())
            wxVTKRenderWindowInteractor.Render()

        # RenderWindowInteractor event observers
        self.AddObserver("MouseMoveEvent", MoveCursor)
        self.AddObserver("MiddleButtonPressEvent", middleClickMouse)
        self.AddObserver("RightButtonPressEvent", rightClickMouse)

        #----------------------------------------------------------------------
        # wx Event Handling
        def OnKeyPress(event):
            # TODO: Add key handling, disable RenderWindowInteractor shortcuts
            key_code = event.GetKeyCode()

            # Press 'S' to save configuration
            if key_code == ord('S') or key_code == ord('s'):
                newElectrode.SaveConfiguration(elec_map_fname)

            # Press 'O' to open/load configuration
            if key_code == ord('O') or key_code == ord('o'):
                newElectrode.channelActors.InitTraversal()
                for idx in range(newElectrode.channelActors.GetNumberOfItems()):
                    nextActor = newElectrode.channelActors.GetNextItem()
                    ren.RemoveActor(nextActor)

                newElectrode.LoadConfiguration(elec_map_fname)

                newElectrode.channelActors.InitTraversal()
                for idx in range(newElectrode.channelActors.GetNumberOfItems()):
                    nextActor = newElectrode.channelActors.GetNextItem()
                    ren.AddViewProp(nextActor)

            # Press 'UP ARROW' to increase opacity of cortex surface
            if key_code == wx.WXK_UP:
                myPatient.SetOpacityUp()

            # Press 'DOWN ARROW' to decrease opacity of cortex surface
            if key_code == wx.WXK_DOWN:
                myPatient.SetOpacityDown()

            # Press 'RIGHT ARROW' to start the real time simulation
            if key_code == wx.WXK_RIGHT:
                startTime = time.clock()
                tic = time.clock()
                toc = time.clock()
                while(1):
                    time.sleep(0.02)
                    tic = time.clock()
                    myElectrodeCT.UpdateGridSurface(simStream.next())
                    self.Render()
                    toc = time.clock()
                    #print "Update Time: %3.2f msec || Elapsed Time: %3f sec" % (100*(toc-tic), time.clock()-startTime)

        parent.Bind(wx.EVT_KEY_UP, OnKeyPress)

if __name__ == '__main__':
    app = wx.App()
    frame = wx.Frame(None, -1, 'NiftyElectrodeMapping (NEM) Tool',\
            wx.DefaultPosition, wx.Size(600, 600))

    # Load the specified brain NIfTI file
    try:
        brain_data_filename = sys.argv[1]
        print brain_data_filename
        raw_data = NIfTI.ReadFile(brain_data_filename).get_data()
        brain_data = vtkImageImportFromArray()
        brain_data.SetDataSpacing((0.9375,  0.9375,  1.5))
        brain_data.SetArray(raw_data.transpose())
    except:
        print 'Could not import brain NIfTI!'
        exit(1)

    # Get electrode mapping filename
    try:
        electrode_map_fname = sys.argv[2]
    except:
        electrode_map_fname = 'backup.csv'
        print 'No mapping configuration file supplied, saving to backup.csv'

    # Load the specified electrode-CT NIfTI file
    try:
        electrode_data_filename = sys.argv[3]
        print electrode_data_filename
        raw_data = NIfTI.ReadFile(electrode_data_filename).get_data()
        elect_data = vtkImageImportFromArray()
        elect_data.SetDataSpacing((0.9375, 0.9375, 1.5))
        elect_data.SetArray(raw_data.transpose())
    except:
        elect_data = None
        print 'Could not import electrode NIfTI!'


    # Real time EEG simulatio
    fname = "/home/akhambhati/Litt_Fourier/ankk/Neuralynx_RT Project/Mayo Data/NEO_SZ_01/NEO_SZ_01_Grid_slice_003.mat"

    """
    20th order Cauer filter, 100 - 500 Hz bandpass, 
    65dB minimum lower/upper stopband attenuation
    0.5dB maximum pass band ripple
    25 Hz lower/upper transition width
    """
    import scipy.signal as sigp
    filt_num, filt_den = sigp.iirdesign(\
            wp = [100./(2713./2), 500./(2713./2)],\
            ws= [75./(2713./2), 525./(2713./2)],\
            gstop= 60, gpass=0.5, ftype='ellip')
    HFOFeat = lambda winSignal: np.sqrt(np.mean(\
            sigp.filtfilt(filt_num, filt_den, winSignal)**2))
    LLFeat = lambda winSignal: np.sqrt(np.mean(np.diff(winSignal)**2))


    temporalGrid = RTSim.EEGElectrode(fname,\
            winSize=0.1, winShift=0.05, trainTime=60,\
            feat=LLFeat)


    canvas = ElectrodeMappingInteractor(frame,\
                                            brain_data,\
                                            temporalGrid,\
                                            electrode_map_fname,\
                                            elect_data)

    frame.Show()
    app.MainLoop()
