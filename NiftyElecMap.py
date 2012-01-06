#!/usr/bin/env python
"""
A tool for plotting and mapping representative electrode channels on cortical
surface renderings.
"""


# Scientific Libraries
import math
import numpy

# File Type Handling
from Volume.Formats import NIfTI
import sys

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
from electrodes import Electrode
import electrodes

# GUI Libraries
import wx


class ElectrodeMappingInteractor(wxVTKRenderWindowInteractor):
    def __init__(self, parent, brain_data):

        #----------------------------------------------------------------------
        # Renderer and Interactor
        wxVTKRenderWindowInteractor.__init__(self, parent, -1,\
                size=parent.GetSize())
        ren = vtk.vtkRenderer()
        self.GetRenderWindow().AddRenderer(ren)

        #----------------------------------------------------------------------
        # Setup Surface Rendering
        myPatient = Cortex(brain_data)

        #----------------------------------------------------------------------
        # Render cortical surface
        ren.AddViewProp(myPatient)

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
        newElectrode = Electrode()

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

if __name__ == '__main__':
    app = wx.App()
    frame = wx.Frame(None, -1, 'NiftyElectrodeMapping (NEM) Tool',\
            wx.DefaultPosition, wx.Size(600, 600))

    # Load the specified NIfTI file
    try:
        brain_data_filename = sys.argv[1]
        print brain_data_filename
        raw_data = NIfTI.ReadFile(brain_data_filename).get_data()
        brain_data= vtkImageImportFromArray()
        brain_data.SetArray(raw_data)
    except:
        print 'Invalid files!'
        exit(1)

    canvas = ElectrodeMappingInteractor(frame, brain_data)

    frame.Show()
    app.MainLoop()
