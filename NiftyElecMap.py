#!/usr/bin/env python
"""
A quick and dirty VTK based viewer for volume data.
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
        # Perform Rendering
        ren.AddViewProp(myPatient)

        # Use the trackball camera for interaction
        style = vtk.vtkInteractorStyleTrackballCamera()
        self.SetInteractorStyle(style)

        #----------------------------------------------------------------------
        # Electrode Channel Representation as Spheres
        sphereSource = vtk.vtkSphereSource()
        sphereRadius = 1 
        sphereSource.SetRadius(sphereRadius)
        sphereSource.SetPhiResolution(20)
        sphereSource.SetThetaResolution(20)

        sphereMapper = vtk.vtkPolyDataMapper()
        sphereMapper.SetInput(sphereSource.GetOutput())
        sphereMapper.ScalarVisibilityOff()

        electChanActor = vtk.vtkActor()
        electChanActor.PickableOff()
        electChanActor.SetMapper(sphereMapper)
        electChanActor.GetProperty().SetColor(0, 1, 0)

        ren.AddViewProp(electChanActor)

        #----------------------------------------------------------------------
        # Electrode Picking
        electrodeCount = 0
        posPicker = vtk.vtkCellPicker()
        myPatient.SurfacePicker()
        posPicker.AddLocator(myPatient.cortexLocator)
        electPicker = vtk.vtkPropPicker()
        electChanActorCollection = vtk.vtkActorCollection()

        def MoveCursor(wxVTKRenderWindowInteractor, events=""):
            self.GetRenderWindow().HideCursor()
            x, y = wxVTKRenderWindowInteractor.GetEventPosition()
            posPicker.Pick(x, y, 0, ren)
            p = posPicker.GetPickPosition()
            electChanActor.SetPosition(p[0] + sphereRadius,\
                                       p[1] + sphereRadius,\
                                       p[2] + sphereRadius)
            electPicker.PickProp(x, y, ren, electChanActorCollection)
            print type(electPicker.GetActor())
            if electPicker.GetActor() is not None:
                electChanActor.GetProperty().SetColor(1, 0, 0)
            else:
                electChanActor.GetProperty().SetColor(0, 1, 0)
            wxVTKRenderWindowInteractor.Render()

        def middleClickMouse(wxVTKRenderWindowInteractor, events=""):
            x, y = wxVTKRenderWindowInteractor.GetEventPosition()
            posPicker.Pick(x, y, 0, ren)
            p = posPicker.GetPickPosition()
            setNewElectrode(p)
            wxVTKRenderWindowInteractor.Render()

        def setNewElectrode(p):
            electChanActor = vtk.vtkActor()
            electChanActor.SetMapper(sphereMapper)
            electChanActor.GetProperty().SetColor(0, 0, 1)
            electChanActor.SetPosition(p[0], p[1], p[2])
            electChanActorCollection.AddItem(electChanActor)
            ren.AddViewProp(electChanActor)

        def rightClickMouse(wxVTKRenderWindowInteractor, events=""):
            x, y = wxVTKRenderWindowInteractor.GetEventPosition()
            electPicker.PickProp(x, y, ren, electChanActorCollection)
            if electPicker.GetActor() is not None:
                ren.RemoveActor(electPicker.GetActor())

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
