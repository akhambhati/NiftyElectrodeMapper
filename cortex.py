#!/usr/bin/env python
"""
VTK-based surface rendering of cortex region
"""


import vtk


class Cortex(vtk.vtkLODActor):
    def __init__(self, brain_data, isoval=34):
        # Setup Surface Rendering

        # Gaussian smoothing of surface rendering for aesthetics
        # Adds significant delay to rendering
        self.cortexSmoother = vtk.vtkImageGaussianSmooth()
        self.cortexSmoother.SetDimensionality(3)
        self.cortexSmoother.SetRadiusFactors(0.5, 0.5, 0.5)
        self.cortexSmoother.SetInput(brain_data.GetOutput())

        # Apply a marching cubes algorithm to extract surface contour with
        # isovalue of 30 (can change to adjust proper rendering of tissue
        self.cortexExtractor = vtk.vtkMarchingCubes()
        self.cortexExtractor.SetInput(self.cortexSmoother.GetOutput())
        self.cortexExtractor.SetValue(0, isoval)
        self.cortexExtractor.ComputeNormalsOn()

        # Map/Paint the polydata associated with the surface rendering
        self.cortexMapper = vtk.vtkPolyDataMapper()
        self.cortexMapper.SetInput(self.cortexExtractor.GetOutput())
        self.cortexMapper.ScalarVisibilityOff()

        # Color the cortex (RGB)
        self.cortexProperty = vtk.vtkProperty()
        self.cortexProperty.SetColor(1, 1, 1)
        self.cortexProperty.SetOpacity(1);

        # Set the actor to adhere to mapped surface and inherit properties
        self.SetMapper(self.cortexMapper)
        self.SetProperty(self.cortexProperty)
        self.cortexExtractor.Update()

    def SurfacePicker(self):
        self.cortexLocator = vtk.vtkCellLocator()
        self.cortexLocator.SetDataSet(self.cortexExtractor.GetOutput())
        self.cortexLocator.LazyEvaluationOn()

    def SetOpacityUp(self):
        opc = self.GetProperty().GetOpacity()
        print opc
        if (opc + 0.1) <= 1.0:
            self.GetProperty().SetOpacity(opc + 0.1)
            self.cortexExtractor.Update()

    def SetOpacityDown(self):
        opc = self.GetProperty().GetOpacity()
        print opc
        if (opc - 0.1) >= 0.0:
            self.GetProperty().SetOpacity(opc - 0.1)
            self.cortexExtractor.Update()

if __name__ == '__main__':
    import wx, sys
    from vtkImageArray import vtkImageImportFromArray
    from Volume.Formats import NIfTI
    from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor

    app = wx.App()
    frame = wx.Frame(None, -1, 'Volume Rendering with VTK', wx.DefaultPosition, wx.Size(600, 600))

    try:
        brain_data_filename = sys.argv[1]
        print brain_data_filename
        raw_data = NIfTI.ReadFile(brain_data_filename).get_data()
        brain_data = vtkImageImportFromArray()
        brain_data.SetArray(raw_data)
    except:
        print 'Invalid files!'

    myWXrwc = wxVTKRenderWindowInteractor(frame, -1, size=frame.GetSize())
    ren = vtk.vtkRenderer()
    myWXrwc.GetRenderWindow().AddRenderer(ren)

    canvas = Cortex(brain_data)
    ren.AddViewProp(canvas)

    style = vtk.vtkInteractorStyleTrackballCamera()
    myWXrwc.SetInteractorStyle(style)

    frame.Show()
    app.MainLoop()
