#!/usr/bin/env python
"""
VTK-based surface rendering of cortex region
"""


import vtk


class Cortex(vtk.vtkLODActor):
    def __init__(self, brain_data, isoval=30):
        # Setup Surface Rendering

        # Apply a marching cubes algorithm to extract surface contour with
        # isovalue of 30 (can change to adjust proper rendering of tissue
        self.cortexExtractor = vtk.vtkMarchingCubes()
        self.cortexExtractor.SetInput(brain_data.GetOutput())
        self.cortexExtractor.SetValue(0, isoval)
        self.cortexExtractor.ComputeNormalsOn()

        # Laplacian smoothing of surface rendering for aesthetics
        # Adds significant delay to rendering
        """
        self.cortexSmoother = vtk.vtkSmoothPolyDataFilter()
        self.cortexSmoother.SetInput(self.cortexExtractor.GetOutput())
        self.cortexSmoother.SetNumberOfIterations(1)
        """

        # Map/Paint the polydata associated with the surface rendering
        self.cortexMapper = vtk.vtkPolyDataMapper()
        self.cortexMapper.SetInput(self.cortexExtractor.GetOutput())
        self.cortexMapper.ScalarVisibilityOff()

        # Color the cortex (RGB)
        self.cortexProperty = vtk.vtkProperty()
        self.cortexProperty.SetColor(1, 1, 1)

        # Set the actor to adhere to mapped surface and inherit properties
        self.SetMapper(self.cortexMapper)
        self.SetProperty(self.cortexProperty)
        self.cortexExtractor.Update()

    def SurfacePicker(self):
        self.cortexLocator = vtk.vtkCellLocator()
        self.cortexLocator.SetDataSet(self.cortexExtractor.GetOutput())
        self.cortexLocator.LazyEvaluationOn()


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