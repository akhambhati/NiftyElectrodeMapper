#!/usr/bin/env python
"""
Class for displaying and storing user-specified, visual electrode channel
information
"""


# Visualization Libraries
import vtk

# File I/O Libraries
import sys
import csv

# Scientific Libraries
import numpy as np


class ElectrodeGrid2(vtk.vtkAssembly):
    def __init__(self, channelDims=(6, 6),\
                       imagingDims=(0.9375, 0.9375, 1.5),\
                       channelDiam=4.0,\
                       channelSpacing=10,\
                       channelHeight=2.3,\
                       electrodeHeight=0.7,\
                       electrodeLength=90,\
                       electrodeWidth=80):

        # Set dimensional properties for the electrode instance
        self.imagingDims = imagingDims
        self.channelDims = channelDims
        self.channelDiam = channelDiam
        self.channelSpacing = channelSpacing
        self.channelHeight = channelHeight
        self.electrodeHeight = electrodeHeight
        self.electrodeLength = electrodeLength
        self.electrodeWidth = electrodeWidth

        chanID = 0
        for widthPos in range(self.channelDims[1]):
            for lengthPos in range(self.channelDims[0]):
                x = widthPos * self.channelSpacing / self.imagingDims[0]
                z = lengthPos * self.channelSpacing / self.imagingDims[0]
                y = 0.0
                newChan = self.__channelActor(x, y, z)
                newChan.PickableOff()
                self.AddPart(newChan)
                chanID = chanID + 1


    def __CreateGridActor(self):
        self.electrodePD = vtk.vtkPolyData()
        self.electrodePD.SetPoints(self.channelPoints)

        electrodeMapper = vtk.vtkPolyDataMapper()
        electrodeMapper.SetInput(self.electrodePD)

        self.SetMapper(electrodeMapper)
        self.GetProperty().SetColor(1, 1, 0)

    def __channelActor(self, pos_x, pos_y, pos_z):
        channelCyl = vtk.vtkCylinderSource()
        channelCyl.SetRadius(self.channelDiam / self.imagingDims[0])
        channelCyl.SetHeight(self.channelHeight / self.imagingDims[2])
        channelCyl.SetResolution(30)

        channelCylMapper = vtk.vtkPolyDataMapper()
        channelCylMapper.SetInputConnection(channelCyl.GetOutputPort())

        channelCylActor = vtk.vtkActor()
        channelCylActor.SetMapper(channelCylMapper)
        channelCylActor.GetProperty().SetColor(1.0, 1.0, 0)
        channelCylActor.SetPosition(pos_x, pos_y, pos_z)

        return channelCylActor


    def UpdateChannelCursor(self, x_c, y_c, z_c, deleteCursor=0):
        """
        Updates the position of the placement cursor as the mouse moves within
        the render window
        """
        self.SetPosition(x_c, y_c, z_c)

    def Pitch(self):
        self.RotateX(5)

    def Roll(self):
        self.RotateY(5)

    def Yaw(self):
        self.RotateZ(5)

    def RegisterGrid(self, brain_surface):

        # Create a vtkPoints containing channel center locations
        self.channelPoints = vtk.vtkPoints()
        self.channelPolyData = vtk.vtkPolyData()

        chanID = 0
        for widthPos in range(self.channelDims[1]):
            for lengthPos in range(self.channelDims[0]):
                x = widthPos * self.channelSpacing / self.imagingDims[0]
                z = lengthPos * self.channelSpacing / self.imagingDims[0]
                y = 0.0
                self.channelPoints.InsertPoint(chanID, x, y, z)
                chanID = chanID + 1
        self.channelPolyData.SetPoints(self.channelPoints)

        transMat = self.GetMatrix()

        allAssemblyTransform = vtk.vtkTransform()
        allAssemblyTransform.SetMatrix(transMat)
        allAssemblyTransform.Update()

        allAssemblyTransFilt = vtk.vtkTransformPolyDataFilter()
        allAssemblyTransFilt.SetInput(self.channelPolyData)
        allAssemblyTransFilt.SetTransform(allAssemblyTransform)
        allAssemblyTransFilt.Update()

        currentChannelPolyData = allAssemblyTransFilt.GetOutput()

        #### MAKE SURE TRANSFORMATION ACTUALLY WORKS.
        """
        pd = allAssemblyTransFilt.GetOutput()

        transformAssembly = vtk.vtkAssembly()
        for chanIdx in range(pd.GetNumberOfPoints()):
            point = [0, 0, 0]
            pd.GetPoint(chanIdx, point)
            regChan = self.__channelActor(point[0], point[1], point[2])
            transformAssembly.AddPart(regChan)

        return transformAssembly
        """
        ####

        icp = vtk.vtkIterativeClosestPointTransform()
        icp.SetSource(currentChannelPolyData)
        icp.SetTarget(brain_surface.GetOutput())
        icp.GetLandmarkTransform().SetModeToRigidBody()
        icp.SetMaximumNumberOfIterations(100)
        icp.StartByMatchingCentroidsOn()
        icp.Modified()
        icp.Update()

        print icp.GetMatrix()

        icpTransformFilter = vtk.vtkTransformPolyDataFilter()
        icpTransformFilter.SetInput(self.channelPolyData)
        icpTransformFilter.SetTransform(icp)
        icpTransformFilter.Update()

        transformedSource = icpTransformFilter.GetOutput()
        channelCollectionTransformed = vtk.vtkActorCollection()
        channelAssemblyTransformed = vtk.vtkAssembly()
        for chanIdx in range(transformedSource.GetNumberOfPoints()):
            point = [0, 0, 0]
            transformedSource.GetPoint(chanIdx, point)
            regChan = self.__channelActor(point[0], point[1], point[2])
            channelCollectionTransformed.AddItem(regChan)
            channelAssemblyTransformed.AddPart(regChan)

        return channelAssemblyTransformed


class ElectrodeGrid(vtk.vtkAssembly):
    """
    Actor class for rendering an electrode plane with specified channel dims.
    All electrode property units are in mm.
    """
    def __init__(self, channelDims=(8, 8),\
                       imagingDims=(0.9375, 0.9375, 1.5),\
                       channelDiam=4.0,\
                       channelSpacing=10,\
                       electrodeHeight=0.7,\
                       electrodeLength=90,\
                       electrodeWidth=80):

        # Set dimensional properties for the electrode instance
        self.imagingDims = imagingDims
        self.channelDims = channelDims
        self.channelDiam = channelDiam
        self.channelSpacing = channelSpacing
        self.electrodeHeight = electrodeHeight
        self.electrodeLength = electrodeLength
        self.electrodeWidth = electrodeWidth

        chanCount = 1
        self.channelActors = vtk.vtkActorCollection()
        for widthPos in range(self.channelDims[1]):
            for lengthPos in range(self.channelDims[0]):
                newChan = self.__channelActor()
                newChan.SetPosition(\
                        widthPos * self.channelSpacing / self.imagingDims[0],\
                        0,\
                        lengthPos * self.channelSpacing / self.imagingDims[0]\
                        )
                newChan.PickableOff()
                self.channelActors.AddItem(newChan)
                self.AddPart(newChan)
                chanCount = chanCount + 1
        self.PickableOff()

    def __channelActor(self):
        channelCyl = vtk.vtkCylinderSource()
        channelCyl.SetRadius(self.channelDiam / self.imagingDims[0])
        channelCyl.SetHeight(self.electrodeHeight / self.imagingDims[2])
        channelCyl.SetResolution(30)

        channelCylMapper = vtk.vtkPolyDataMapper()
        channelCylMapper.SetInputConnection(channelCyl.GetOutputPort())

        channelCylActor = vtk.vtkActor()
        channelCylActor.SetMapper(channelCylMapper)
        channelCylActor.GetProperty().SetColor(1.0, 1.0, 0)

        return channelCylActor

    def __electrodeBorders(self):
        cornerTopLeft = (0,\
                         self.electrodeLength / self.imagingDims[0],\
                         0)
        cornerTopRight = (self.electrodeWidth / self.imagingDims[0],\
                          self.electrodeLength / self.imagingDims[0],\
                          0)
        cornerBotRight = (self.electrodeWidth / self.imagingDims[0],\
                          0,\
                          0)
        cornerBotLeft = (0, 0, 0)

        borderPoints = vtk.vtkPoints()
        borderPoints.InsertNextPoint(cornerBotLeft)
        borderPoints.InsertNextPoint(cornerTopLeft)
        borderPoints.InsertNextPoint(cornerTopRight)
        borderPoints.InsertNextPoint(cornerBotRight)
        borderPoints.InsertNextPoint(cornerBotLeft)

        borderLines = vtk.vtkPolyLine()
        borderLines.GetPointIds().SetNumberOfIds(5)
        borderLines.GetPointIds().SetId(0, 0)
        borderLines.GetPointIds().SetId(1, 1)
        borderLines.GetPointIds().SetId(2, 2)
        borderLines.GetPointIds().SetId(3, 3)
        borderLines.GetPointIds().SetId(4, 4)

        cells = vtk.vtkCellArray()
        cells.InsertNextCell(borderLines)

        borderPolyData = vtk.vtkPolyData()
        borderPolyData.SetPoints(borderPoints)
        borderPolyData.SetLines(cells)

        borderMapper = vtk.vtkPolyDataMapper()
        borderMapper.SetInput(borderPolyData)

        borderActor = vtk.vtkActor()
        borderActor.SetMapper(borderMapper)
        borderActor.GetProperty().SetColor(1.0, 0.0, 0.5)

        return borderActor

    def UpdateChannelCursor(self, x_c, y_c, z_c, deleteCursor=0):
        """
        Updates the position of the placement cursor as the mouse moves within
        the render window
        """
        self.SetPosition(x_c, y_c, z_c)

    def Pitch(self):
        self.RotateX(5)

    def Roll(self):
        self.RotateY(5)

    def Yaw(self):
        self.RotateZ(5)


class CTElectrode(vtk.vtkLODActor):
    """
    Actor class for rendering CT-based electrodes
    """
    def __init__(self, electrode_data):
        """
        Setup Surface Rendering
        """

        # Apply a discrete marching cubes algorithm to extract segmented
        # surface contours with incremental values
        self.electrodeExtractor = vtk.vtkDiscreteMarchingCubes()
        self.electrodeExtractor.SetInput(electrode_data.GetOutput())
        self.electrodeExtractor.GenerateValues(1, 0, 42)
        """        self.electrodeExtractor.GenerateValues(1,\
                np.min(electrode_data.GetArray()),\
                np.max(electrode_data.GetArray()))
        """
        self.electrodeMapper = vtk.vtkPolyDataMapper()
        self.electrodeMapper.SetInputConnection(\
                self.electrodeExtractor.GetOutputPort())
        self.electrodeMapper.ScalarVisibilityOff()

        self.electrodeProperty = vtk.vtkProperty()
        self.electrodeProperty.SetColor(1.0, 0.5, 0.0)

        self.SetMapper(self.electrodeMapper)
        self.SetProperty(self.electrodeProperty)
        self.electrodeExtractor.Update()

class MappedElectrode:
    """
    Class for mapping and rendering mapped electrode
    """
    global sphereRadius
    sphereRadius = 1

    def __init__(self):
        """
        Set a particular cursor formapping channels for this electrode, and
        initialize actor collections and actor property lists
        """
        # Set the electrode placement cursor
        self.channelCursor = self.__CreateChannelRepresentationActor()
        self.channelCursor.GetProperty().SetColor(0, 1, 0)
        self.channelCursor.PickableOff()
        # Instantiate an collection for positioned channel actors
        self.channelActors = vtk.vtkActorCollection()
        # Instantiate a dictionary list to hold information mapped channel info
        self.channelInfo = []
        self.channelProperties = ['channelID',\
                                  'x_coor',\
                                  'y_coor',\
                                  'z_coor',\
                                  'red_px',\
                                  'grn_px',\
                                  'blu_px']

    def __CreateChannelRepresentationActor(self):
        """
        Returns a sphere actor that serves as a channel marker and
        as a placement cursor
        """

        # Create the source, which creates internal array representation
        sphereSource = vtk.vtkSphereSource()
        sphereSource.SetRadius(sphereRadius)
        sphereSource.SetPhiResolution(20)
        sphereSource.SetThetaResolution(20)

        # map the source to polygonal data
        sphereMapper = vtk.vtkPolyDataMapper()
        sphereMapper.SetInput(sphereSource.GetOutput())
        sphereMapper.ScalarVisibilityOff()

        # take the polygonal data and create a visual actor for it
        sphereActor = vtk.vtkActor()
        sphereActor.SetMapper(sphereMapper)

        return sphereActor

    def UpdateChannelCursor(self, x_c, y_c, z_c, deleteCursor=0):
        """
        Updates the position of the placement cursor as the mouse moves within
        the render window
        """

        self.channelCursor.SetPosition(x_c + sphereRadius, \
                                       y_c + sphereRadius, \
                                       z_c + sphereRadius)

        # Handles coloring the cursor when the mouse is above an existing chan
        if deleteCursor == 0:
            self.channelCursor.GetProperty().SetColor(0, 1, 0)
        else:
            self.channelCursor.GetProperty().SetColor(1, 0, 0)

    def AddChannelRepresentationActor(self, x_coor,\
                                            y_coor,\
                                            z_coor,\
                                            red_px=0.0,\
                                            grn_px=0.0,\
                                            blu_px=1.0):
        """
        Add a new channel to the electrode configration with spatial position,
        size, and color
        """

        # Initialize the new channel actor
        newChannelActor = self.__CreateChannelRepresentationActor()
        newChannelActor.GetProperty().SetColor(red_px, grn_px, blu_px)
        newChannelActor.SetPosition(x_coor + sphereRadius,\
                                    y_coor + sphereRadius,\
                                    z_coor + sphereRadius)

        # Add actors to running collection of actors in the scene
        self.channelActors.AddItem(newChannelActor)
        return newChannelActor

    def RemoveChannelRepresentationActor(self, removedActor):
        """
        Remove a channel from the collection of channel actors
        """
        self.channelActors.RemoveItem(removedActor)

    def SetElectrodeColor(self, r, g, b):
        """
        Set the color for a channel actor
        """
        tempProperty = vtk.vtkProperty()
        tempProperty.SetColor(r, g, b)
        self.channelActors.ApplyProperty(tempProperty)

    def __ParseChannelProperties(self):
        """
        Populate a list of dictionaries containing properties of all the
        channel actors in the current configuration
        """

        # Initialize the channel info list that will contain dictionaries
        self.channelInfo = []

        # Add an entry to the channelProperties list for keeping track of chans
        self.channelActors.InitTraversal()
        for idx in range(self.channelActors.GetNumberOfItems()):
            nextActor = self.channelActors.GetNextActor()
            pos = nextActor.GetPosition()
            col = nextActor.GetProperty().GetColor()

            # Add a new dictionary to the channelInfo list
            self.channelInfo.append({})
            validActor = True
            for prop in self.channelProperties:
                if prop == 'channelID':
                    self.channelInfo[-1][prop] = len(self.channelInfo) - 1
                elif prop == 'x_coor':
                    self.channelInfo[-1][prop] = pos[0] - sphereRadius
                elif prop == 'y_coor':
                    self.channelInfo[-1][prop] = pos[1] - sphereRadius
                elif prop == 'z_coor':
                    self.channelInfo[-1][prop] = pos[2] - sphereRadius
                elif prop == 'red_px':
                    self.channelInfo[-1][prop] = col[0]
                elif prop == 'grn_px':
                    self.channelInfo[-1][prop] = col[1]
                elif prop == 'blu_px':
                    self.channelInfo[-1][prop] = col[2]
                else:
                    validActor = False

            # Failed parses are popped from the list of dictionaries
            if validActor is False:
                self.channelInfo.pop()

    def __GroupChannelAdd(self):
        """
        Populate the electrode's actor collection with channel properties from
        list-dictionary
        """

        # Reset the actor collection for adding the new actors.
        self.channelActors = vtk.vtkActorCollection()

        # Run through each channel and add a new actor for each loaded chan
        # Info was stored as strings, needs to be converted back to float
        for chan in self.channelInfo:
            nextLoadActor = self.AddChannelRepresentationActor(\
                    float(chan[self.channelProperties[1]]),\
                    float(chan[self.channelProperties[2]]),\
                    float(chan[self.channelProperties[3]]),\
                    float(chan[self.channelProperties[4]]),\
                    float(chan[self.channelProperties[5]]),\
                    float(chan[self.channelProperties[6]]))

    def SaveConfiguration(self, fname):
        f = open(fname, 'wb')
        try:
            writer = csv.writer(f)
            writer.writerow(self.channelProperties)
            self.__ParseChannelProperties()
            for chan in self.channelInfo:
                writeList = []
                for i, j in enumerate(self.channelProperties):
                    writeList.append(chan[self.channelProperties[i]])
                writer.writerow(writeList)
            print "Configuration saved to " + fname
        finally:
            f.close()

    def LoadConfiguration(self, fname):
        f = open(fname, 'rb')
        try:
            reader = csv.reader(f)
            self.channelInfo = []
            firstRec = True
            for fields in reader:
                if firstRec:
                    fieldNames = fields
                    if sorted(fieldNames) == sorted(self.channelProperties):
                        firstRec = False
                        print "CSV validity check passed."
                else:
                    self.channelInfo.append({})
                    for i, j in enumerate(fields):
                        print i, j
                        self.channelInfo[-1][fieldNames[i]] = j
            if firstRec:
                print "CSV validity check failed."
            else:
                self.__GroupChannelAdd()
        finally:
            f.close()
