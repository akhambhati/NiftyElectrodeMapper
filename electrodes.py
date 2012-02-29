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
import random



class CTElectrode(vtk.vtkLODActor):
    """
    Actor class for rendering CT-based electrodes
    """
    def __init__(self, electrode_data):

        def __pca(A):
            """
            Principal Component Analaysis
            """
            # Get Dimensions
            num_data, dim = A.shape

            # Center data
            mean_A = A.mean(axis=0)
            for i in range(num_data):
                A[i] -= mean_A

            u,s,v = np.linalg.svd(A)

            normPos = s.argmin()
            return v[:, normPos]

        """
        Setup Surface Rendering
        """

        # Apply a discrete marching cubes algorithm to extract segmented
        # surface contours with incremental values
        self.electrodeExtractor = vtk.vtkDiscreteMarchingCubes()
        self.electrodeExtractor.SetInput(electrode_data.GetOutput())
        self.electrodeExtractor.GenerateValues(1,\
                np.min(electrode_data.GetArray()),\
                np.max(electrode_data.GetArray()))

        self.electrodeMapper = vtk.vtkPolyDataMapper()
        self.electrodeMapper.SetInputConnection(\
                self.electrodeExtractor.GetOutputPort())
        self.electrodeMapper.ScalarVisibilityOff()

        self.electrodeProperty = vtk.vtkProperty()
        self.electrodeProperty.SetColor(1.0, 0.5, 0.0)

        self.SetMapper(self.electrodeMapper)
        self.SetProperty(self.electrodeProperty)
        self.electrodeExtractor.Update()

        self.grid = vtk.vtkAssembly()
        self.electrodePolyData = vtk.vtkPolyData()
        electrodePoints = vtk.vtkPoints()
        ePointMat = np.mat(np.zeros(shape = (36, 3)))
        chanIdx = 0
        for segLabel in np.unique(electrode_data.GetArray()):
            # This is hacked for MAYO34 in order to show just the 6x6 grid
            sphere = vtk.vtkSphereSource()
            sphere.SetRadius(2)

            sphereMap = vtk.vtkPolyDataMapper()
            sphereMap.SetInput(sphere.GetOutput())

            if (segLabel > 0 and segLabel < 37):
                x,y,z = np.nonzero(electrode_data.GetArray() == segLabel)
                nx = 0.9375*np.mean(z)
                ny = 0.9375*np.mean(y)
                nz = 1.5*np.mean(x)

                ePointMat[chanIdx] = [nx, ny, nz]
                chanIdx = chanIdx + 1
        normalVector = __pca(ePointMat.copy())

        for chan in range(chanIdx):
            xx = ePointMat[chan, 0]
            yy = ePointMat[chan, 1]
            zz = ePointMat[chan, 2]

            sphereActor = vtk.vtkOpenGLActor()
            sphereActor.SetMapper(sphereMap)
            sphereActor.GetProperty().SetColor(1.0, 0.0, 1.0)
            sphereActor.SetPosition(xx - 0 * normalVector[0],\
                                    yy - 0 * normalVector[1],\
                                    zz - 0 * normalVector[2])

            electrodePoints.InsertNextPoint(xx - 0 * normalVector[0],\
                                            yy - 0 * normalVector[1],\
                                            zz - 0 * normalVector[2])
            self.grid.AddPart(sphereActor)

        self.electrodePolyData.SetPoints(electrodePoints)
        self.GridSurface()

    def GridSurface(self):
        #----------------------------------------------------------------------
        # ELECTRODE CONTACT DELAUNAY TRIANGULATION ----------------------------
        deln = vtk.vtkDelaunay3D()
        deln.SetInput(self.electrodePolyData)
        deln.SetTolerance(0.01)
        tmapper = vtk.vtkTextureMapToSphere()
        tmapper.SetInputConnection(deln.GetOutputPort())
        #tmapper.PreventSeamOn()
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputConnection(tmapper.GetOutputPort())

        # TEST TEXTURE PART
        atext = vtk.vtkOpenGLTexture()

        self.contactScalarData = vtk.vtkDoubleArray()
        self.contactScalarData.SetNumberOfComponents(1)
        for tupleID in range(self.electrodePolyData.GetNumberOfPoints()):
            self.contactScalarData.InsertNextValue(random.uniform(0, 1))

        self.UpdateGridSurface()
        atext.SetInput(self.electrodePolyData)
        self.triangulation = vtk.vtkOpenGLActor()
        self.triangulation.SetMapper(mapper)
        self.triangulation.SetTexture(atext)
        self.triangulation.GetProperty().SetOpacity(1)

    def UpdateGridSurface(self):
        for tupleID in range(self.contactScalarData.GetNumberOfTuples()):
            self.contactScalarData.SetTuple1(tupleID, random.uniform(0,1))
        self.electrodePolyData.GetPointData().SetScalars(self.contactScalarData)
        self.electrodePolyData.Modified()


class MappedElectrode:
    """
    Class for mapping and rendering mapped electrode
    """
    global sphereRadius
    sphereRadius = 1

    def __init__(self):
        """
        Set a particular cursor for mapping channels for this electrode, and
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
