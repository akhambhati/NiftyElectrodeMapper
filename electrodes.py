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

class Electrode:
    global sphereRadius
    sphereRadius = 1
    def __init__(self):
        #----------------------------------------------------------------------
        # Set the electrode placement cursor
        self.channelCursor = self.__CreateChannelRepresentationActor()
        self.channelCursor.GetProperty().SetColor(0, 1, 0)
        self.channelCursor.PickableOff()

        # Instantiate an collection for positioned channel actors
        self.channelActors = vtk.vtkActorCollection()

    def __CreateChannelRepresentationActor(self, sphereRadius = 1):
        sphereSource = vtk.vtkSphereSource()
        sphereSource.SetRadius(sphereRadius)
        sphereSource.SetPhiResolution(20)
        sphereSource.SetThetaResolution(20)

        sphereMapper = vtk.vtkPolyDataMapper()
        sphereMapper.SetInput(sphereSource.GetOutput())
        sphereMapper.ScalarVisibilityOff()

        sphereActor = vtk.vtkActor()
        sphereActor.SetMapper(sphereMapper)

        return sphereActor


    def UpdateChannelCursor(self, x_c, y_c, z_c, deleteCursor = 0, sphereRadius = 1):
        self.channelCursor.SetPosition(x_c + sphereRadius, \
                                       y_c + sphereRadius, \
                                       z_c + sphereRadius)
        if deleteCursor == 0:
            self.channelCursor.GetProperty().SetColor(0, 1, 0)
        else:
            self.channelCursor.GetProperty().SetColor(1, 0, 0)

    def AddChannelRepresentationActor(self, x, y, z, sphereRadius = 1):
        newChannelActor = self.__CreateChannelRepresentationActor()
        newChannelActor.GetProperty().SetColor(0, 0, 1)
        newChannelActor.SetPosition(x + sphereRadius,\
                                    y + sphereRadius,\
                                    z + sphereRadius)
        self.channelActors.AddItem(newChannelActor)
        return newChannelActor

    def RemoveChannelRepresentationActor(self, removedActor):
        self.channelActors.RemoveItem(removedActor)

    def SetElectrodeColor(self, r, g, b):
        tempProperty = vtk.vtkProperty()
        tempProperty.SetColor(r, g, b)
        self.channelActors.ApplyProperty(tempProperty)

    def SaveConfiguration(self, fname):
        f = open(fname, 'wb')
        try:
            writer = csv.writer(f)
            writer.writerow( ('channelID',\
                              'x_coor',\
                              'y_coor',\
                              'z_coor',\
                              'red_px',\
                              'grn_px',\
                              'blu_px') )

            self.channelActors.InitTraversal()
            for actorIdx in range(self.channelActors.GetNumberOfItems()):
                nextActor = self.channelActors.GetNextActor()
                pos = nextActor.GetPosition()
                col = nextActor.GetProperty().GetColor()
                writer.writerow( (actorIdx,\
                                    pos[0],\
                                    pos[1],\
                                    pos[2],\
                                    col[0],\
                                    col[1],\
                                    col[2]) )
        finally:
            f.close()

    def LoadConfiguration(self, fname):
        f = open(fname, 'rb')
        try:
            reader = csv.reader(f)
            for row in reader:
                print row[0] + ' , ' + row[1] 
        finally:
            f.close()


