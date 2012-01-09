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

    def UpdateChannelCursor(self, x_c, y_c, z_c, deleteCursor = 0):
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
                                            red_px = 0,\
                                            grn_px = 0,\
                                            blu_px = 1):
        """
        Add a new channel to the electrode configuration with spatial position,
        size, and color
        """

        # Initialize the new channel actor
        newChannelActor = self.__CreateChannelRepresentationActor()
        newChannelActor.GetProperty().SetColor(red_px, grn_px, blu_px)
        newChannelActor.SetPosition(x_coor + sphereRadius,\
                                    y_coor + sphereRadius,\
                                    z_coor + sphereRadius)

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
            self.channelActors.AddItem(nextLoadActor)

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
            print "Configuration saved to "  + fname
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
