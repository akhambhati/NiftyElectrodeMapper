"""
Routines related to the processing of location coordinates of volume data sets
"""


import numpy as np


def centroid(x, y, z):
    """
    Calculate the centroid for a particular set of three-dimensional volume 
    coordinates
    """

    return (np.mean(x), np.mean(y), np.mean(z))
