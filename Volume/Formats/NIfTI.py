'''
A module for handling NIfTI data types using the Nibabel package. 

Contains methods for Reading and Writing NIfTI files from a particular 
array format.
'''

def ReadFile(data_path):
    '''
    Import a CT with electrode channels that have been individually
    segmented, such that a particular label corresponds to a single
    channel
    '''

    if data_path.endswith('nii') or data_path.endswith('nii.gz'):
        try:
         import nibabel as nib
        except ImportError:
            print "Please install PyNIfTI to load CT data of this type"
        try:
            nim = nib.load(data_path)
            return nim
        except ImportError:
            print "Could not read the specified NIfTI file"
    else:
        print "Please input data in NIfTI (*.nii/*.nii.gz) format"
