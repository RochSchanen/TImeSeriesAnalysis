#!/usr/bin/python3
# file: ielib.py
# author: Roch Schanen
# date: 2024 06 21
# content:
# repository:

#####################################################################
#                                                          PACKAGES #
#####################################################################

# numpy         from "https://numpy.org/"

#####################################################################
#                                                           IMPORTS #
#####################################################################

def displayFileSize(fp):
    from os import stat
    fs = stat(fp).st_size
    print(f'File Size is {fs:.0f} Bytes.')
    print(f'File Size is {fs/1024:.0f} KB.')
    print(f'File Size is {fs/1024/1024:.0f} MB.')
    return fs

def importCSV(fp, skip = 4):

    from numpy import loadtxt
    data = loadtxt(fp, delimiter=',', skiprows = skip)
    # 12s for 15055399 data points

    # from numpy import genfromtxt
    # data = genfromtxt(fp, delimiter=',', skip_header = skip)
    # 36s for 15055399 data points

    print(f"loaded {data[:,0].size} data points.")    

    from sys import getsizeof
    print(f"memory usage = {getsizeof(data)} Bytes.")

    return data
