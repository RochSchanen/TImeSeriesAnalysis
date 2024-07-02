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

# compute 3-decades 
def vBin(value):
    from numpy import floor, log
    value = abs(value)
    C, S = (1E+00, "") if float(value) == 0.0 else {
         0: (1024**(-0),  ""),
        +1: (1024**(-1), "K"),
        +2: (1024**(-2), "M"),
        +3: (1024**(-3), "G"),
        +4: (1024**(-4), "T"),
        +5: (1024**(-5), "P"),
    }[int(floor(log(value)/(10*log(2))))]
    return C, S

def fBin(value, format = ".3f"):
    C, S = vBin(value)
    if S == "" : format = ".0f"
    return f"{value*C:{format}}{S}B"

def displayFileSize(fp):
    from os import stat
    fs = stat(fp).st_size
    print(f'File Size is {fBin(fs)}')
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
    print(f"memory usage = {fBin(getsizeof(data))}")

    return data
