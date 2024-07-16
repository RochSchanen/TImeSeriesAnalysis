#!/usr/bin/python3
# file: tools.py
# author: Roch Schanen
# created: 2024 07 10
# content:

from sys import exit

_C0, _C1, _LL = [], [], ""

def loadBlock(filepath, blocknumber):
    HEADER_SIZE, BLOCK_SIZE = 0, 1024
    global _C0, _C1, _LL
    if blocknumber == 0:
        _C0, _C1, _LL = [], [], ""
        fh = open(filepath, "rb")
        HEADER_SIZE = 4
    fb = fh.read(BLOCK_SIZE)
    if not fb: return 0
    ft = fb.decode('utf-8')
    ft = f"{ll}{ft}"
    L = ft.split("\n")[HEADER_SIZE:]
    for l in L[:-1]:
        t, v = l.split(",")
        _C0.append(float(t))
        _C1.append(float(v))
    _LL = L[-1]
    return blocknumber+1

_FS, _FE, _FN = 0.0, framelength, 0
def loadFrame(filepath, framestart, framelength, period):
    
    global _FS, _FE, _FN

    if float(framestart) == 0.0:
        _FS, _FE, _FN = 0.0, framelength, 0
        _FN = loadBlock(filepath, _FN)
        if not _FN: return None

    while not _C0[-1] > _FE:
        _FN = loadBlock(filepath, _FN)
        if not _FN: return None

    return data

def framecount():
    return _FN

def fix(value, period):
    i = round(value/period)
    return i*period




# declare storage lists
TIME, PERIOD, AMPLITUDE = [], [], []

# clear obsolete data (clip head)
T, V = T[j_start:], V[j_start:]

# shift frame by an integer number of cycles
f_start += fix(f_interval, period)
f_stop  += fix(f_interval, period)



# setup frame boundaries
f_start = t_first
f_stop  = t_first + fix(f_width, period)


