#!/usr/bin/python3
# file: tools.py
# author: Roch Schanen
# created: 2024 07 10
# content:

#####################################################################

# global variables for "loadBlock"
_FH = None
_C0, _C1, = [], []
_LL, _NB  = "", 0

def loadBlock(  FILEPATH = "",       # normally empty
                BLOCK_SIZE  = 1<<20, # default is 1MB
                HEADER_SIZE = 4,     # default is 4 lines
        ): # load one block from file

    # _C0 and _C1 are the current storage
    # _NB the number of blocks already read
    # _LL the last, most likely partial, string
    # _FH the fiile handle
    global _FH, _C0, _C1, _NB, _LL

    # no header as default
    hsz = 0
    
    if FILEPATH: # RESET FILE
        if _FH: _FH.close()
        _FH = open(FILEPATH, "rb")
        # RESET global variables
        _C0, _C1 = [], []
        _LL, _NB = "", 0
        # skip file header
        hsz = HEADER_SIZE
    
    # get the block
    fb = _FH.read(BLOCK_SIZE)
    if not fb: # no more block available
        return False
    
    # convert to text
    ft = fb.decode('utf-8')
    # catenate with previous string part
    ft = f"{_LL}{ft}"

    # scan through text and collect numerics    
    L = ft.split("\n")[hsz:]
    for l in L[:-1]:
        t, v = l.split(",")
        _C0.append(float(t))
        _C1.append(float(v))
    
    # buffer last string part
    _LL = L[-1]

    # done
    _NB += 1
    return True

def getBlockCount():
    return _NB

#####################################################################

# global variables for "loadFrame"
_NF = 0

def loadFrame(  framestart, # start of frame in S
                framestop,  # frame width in S
                FILEPATH = "",  # normally empty
        ): # load one frame from file
    
    # _NF frame number
    global _NF, _C0, _C1

    if FILEPATH: # RESET FILE
        loadBlock(FILEPATH)
        # RESET global variables
        _NF = 0

    # load blocks until the end of frame
    # is reached or no more blocks is available
    while not _C0[-1] > framestop:
        if not loadBlock():
            return None

    from numpy import array
    # convert to numpy array
    C = array(_C0)
    
    from numpy import searchsorted
    # find frame span indices
    js = searchsorted(C, framestart)
    je = searchsorted(C, framestop)

    # convert and coerce data to request
    C0 = array(_C0[js:je])
    C1 = array(_C1[js:je])

    # clear obsolete buffer parts
    _C0, _C1 = _C0[js:], _C1[js:]

    # done
    _NF += 1
    return C0, C1

def getFrameCount():
    return _NF

#####################################################################

def getC0():
    from numpy import array
    return array(_C0)

def getC1():
    from numpy import array
    return array(_C1)

#####################################################################

def closeFile():
    _FH.close()