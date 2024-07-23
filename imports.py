#!/usr/bin/python3
# file = imports.py
# author = Roch Schanen
# created = 2024 07 22
# content = class & function to import data from files

###############################################################################

""" The dataBlockStream class splits large size file (typically 1Gb or more) in
a stream of smaller size blocks of data. Each block is converted into a set of
text lines. Each lines is scanned converted into data. The text line to data
convertion is provided by the user. see below: scanf_torsion() and
scanf_levitation(). """

class dataBlockStream():

    def __init__(self,
            data,                   # data storage list
            scanf,                  # scanning function
            filepath,               # path to the data file
            blocksize  = 1<<20,     # [bytes] default is 1MB
            headersize = 0,         # [lines]
            ):
        # copy references
        self.dt = data   # reference to data storage
        self.sf = scanf  # reference to scan function
        # record parameters
        self.bs = blocksize
        self.hs = headersize
        # initialise locals
        self.fh = None
        # initialise file
        self.resetstream(filepath)
        # done
        return        

    def resetstream(self, fp = None):
        # initialise locals
        self.ll = ""    # last line buffer
        self.nb = 0     # number of blocks
        # instanciate file
        if fp: self.fh = open(fp, "rb")
        # rewind file pointer
        self.fh.seek(0)
        # done
        return        

    def loadblock(self):
        # select header size
        hs = 0 if self.nb else self.hs
        # get new block
        fb = self.fh.read(self.bs)
        # no more block available
        if not fb:
            return False
        # convert to text and catenate last line
        text = f"{self.ll}{fb.decode('utf-8')}"
        # split into a list of text lines
        lines = text.split("\n")[hs:]
        # scan through lines (skip last)
        for l in lines[:-1]:
            self.sf(self.dt, l)
        # buffer last string
        self.ll = lines[-1]
        # if last block, add last string
        if len(fb) < self.bs:
            # check for empty line
            # (if the file ends with char '\n')
            if self.ll:
                # last convertion
                self.sf(self.dt, self.ll)
        # one more block loaded
        self.nb += 1
        # done
        return True

    def getblockcount(self):
        return self.nb

    # def seekblock(self, n):
    #     from os import stat
    #     # get file size
    #     fs = stat(fp).st_size
    #     # check block boundary
    #     if self.bs*n > fs:
    #         return False
    #     self.fh.seek(self.bs*n)
    #     # done
    #     return True

###############################################################################

def scanf_torsion(data, line):
    # two columns separated by a coma
    t, v = line.split(',')
    # convert and store data
    data.append([float(t), float(v)])
    # done
    return

def scanf_levitation(data, line):
    # record only the first 3 columns (separated by tabulations)
    t, x, y = line.split('\t')[:3]
    # convert and store data
    data.append([float(t), float(x), float(y)])
    # done
    return

# the user can generate more scan function. They require two input
# parameters: a data storage list and the line to scan. The function
# should append the new converted data from the string to the storage
# list. It is recommanded to optimize the code for speed.

###############################################################################

""" The dataFrameStream class splits large size file (typically 1Gb or more) in
a stream of frames of data. The frame boundaries are selected from a column
which data must be in increasing order (for example the time of a measurement).
Successive frames are meant to be loaded in a progessive order. Obsolete data
are discarded and new blocks are loaded automatically as new frames are
requested. """

class dataFrameStream():

    def __init__(self,
            data,                   # data storage list
            scanf,                  # scanning function
            filepath,               # path to the data file
            column = 0,             # monotically increasing column number
            blocksize  = 1<<20,     # [bytes] default is 1MB
            headersize = 0,         # [lines]
            ):
        # initialise locals
        self.nf = 0
        # record parameters
        self.cn = column
        # initialise dataBlockStream
        self.db = dataBlockStream(
                data,
                scanf, 
                filepath,
                blocksize,
                headersize,
            )
        # load the first block
        self.db.loadblock()
        # done
        return

    def loadframe(self,
            framestart, 
            framestop,
            ):
        # columns of last data line
        cs = self.db.dt[-1]
        # load blocks until the end of frame is found
        # or no more blocks is available.
        while not cs[self.cn] > framestop:
            # get next available block
            if not self.db.loadblock():
                # last block reached before end of the frame
                return None
            # select last line data
            cs = self.db.dt[-1]
        # convert data to numpy array
        from numpy import array
        D = array(self.dt)
        # find frame span indices
        from numpy import searchsorted
        js = searchsorted(D[:, self.cn], framestart)
        je = searchsorted(D[:, self.cn], framestop)
        # coerce data to request
        D = D[js:je, :]
        # clear obsolete buffer parts
        del self.dt[:js]
        # one more frame loaded
        self.nf += 1
        # done
        return D

###############################################################################

if __name__ == "__main__":
    
    SELECT = {
            1: "TEST DATABLOCKSTREAM",
            2: "TEST DATAFRAMESTREAM",
     }[2]

###############################################################################
    
    if SELECT == "TEST DATAFRAMESTREAM":

        print(SELECT)

        DATA = []

        d = dataFrameStream(
            DATA, 
            scanf_torsion,
            ".data/Recording 0.csv",    # small test file
            column = 0,                 # time column default is 0
            blocksize = 1024,           # try small blocks
            headersize = 4,             # 4 lines to skip
            )

        D = d.loadframe(0.0, 0.010)

        print(D)


###############################################################################

    if SELECT == "TEST DATABLOCKSTREAM":

        print(SELECT)

        DATA = []

        d = dataBlockStream(
            DATA, 
            scanf_torsion,
            ".data/Recording 0.csv",    # small test file
            blocksize = 1024,           # try small blocks
            headersize = 4,             # 4 lines to skip
            )

        while d.loadblock():
            print(f"- load BLOCK{d.getblockcount()}")

        for i, l in enumerate(DATA):
            print(f"{i}: {l[0]:.3E}, {l[1]:.3E}")

###############################################################################
