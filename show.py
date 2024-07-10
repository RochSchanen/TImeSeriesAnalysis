#!/usr/bin/python3
# file: show.py
# author: Roch Schanen
# date: 2024 06 21
# content:
# repository:

# BLOCK_SIZE = 1<<20 # 1MB
BLOCK_SIZE = 1<<10 # 1KB

#################
### FILE NAME ###
#################

fn = [
    r'Recording 3',
    r'Ringdown_1325_600mK',
    f'Ringdown_1325mV_750mK',
    ]

# default (development)
fp = f".data/{fn[0]}.csv"

# file from argument
from sys import argv, exit
from os.path import exists as validPath
if len(argv)>1:
    if not validPath(argv[1]):
        print(f"file '{argv[1]}' not found.")
        print(f"exiting...")
        exit()
    fp = argv[1]

##################
### LAST BLOCK ###
##################

# read last block
fh = open(fp, "rb")
fh.seek(-BLOCK_SIZE, 2)
fb = fh.read()
fh.close()

# convert binary to text
ft = fb.decode('utf-8')

# convert text to values
T, V = [], []
for l in ft.split("\n")[1:]:
    if l == "": break
    t, v = l.split(",")
    T.append(float(t))
    V.append(float(v))

t_last = T[-1]

###################
### FIRST BLOCK ###
###################

# read first block
fh = open(fp, "rb")
fb = fh.read(BLOCK_SIZE)

# convert binary to text
ft = fb.decode('utf-8')

# convert text to values:
# skip first four lines
# the last line is most likely truncated
T, V, L = [], [], ft.split("\n")[4:]
for l in L[:-1]:
    t, v = l.split(",")
    T.append(float(t))
    V.append(float(v))
# buffer last line
ll = L[-1]

t_first = T[0]
t_delta = T[1]-T[0]

#####################
### DISPLAY INFOS ###
#####################

# display file infos
from figlib import fEng
from ielib import getFileSize, fBin
fs = getFileSize(fp)
print(f"""--- file info
file size       = {fBin(fs)}
time start      = {fEng(t_first, "07.3f")}S
time end        = {fEng(t_last, "07.3f")}S
time intervals  = {fEng(t_delta)}S
---""")

nb = 1

###################
### NEXT BLOCKS ###
###################

if len(argv)>2:
    t_length = float(argv[2])
    if len(argv)>3:
        t_first = float(argv[3])

    print(f"time origin = {fEng(t_first)}S")
    print(f"time length = {fEng(t_length)}S")

    # load blocks until time start
    while T[-1] < t_first:
        # release memory
        T, V = [], []
        # read next block
        fb = fh.read(BLOCK_SIZE)
        # convert binary to text
        ft = fb.decode('utf-8')
        # catenate buffer
        ft = f"{ll}{ft}"
        # convert text to values
        T, V, L = [], [], ft.split("\n")
        for l in L[:-1]:
            t, v = l.split(",")
            T.append(float(t))
            V.append(float(v))
        # buffer last line
        ll = L[-1]

    t_last = t_first + t_length

    # load blocks until time stop
    while T[-1] < t_last:
        # increment block number
        nb += 1
        # read next block
        fb = fh.read(BLOCK_SIZE)
        # convert binary to text
        ft = fb.decode('utf-8')
        # catenate buffer
        ft = f"{ll}{ft}"
        # convert text to values
        L = ft.split("\n")
        for l in L[:-1]:
            t, v = l.split(",")
            T.append(float(t))
            V.append(float(v))
        # buffer last line
        ll = L[-1]

    # convert list to numpy array
    from numpy import array, searchsorted
    T, V = array(T), array(V)
    
    # coerce to boundaries
    i = searchsorted(T, t_first)
    j = searchsorted(T, t_last)+1
    T, V = T[i:j], V[i:j]

    # imports done
    fh.close()

    ######################
    ### EXPORT DISPLAY ###
    ######################

    from figlib import Document, stdFig, headerText
    # get output filepath parameters
    from os.path import split as splitPath
    p, n = splitPath(fp)
    from os.path import splitext
    n, e = splitext(n)
    D = Document()
    D.opendocument(f".outputs/{n}.pdf")
    fg, ax = stdFig(f"FIG_SHOW", T, "S", "Time", V, "V", "Signal")
    headerText(fg, f"""
        --- time subset ---
    time start      = {fEng(t_first)}S
    time stop       = {fEng(t_last)}S
    time intervals  = {fEng(t_delta)}S
        --- data infos ---
    data points number = {fEng(len(T))}
    block number = {fEng(nb)}
    block size   = {fEng(BLOCK_SIZE)}
    """)
    D.exportfigure(f"FIG_SHOW")
    D.closedocument()
