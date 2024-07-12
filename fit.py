#!/usr/bin/python3
# file: fit.py
# author: Roch Schanen
# created: 2024 07 10
# content:

BLOCK_SIZE = 1<<20 # 1MB

HEADER_SIZE = 4 # file header lines to skip

#############
### USAGE ###
#############

from sys import argv, exit

# all arguments have fixed positions
# all arguments are required

if len(argv)<2:
    print("""    --- usage ---
    > python3 fit.py "source.csv" "destination.dat" "frequency" "span" "resolution"
    - frequency: estimated frequency [Hz] (used to compute the pre-fit parameters).
    - span: the time window used to fit the harmonic signal [S].
    - resolution: the time by which the window is shifted between harmonic fits [S].
    The resolution is recommanded to be similar or less than the value of the span.
    """)
    exit()

###################
### SOURCE FILE ###
###################

from sys import argv, exit
from os.path import exists as validPath

if not validPath(argv[1]):
    print(f"file '{argv[1]}' not found.")
    print(f"exiting...")
    exit()

fpi = argv[1]

###################
### OUTPUT FILE ###
###################

from os.path import split as splitPath
p, n = splitPath(argv[2])

if p.strip() == "":
    p = f".outputs"

from os.path import isdir as validDir
if not validDir(p):
    print(f"directory '{p}' not found.")
    print(f"exiting...")
    exit()

if n.strip() == "":
    n = splitPath(argv[1])[1]

fpo = ""

########################
### OTHER PARAMETERS ###
########################

f0 = float(argv[3]) # approx. frequency
sp = float(argv[4]) # window/frame span
rs = float(argv[5]) # step resolution

##################
### LAST BLOCK ###
##################

# read last block
fh = open(fpi, "rb")
# only load 1KB block: one data line
# is approximately 30 characters. we
# only need the last line.
fh.seek(-1024, 2)
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

# GET LAST DATA TIME VALUE
t_last = T[-1]

######################################
###################### FIRST FRAME ###
######################################

###################
### FIRST BLOCK ###
###################

nb = 0 # block number

# read first block
fh = open(fpi, "rb")
# use default block size
fb = fh.read(BLOCK_SIZE)
# convert binary to text
ft = fb.decode('utf-8')
# convert text to values:
# skip first four lines
# the last line is most likely truncated
T, V, L = [], [], ft.split("\n")[HEADER_SIZE:]
for l in L[:-1]:
    t, v = l.split(",")
    T.append(float(t))
    V.append(float(v))
# buffer the last line
ll = L[-1]

# GET FIRST DATA TIME VALUE
t_first = T[0]

# GET DATA TIME INTERVAL
t_delta = T[1]-T[0]

###################
### NEXT BLOCKS ###
###################

# first frame boundaries
f_start, f_stop = t_first, t_first + sp

# load blocks until frame length is reached
while T[-1] < f_stop:
    # increment block number
    nb += 1
    # read next block
    fb = fh.read(BLOCK_SIZE)
    # convert binary to text
    ft = fb.decode('utf-8')
    # catenate buffered line
    ft = f"{ll}{ft}"
    # convert text to values
    L = ft.split("\n")
    for l in L[:-1]:
        t, v = l.split(",")
        T.append(float(t))
        V.append(float(v))
    # buffer the last line
    ll = L[-1]


# compute an approximate two and a half times the
# period length for the first frame length
f_length = 2.50 / f0
