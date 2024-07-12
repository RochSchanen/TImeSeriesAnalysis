#!/usr/bin/python3
# file: fit.py
# author: Roch Schanen
# created: 2024 07 10
# content:

# default block size is 1MB
BLOCK_SIZE = 1<<20

# excel file header length (in number of lines)
HEADER_SIZE = 4

# minimum number of cycles for pre-fit estimations
PRE_CYCLES = 2.5

# time window for singal filtering
FILTER_LENGTH  = 0.25 # in one cycle units

# polynomial order for singal filtering
FILTER_ORDER = 2

#############
### USAGE ###
#############

from sys import argv, exit
# all arguments have fixed positions
# all arguments are required
if len(argv)<2:
    print("""    --- usage ---
    > python3 fit.py 
        1) "sour.csv" 
        2) "dest.dat"
        3) "frequency"
        4) "width"
        5) "interval"
    1) source file path
    2) destination file path
    3) reference frequency [Hz]
    4) window size used to fit harmonics [S]
    5) window center intervals [S]
    """)
    exit()

###################
### SOURCE FILE ###
###################

from sys import argv, exit
from os.path import exists as validPath
# check path
if not validPath(argv[1]):
    print(f"file '{argv[1]}' not found.")
    print(f"exiting...")
    exit()
# build input path
fpi = argv[1]

###################
### OUTPUT FILE ###
###################

from os.path import split as splitPath
p, n = splitPath(argv[2])
# coerce path
if p.strip() == "":
    p = f".outputs"
# check path
from os.path import isdir as validDir
if not validDir(p):
    print(f"directory '{p}' not found.")
    print(f"exiting...")
    exit()
# coerce name
if n.strip() == "":
    n = splitPath(argv[1])[1]
# build output filename & path
fpo = "fpo=''"

########################
### OTHER PARAMETERS ###
########################

f_reference = float(argv[3])
w_width     = float(argv[4])
w_interval  = float(argv[5])

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

# block number
nb = 0
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
f_start, f_stop = t_first, t_first + w_width

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

######################
### PRE-FIT VALUES ###
######################

# compute frame length [S]
f_length = PRE_CYCLES / f_reference

from numpy import array
# convert list to numpy array
X, Y = array(T), array(V)

from numpy import searchsorted
# search boundary indices
i = searchsorted(X, 0)
j = searchsorted(X, f_length) + 1

# coerce to boundaries
X, Y = X[i:j], Y[i:j]

from figlib import stdfig
# plot raw data
fg, ax = stdfig(f"FIG_PRE_",
    "Time", "S", X,
    "Signal", "V", Y)

from figlib import stdplot
stdplot(f"FIG_PRE_", X, Y, fg.colors["gray"])

# define fitting function
def ff(t, p, w, h):
    from numpy import pi, sin
    x = 2*pi*(t-p)/w
    y = h*sin(x)
    return y

# compute running average to help find signal zeros using a
# Savitzky-Golay filter: second degree polynomial fit over
# a quarter of the period length. See documentation at
# "https://docs.scipy.org/doc/scipy/reference/...
# ...generated/scipy.signal.savgol_filter.html"
from scipy.signal import savgol_filter
# filter signal to facilitate the 
# pre-fit parameters estimations
YF = savgol_filter(Y,
    # filter length in points:
    int(1.0 / f_reference * FILTER_LENGTH / t_delta),
    # filter polynomial order:
    FILTER_ORDER)

# display filtered signal
stdplot(f"FIG_PRE_", X, YF, "-.w")

# find the down zero crossings
I = YF > 0 # sign as boolean
J = (I[:-1] & ~I[1:]).nonzero()[0]

# diplay zeros
stdplot(f"FIG_PRE_", X[J], YF[J], 'ko',
    markerfacecolor = 'white',
    markersize = 8,
    )

# collect the first two zeros only, 1 & 2.
i1, i2 = J[0], J[1]

# compute pre-fitting parameters:
t1, t2 = T[i1], T[i2]       # time positions [S]
m1, m2 = min(Y), max(Y)     # signal extrema [V]

# centre position, period, and amplitude
P = [(t1+t2)/2.0, t2-t1, (m2-m1)/2.0]

from scipy.optimize import curve_fit as fit
# fit PRE_CYCLES
P, C = fit(ff, X, Y, p0 = P)
stdplot(f"FIG_PRE_", X, ff(X, *P), "-k")

# get parameters:
phase, period, amplitude = P 

# export results to figure
from figlib import headerText, fEng
headerText(fg, f""" 
Phase     = {fEng(phase)}S
Period    = {fEng(period)}S
Amplitude = {fEng(amplitude)}V
""")

from figlib import Document
# create pdf
D = Document()
D.opendocument(".outputs/FIG_PRE_.pdf")
D.exportfigure(f"FIG_PRE_")
D.closedocument()

#######################
### FIT FIRST FRAME ###
#######################

# convert list to numpy array
X, Y = array(T), array(V)

# search boundary indices
i = searchsorted(X, f_start)
j = searchsorted(X, f_stop) + 1

# coerce to boundaries
X, Y = X[i:j], Y[i:j]

fg, ax = stdfig(f"FIG_FRAME0_",
    "Time", "S", X,
    "Signal", "V", Y)

stdplot(f"FIG_FRAME0_", X, Y, fg.colors["gray"])

# fit first frame using pre-fit values
P, C = fit(ff, X, Y, p0 = P)
stdplot(f"FIG_FRAME0_", X, ff(X, *P), "-k")

# get parameters:
phase, period, amplitude = P 

# export results to figure
headerText(fg, f""" 
Phase     = {fEng(phase)}S
Period    = {fEng(period)}S
Amplitude = {fEng(amplitude)}V
""")

# create pdf
D = Document()
D.opendocument(".outputs/FIG_FRAME0_.pdf")
D.exportfigure(f"FIG_FRAME0_")
D.closedocument()
