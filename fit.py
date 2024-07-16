#!/usr/bin/python3
# file: fit.py
# author: Roch Schanen
# created: 2024 07 10
# content:

from time import time
# time measurement this script execution
time_start = time()

#############
### DEBUG ###
#############

# debug enabled flag
DEBUG = False
# stop after "DEBUG_BREAK" frames if DEBUG is True
DEBUG_BREAK = 0 # no break
DEBUG_BREAK = 5
# list of frame to display
debug_frame_list = [
    # 1, 2, 3, 4, 5,
    0, 1, 10, 100, 200,
    ] if DEBUG else []

##############
### CONFIG ###
##############

""" To preserve memory usage, the file is
read in successive blocks. The data are
collected and disgarded immediately after
analysis. Therefore, the code can run
through any file size. """

# A default block size of 1MB seems reasonable
BLOCK_SIZE = 1<<20

# excel file header length
# (This is the number of lines to skip before
# the data values become available)
HEADER_SIZE = 4

""" a few cycles of the signal are loaded first
to guess the initial fitting parameters. The signal
is filtered using a Savitzky-Golay filter to obtain
clean zeros crossing. From the zeros values, we
compute the approximate phase and period of the
signal and the extrema gives us an approximate
value of the amplitude. """

# number of cycles for pre-fit estimations
# (The value of 2.5 insures that at leat two
# zeroes are found within the time interval) 
PRE_CYCLES = 2.5

# window size of the Savitzky-Golay filter
# (the size is given in cycle units: a value
# of 1.0 correspond to one cycle)
FILTER_LENGTH  = 0.25

# polynomial order of the Savitzky-Golay filter
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
    5) window centre intervals [S]
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
fp, fn = splitPath(argv[2])
# coerce path
if fp.strip() == "":
    fp = f".outputs"
# check path
from os.path import isdir as validDir
if not validDir(fp):
    print(f"directory '{fp}' not found.")
    print(f"exiting...")
    exit()
# coerce name
if fn.strip() == "":
    fn = splitPath(argv[1])[1]
# split file name
from os.path import splitext
fn, fe = splitext(fn)

########################
### OTHER PARAMETERS ###
########################

# approx. signal frequency:
s_frequency = float(argv[3])

# force duration to be an integer
# multiple of the signal period.
# this is used to fix the phase
# fitting parameter.
def fix(value, period):
    i = round(value/period)
    return i*period

# width of the fitting frames
f_width     = float(argv[4])

# intervals between the frames
f_interval  = float(argv[5])

######################
### EXPORT RESULTS ###
######################

from figlib import Document
D = Document()
if DEBUG:
    D.opendocument(f"{fp}/{fn}.DEBUG_FIGS.pdf")
else:
    D.opendocument(f"{fp}/{fn}.PREFIT.pdf")

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

##################################
###################### FRAME 0 ###
##################################

nf = 0

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

# compute frame length [S]
f_length = PRE_CYCLES / s_frequency

# compute frame boundaries
f_start = t_first
f_stop  = t_first + f_length

# load blocks until frame length is reached
while not T[-1] > f_stop:
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

from numpy import array
# convert list to numpy array
X, Y = array(T), array(V)

from numpy import searchsorted
# search boundary indices
j_length = searchsorted(X, f_length)

# coerce to boundaries (clip tail)
X, Y = X[:j_length], Y[:j_length]

from figlib import stdfig
# plot raw data
fg, ax = stdfig(f"FIG_PRE_",
    "Time", "S", X,
    "Signal", "V", Y)

from figlib import stdplot
stdplot(f"FIG_PRE_", X, Y, fg.colors["grey"])

# define fitting function
def ff(t, p, w, h):
    from numpy import pi, sin
    x = 2*pi*(t-p)/w
    y = h*sin(x)
    return y

from scipy.signal import savgol_filter
# Savitzky-Golay filter
YF = savgol_filter(Y,
    # filter length (in points)
    int(1.0 / s_frequency * FILTER_LENGTH / t_delta),
    # filter polynomial order
    FILTER_ORDER)

# display filtered signal
stdplot(f"FIG_PRE_", X, YF, "-.w")

# find the down zero crossings
I = YF > 0 # sign as boolean
J = (I[:-1] & ~I[1:]).nonzero()[0]

# display zeros
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

# strip multiline string
def strip(mls):
    r = ""
    for l in mls.split("\n"):
        r = f"{r}{l.strip()}\n"
    return r

# export results to figure
from figlib import headerText, fEng
headerText(fg, strip(f"""
    --- INITIAL GUESSES ---
    Phase     = {fEng((t1+t2)/2.0)}S
    Period    = {fEng(t2-t1)}S
    Amplitude = {fEng((m2-m1)/2.0)}V

    --- PRE-FIT ---
    Phase     = {fEng(phase)}S
    Period    = {fEng(period)}S
    Amplitude = {fEng(amplitude)}V
    """))

# export
D.exportfigure(f"FIG_PRE_")

######################################
###################### NEXT FRAMES ###
######################################

####################
### PREPARE LOOP ###
####################

# declare storage lists
TIME, PERIOD, AMPLITUDE = [], [], []

# setup frame boundaries
f_start = t_first
f_stop  = t_first + fix(f_width, period)

while True:

    ###################
    ### NEXT BLOCKS ###
    ###################

    # load blocks until frame length is reached
    while not T[-1] > f_stop:
        # read next block
        fb = fh.read(BLOCK_SIZE)
        # end-of-file reached
        if not fb: break
        # increment block number
        nb += 1
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
    if not fb: break

    #################
    ### FIT FRAME ###
    #################

    # convert list to numpy array
    X, Y = array(T), array(V)

    # search boundary index
    j_start  = searchsorted(X, f_start)
    j_stop   = searchsorted(X, f_stop)

    # coerce to boundaries
    X, Y = X[j_start:j_stop], Y[j_start:j_stop]

    # fit frame
    P, C = fit(ff, X-X[0], Y, p0 = P)

    # get parameters:
    phase, period, amplitude = P 

    # debug display
    if nf in debug_frame_list:
        fg, ax = stdfig(f"FIG_FRAME{nf}_",
            "Time", "S", X,
            "Signal", "V", Y)
        # add raw data
        stdplot(f"FIG_FRAME{nf}_", X, Y, fg.colors["grey"])
        # add fit
        stdplot(f"FIG_FRAME{nf}_", X, ff(X-X[0], *P), "-k")
        # export results to figure
        headerText(fg, strip(f"""
        --- FRAME{nf} ---
        Phase     = {fEng(phase)}S
        Period    = {fEng(period)}S
        Amplitude = {fEng(amplitude)}V
        """))
        # export
        D.exportfigure(f"FIG_FRAME{nf}_")

    # collect data
    TIME.append((X[0]+X[-1])/2.0)
    PERIOD.append(period)
    AMPLITUDE.append(amplitude)

    # debug break
    if DEBUG:
        if DEBUG_BREAK:
            if nf > DEBUG_BREAK: break

    ##########################
    ### PREPARE NEXT FRAME ###
    ##########################

    # update frame number
    nf += 1

    # clear obsolete data (clip head)
    T, V = T[j_start:], V[j_start:]

    # shift frame by an integer number of cycles
    f_start += fix(f_interval, period)
    f_stop  += fix(f_interval, period)

    #################
    ### LOOP ENDS ###
    #################

######################################
################## DISPLAY RESULTS ###
######################################

# convert list to numpy arrays
TIME = array(TIME)
PERIOD = array(PERIOD)
AMPLITUDE = array(AMPLITUDE)

# export results
from numpy import savez
savez(f"{fp}/{fn}.fit.npz",
    DATA_TIME = TIME,
    DATA_AMPLITUDE = AMPLITUDE,
    DATA_PERIOD = PERIOD,    
    )

if DEBUG:
    stdfig(f"FIG_A", "Time", "S", TIME, "Amplitude", "V", AMPLITUDE)
    stdplot(f"FIG_A", TIME, AMPLITUDE, fg.colors["grey"])
    stdfig(f"FIG_F", "Time", "S", TIME, "Frequency", "Hz", 1.0/PERIOD)
    stdplot(f"FIG_F", TIME, 1/PERIOD, fg.colors["grey"])
    D.exportfigure(f"FIG_A")
    D.exportfigure(f"FIG_F")

############
### DONE ###
############

# export
D.closedocument()

time_end = time()
print(f"process duration = {fEng(time_end-time_start)}S")
print(f"frame processed = {nf}")
print(f"block loaded = {nb}")
