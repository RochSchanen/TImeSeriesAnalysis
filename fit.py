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
# DEBUG_BREAK = 0 # means no break
DEBUG_BREAK = 5

# list of frame to display
debug_frame_list = [
    1, 2, 3, 4, 5,
    # 0, 1, 10, 100, 200,
    ] if DEBUG else []

debug_frame_list = [1, 10, 30, 100]

##############
### CONFIG ###
##############

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
D.opendocument(f"{fp}/{fn}.fit.pdf")

#######################
### GET FIRST BLOCK ###
#######################

from tools import loadBlock
from tools import loadFrame
from tools import getC0

# load first block
loadBlock(fpi)

# get time vector from first block
T = getC0()

# GET START TIME
t_first = T[0]

# GET TIME INTERVAL
t_delta = T[1]-T[0]

#####################
### GET PRE-FRAME ###
#####################

# compute frame length [S]
f_length = PRE_CYCLES / s_frequency

# compute frame boundaries
f_start = t_first
f_stop  = t_first + f_length

X, Y = loadFrame(f_start, f_stop)

##########################
### GET PRE-FIT VALUES ###
##########################

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
t1, t2 = X[i1], X[i2]       # time positions [S]
m1, m2 = min(Y), max(Y)     # signal extrema [V]

# centre position, period, and amplitude
P = [(t1+t2)/2.0, t2-t1, (m2-m1)/2.0]

########################
### RE-FIT PRE-FRAME ###
########################

from scipy.optimize import curve_fit as fit
P, C = fit(ff, X, Y, p0 = P)
stdplot(f"FIG_PRE_", X, ff(X, *P), "-k")

# get parameters:
phase, period, amplitude = P 

#################################
### EXPORT PRE-FIT PARAMETERS ###
#################################

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

####################
### PREPARE LOOP ###
####################

from tools import getFrameCount

# declare storage lists
TIME, PERIOD, AMPLITUDE = [], [], []

# setup frame boundaries
f_start = t_first
f_stop  = t_first + f_width

while True:

    #################
    ### GET FRAME ###
    #################

    data = loadFrame(f_start, f_stop)
    if data is None: break
    X, Y = data

    #################
    ### FIT FRAME ###
    #################

    YF = savgol_filter(Y, int(1.0 / s_frequency * 0.25 / t_delta), 1)

    # fit frame
    P, C = fit(ff, X-X[0], YF, p0 = P)

    # get parameters:
    phase, period, amplitude = P 

    # debug display
    nf = getFrameCount()
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

    #########################
    ### PREPARE NEXT LOOP ###
    #########################

    # shift frame by an integer number of cycles
    f_start += fix(f_interval, period)
    f_stop  += fix(f_interval, period)
    
######################################
################## DISPLAY RESULTS ###
######################################

from numpy import array

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

from tools import getBlockCount

# export
D.closedocument()

time_end = time()
print(f"process duration = {fEng(time_end-time_start)}S")
print(f"frame processed = {getFrameCount()}")
print(f"block loaded = {getBlockCount()}")
