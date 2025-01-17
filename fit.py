#!/usr/bin/python3
# file: fit.py
# author: Roch Schanen
# created: 2024 07 10
# content:

# standard libraries
from time import time

from sys import argv, exit
from sys import argv, exit

from os.path import exists as validPath
from os.path import split as splitPath
from os.path import isdir as validDir
from os.path import splitext

# common libraries
from numpy import pi, sin
from numpy import savez
from numpy import array

from scipy.optimize import curve_fit as fit
from scipy.signal import savgol_filter

# local libraries
from tools import *

from figlib import Document
from figlib import stdfig
from figlib import stdplot
from figlib import headerText
from figlib import fEng

#####################################################################
#####################################################################
#####################################################################

# record code duration
time_start = time()

#############
### DEBUG ###
#############

DEBUG_BREAK = 0 # 0 means no break
debug_frame_list = [1, 10, 30, 100]
debug_frame_list = []

#############
### USAGE ###
#############

# all arguments have fixed positions
# all arguments are required
if len(argv)<2:
    print("""    --- usage ---
    > python3 fit.py
        1) "sour.csv" source file path
        2) "dest.dat" destination file path
        3) "frequency" reference frequency [Hz]
        4) "width" window size used to fit the harmonic function [S]
        5) "interval" window shift intervals [S]
    """)
    exit()

###################
### SOURCE FILE ###
###################

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

fp, fn = splitPath(argv[2])
# coerce path
if fp.strip() == "":
    fp = f".outputs"
# check path
if not validDir(fp):
    print(f"directory '{fp}' not found.")
    print(f"exiting...")
    exit()
# coerce name
if fn.strip() == "":
    fn = splitPath(argv[1])[1]
# split file name
fn, fe = splitext(fn)

########################
### OTHER PARAMETERS ###
########################

# approx. signal frequency:
s_frequency = float(argv[3])
# width of the fitting frames
f_width     = float(argv[4])
# intervals between the frames
f_interval  = float(argv[5])

############################################
### OPEN DOCUMENT FOR DISPLAYING RESULTS ###
############################################

D = Document()
D.opendocument(f"{fp}/{fn}.fit.pdf")

#####################################################################
#####################################################################
#####################################################################

############################
### GET FIRST DATA BLOCK ###
############################

# load first block
loadBlock(fpi)
# get time vector from first block
T = getC0()
# START TIME:
t_first = T[0]
# TIME INTERVAL:
t_delta = T[1]-T[0]

#####################
### PRE-FIT FRAME ###
#####################

""" a few cycles of the signal are loaded first
to guess the initial fitting parameters. The signal
is filtered using a Savitzky-Golay filter to obtain
clean zeros crossing. From the zeros values, we
compute the approximate phase and period of the
signal and the extrema gives us an approximate
value of the amplitude. Then a fit is evaluated
on these first few cycles as check """

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

# compute frame length [S]
f_length = PRE_CYCLES / s_frequency

# compute frame boundaries
f_start = t_first
f_stop  = t_first + f_length

# load the frame
X, Y = loadFrame(f_start, f_stop)

#####################
### PRE-FIT PLOTS ###
#####################

# plot raw data
fg, ax = stdfig(f"FIG_PRE_", "Time", "S", X, "Signal", "V", Y)
stdplot(f"FIG_PRE_", X, Y, fg.colors["grey"])

# filter data using a Savitzky-Golay filter
YF = savgol_filter(Y,
    # filter length (in points)
    int(1.0 / s_frequency * FILTER_LENGTH / t_delta),
    # filter polynomial order
    FILTER_ORDER)
stdplot(f"FIG_PRE_", X, YF, "-.w")

# find the down zero crossings
I = YF > 0 # sign as boolean
J = (I[:-1] & ~I[1:]).nonzero()[0]
stdplot(f"FIG_PRE_", X[J], YF[J], 'ko',
    markerfacecolor = 'white',
    markersize = 8,
    )

# define fitting function (simple harmonic)
def ff(t, p, w, h):
    x = 2*pi*(t-p)/w
    y = h*sin(x)
    return y

# COMPUTE PRE-FITTING PARAMETERS:
i1, i2 = J[0], J[1] # collect the first two zeros
t1, t2 = X[i1], X[i2]   # [S]
m1, m2 = min(Y), max(Y) # [V]
#  [centre position, period, amplitude]
P = [(t1+t2)/2.0, t2-t1, (m2-m1)/2.0]

# fit frame
P, C = fit(ff, X, Y, p0 = P)
stdplot(f"FIG_PRE_", X, ff(X, *P), "-k")

# FIT PARAMETERS:
phase, period, amplitude = P 

# display text results
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

# figure done
D.exportfigure(f"FIG_PRE_")

#####################################################################
#####################################################################
#####################################################################

####################
### PREPARE LOOP ###
####################

# declare storage lists
TIME, PERIOD, AMPLITUDE = [], [], []
# setup frame boundaries
f_start = t_first
f_stop  = t_first + f_width

while True:

    #######################
    ### LOAD NEXT FRAME ###
    #######################

    data = loadFrame(f_start, f_stop)
    if data is None: break
    X, Y = data

    ####################
    ### FILTER FRAME ###
    ####################

    Y = savgol_filter(Y, int(1.0 / s_frequency * 0.25 / t_delta), 1)

    #################
    ### FIT FRAME ###
    #################

    # fit frame
    P, C = fit(ff, X-X[0], Y, p0 = P)

    # get parameters:
    phase, period, amplitude = P 

    ############################
    ### EXTRA FRAMES DISPLAY ###
    ############################

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

    ###################
    ### RECORD DATA ###
    ###################

    # collect data
    TIME.append((X[0]+X[-1])/2.0)
    PERIOD.append(period)
    AMPLITUDE.append(amplitude)


    #########################
    ### PREPARE NEXT LOOP ###
    #########################

    # debug break
    if DEBUG_BREAK:
        if nf > DEBUG_BREAK: break

    # shift frame by an integer number of cycles
    # the fix is important as it allows to keep
    # the initial guess parameter "phase" valid
    f_start += fix(f_interval, period)
    f_stop  += fix(f_interval, period)

#####################################################################
#####################################################################
#####################################################################

#############################
### SAVE RESULTS ###
#############################

# convert list to numpy arrays
TIME = array(TIME)
PERIOD = array(PERIOD)
AMPLITUDE = array(AMPLITUDE)

# export results
savez(f"{fp}/{fn}.fit.npz",
    DATA_TIME = TIME,
    DATA_AMPLITUDE = AMPLITUDE,
    DATA_PERIOD = PERIOD,    
    )

#############################
### QUICK RESULTS DISPLAY ###
#############################

stdfig(f"FIG_A", "Time", "S", TIME, "Amplitude", "V", AMPLITUDE)
stdplot(f"FIG_A", TIME, AMPLITUDE, fg.colors["grey"])
stdfig(f"FIG_F", "Time", "S", TIME, "Frequency", "Hz", 1.0/PERIOD)
stdplot(f"FIG_F", TIME, 1/PERIOD, fg.colors["grey"])
D.exportfigure(f"FIG_A")
D.exportfigure(f"FIG_F")

############
### DONE ###
############

# export document
D.closedocument()

# display duration and basic infos
time_end = time()
print(f"process duration = {fEng(time_end-time_start)}S")
print(f"frame processed = {getFrameCount()}")
print(f"block loaded = {getBlockCount()}")
