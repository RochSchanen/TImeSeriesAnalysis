#!/usr/bin/python3
# file: psd.py
# author: Roch Schanen
# created: 2024 07 10
# content:

# standard libraries
from sys import argv, exit

from os.path import exists as validPath
from os.path import split as splitPath
from os.path import isdir as validDir
from os.path import splitext

from time import time

# common libraries
from numpy import pi
from numpy import sin
from numpy import sqrt
from numpy import square
from numpy import arctan2
from numpy import savez
from numpy import array
from numpy import diff
from numpy import where
from numpy import any as np_any

from scipy.optimize import curve_fit as fit
from scipy.signal import savgol_filter

# local libraries
from tools import *

from figlib import Document
from figlib import stdFigure

#####################################################################
#####################################################################
#####################################################################

# record code duration
time_start = time()

#############
### DEBUG ###
#############

DEBUG_BREAK = 0 # 0 means no break
# debug_frame_list = [1, 10, 30, 100]
debug_frame_list = [1, 10]

#############
### USAGE ###
#############

# all arguments have fixed positions
# all arguments are required
if len(argv)<2:
    print("""    --- usage ---
    > python3 psd.py
        1) "sour.csv": source file path
        2) "dest.dat": destination file path
        3) "frequency": reference frequency [Hz]
        4) "time_constant": low pass filters time constant [S]
        5) "window size" data frame lengths [S]
    """)
    exit()

###################
### SOURCE FILE ###
###################

# check file path validity
if not validPath(argv[1]):
    print(f"file '{argv[1]}' not found.")
    print(f"exiting...")
    exit()
# copy and continue
fpi = argv[1]

print(f"processing '{fpi}'")

###################
### OUTPUT FILE ###
###################

fp, fn = splitPath(argv[2])
# coerce path (use default output directory)
if fp.strip() == "":
    fp = f".outputs"
# check path validity
if not validDir(fp):
    print(f"directory '{fp}' not found.")
    print(f"exiting...")
    exit()
# coerce name (copy input filename)
if fn.strip() == "":
    fn = splitPath(fpi)[1]
# split between file name and file extension
fn, fe = splitext(fn)

print(f"all outputs go to '{fp}\\'")

########################
### OTHER PARAMETERS ###
########################

# reference frequency [Hz]
r_frequency = float(argv[3])
# time constant of the low pass filter [S]
t_constant = float(argv[4])
# frames length (one point computed per frame)
f_width = float(argv[5])
print(f"{f_width}")

############################################
### OPEN DOCUMENT FOR DISPLAYING RESULTS ###
############################################

D = Document()
# use a ".psd." extension for all "psd.py" outputs files
D.opendocument(f"{fp}/{fn}.psd.pdf")

print(f"quick view exported to '{fp}/{fn}.psd.pdf'")

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

print(f"data start     = {fEng(t_first)}S")
print(f"data intervals = {fEng(t_delta)}S")

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
on these first few cycles and used to compute the
initial low pass filter levels and the reference
phase (PSDX maximum, PSDY minimum). This insures 
clean start values for the filters. """

# number of cycles for pre-fit estimations
# (The value of 2.5 insures that at leat two
# zeroes are found within the time interval) 
PRE_CYCLES = 2.5
PRE_CYCLES = 12.5

# window size of the Savitzky-Golay filter
# (the size is given in cycle units: a value
# of 1.0 correspond to one cycle)
FILTER_LENGTH  = 0.25

# polynomial order of the Savitzky-Golay filter
FILTER_ORDER = 2 # was empirically determined

# compute frame length [S]
f_length = PRE_CYCLES / r_frequency

# compute frame boundaries
f_start = t_first
f_stop  = t_first + f_length

# load the frame
X, Y = loadFrame(f_start, f_stop)

#####################
### PRE-FIT PLOTS ###
#####################

# plot raw data
fg = stdFigure(f"FIG_PRE_", "Time", "S", "Signal", "V")
fg.plot(X, Y, fg.color["grey"])

# filter data using a Savitzky-Golay filter
# (the filter insure a clean single zero crossing
# per period, avoiding back crossing due to poor
# signal-to-noise ratio)
YF = savgol_filter(Y,
    # filter length (in points)
    int(1.0 / r_frequency * FILTER_LENGTH / t_delta),
    # filter polynomial order
    FILTER_ORDER)
fg.plot(X, YF, "-.w")

# find the down zero crossings
I = YF > 0 # sign as boolean
J = (I[:-1] & ~I[1:]).nonzero()[0]
fg.plot(X[J], YF[J], 'ko',
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
fg.plot(X, ff(X, *P), "-k")

# FIT PARAMETERS:
phase, period, amplitude = P 

# display text results
fg.header(strip(f"""
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

#################################################################
### THE FOLLOWING FUNCTION EMULATES THE DIGITAL LOCK-IN SR830 ###
#################################################################

def PSD_RMS(PSD_FRQ, PSD_TMC, PSD_NUM, T, V, LPFS = None):
    # - PSD_FRQ is frequency
    # - PSD_TMC is timer constant
    # - PSD_NUM is number of low pass filters:
    # 1 ->  -6dB (wait  5 PSD_TMC)
    # 2 -> -12dB (wait  7 PSD_TMC)
    # 3 -> -18dB (wait  9 PSD_TMC)
    # 4 -> -24dB (wait 10 PSD_TMC)
    # - T is the time vector in units of Seconds
    # - V is the signal vector in units of Volts
    # - use LPFS to initialise the low pass filter levels

    # the time sampling intervals must be fixed
    SAMPLING = T[1]-T[0]
    # imports
    from numpy import pi, sin, cos, zeros, sqrt
    # compute phase vector
    P = 2.0*pi*PSD_FRQ*T
    # compute reference vectors
    SIN, COS = sin(P), cos(P)
    # compute digital low pass filter alpha value
    PSD_ALPH = SAMPLING / (SAMPLING + PSD_TMC)
    # reserve memory for the low pass filter outputs
    VX = zeros((len(T), PSD_NUM))
    VY = zeros((len(T), PSD_NUM))
    # setup initial low pass filter levels (defaults to zeros)
    if LPFS is not None: VX[0, :], VY[0, :] = LPFS
    # sweep through references and the signal vector V
    for i, (s, c, v) in enumerate(zip(SIN, COS, V)):
        # use start up value on the first iteration
        k = 0 if i == 0 else i-1
        # compute and record the first stage detection output
        VX[i, 0] = VX[k, 0] + (s*v - VX[k, 0])*PSD_ALPH
        VY[i, 0] = VY[k, 0] + (c*v - VY[k, 0])*PSD_ALPH
        # compute and record the next low pass filters state
        for j in range(1, PSD_NUM):
            VX[i, j] = VX[k, j] + (VX[i, j-1]-VX[k, j])*PSD_ALPH
            VY[i, j] = VY[k, j] + (VY[i, j-1]-VY[k, j])*PSD_ALPH
    # record low pass filters final levels
    LPFS = VX[-1, :], VY[-1, :]
    # keep last low pass filter data set
    VXN, VYN = sqrt(2)*VX[:,-1], sqrt(2)*VY[:,-1]
    # done: return the last low pass filters data sets
    # with the last low pass filter current levels
    return VXN, VYN, LPFS

#####################################################################
#####################################################################
#####################################################################

####################
### PREPARE LOOP ###
####################

# declare storage lists
TIME, PSDX, PSDY = [], [], []
# setup frame boundaries
f_start = t_first
f_stop  = t_first + f_width
# PSD reference frequency
REF_FREQ = 1.0/(period)
# REF_FREQ = 95.89
# PSD SLOPE (-6dB times the number of low pass filters)
PSD_NUM = 2
# PSD TIME CONSTANT in seconds, adjusted (see PSD_RMS comments).
# PSD_TAU = t_constant / {1:3, 2:5, 3:7, 4:10}[PSD_NUM]
PSD_TAU = t_constant

# INITIAL PSD LEVELS
LPFS = ([amplitude/2.0]*PSD_NUM, [0.0]*PSD_NUM)

print(f"reference frequency = {fEng(REF_FREQ)}Hz")
print(f"time constant       = {fEng(PSD_TAU)}S")
print(f"window size         = {fEng(f_width)}S")

while True:

    #######################
    ### LOAD NEXT FRAME ###
    #######################

    data = loadFrame(f_start, f_stop)
    
    # break when end-of-file is reached
    if data is None: break
    
    # split time and signal channels
    T, V = data

    ########################
    ### APPLY PSD FILTER ###
    ########################

    # compute PSD outputs (coerce phase to zero)
    X, Y, LPFS = PSD_RMS(
        REF_FREQ, PSD_TAU, PSD_NUM,
        T-phase, V, LPFS)

    # split low pass filter channels
    LPFX, LPFY = LPFS

    ############################
    ### EXTRA FRAMES DISPLAY ###
    ############################

    # frame number is one less than frame count
    nf = getFrameCount()-1

    # export listed frames
    if nf in debug_frame_list:

        fg = stdFigure(f"FIG_FRAME{nf}_", "Time", "S", "Signal", "V")
        fg.plot(T, V, fg.color["grey"])
        fg.plot(T, X, "--b")
        fg.plot(T, Y, "--r")
        
        # export results to figure
        fg.header(strip(f"""
        --- FRAME{nf} ---
        X = {fEng(LPFX[-1])}V
        Y = {fEng(LPFY[-1])}V
        """))

        # export
        D.exportfigure(f"FIG_FRAME{nf}_")

    ###################
    ### RECORD DATA ###
    ###################

    # collect data
    TIME.append(T[-1])
    PSDX.append(LPFX[-1])
    PSDY.append(LPFY[-1])

    #########################
    ### PREPARE NEXT LOOP ###
    #########################

    if DEBUG_BREAK:
        if nf > DEBUG_BREAK:
            break

    # shift frame by an exact integer number of cycles.
    # this fix is important as it keeps the coherence
    # of the phase parameter.
    f_start += f_width
    f_stop  += f_width

#####################################################################
#####################################################################
#####################################################################

#########################################
### COMPUTE RESULTS AND QUICK DISPLAY ###
#########################################
t_constant * {1:3, 2:5, 3:7, 4:10}[PSD_NUM]

# skip = round()

# convert lists to numpy arrays
TIME = array(TIME)
PSDX = array(PSDX)
PSDY = array(PSDY)

# compute PSD amplitude
AMPL = sqrt(square(PSDX)+square(PSDY))

# plot amplitudes in [V]
fg = stdFigure(f"FIG_AMPL", "Time", "S", "Amplitude", "V")
fg.plot(TIME, PSDX, "-.", linewidth = 0.5, color = fg.color["blue"])
fg.plot(TIME, PSDY, "-.", linewidth = 0.5, color = fg.color["red"])
fg.plot(TIME, AMPL, linewidth = 1.0, color = fg.color["black"])
D.exportfigure(f"FIG_AMPL")

# compute signal phase from PSD reference in cycles
PHAS = arctan2(PSDY, PSDX)/2.0/pi
# detect phase jumps: only one of the two array I or J is not all False
I = diff(PHAS, prepend = PHAS[0]) > +0.8
J = diff(PHAS, prepend = PHAS[0]) < -0.8
# unwrap phase into a monotonic function
for i in where(I)[0]: PHAS[i:] -= 1.0
for j in where(J)[0]: PHAS[j:] += 1.0

# plot the continuous signal phase in cycle
fg = stdFigure(f"FIG_PHAS", "Time", "S", "Phase", "Cycles")
fg.plot(TIME, PHAS, fg.color["grey"])

if np_any(TIME[I]):
    fg.plot(TIME[I], PHAS[I], 'ko',
        markersize = 8, markerfacecolor = 'white')

if np_any(TIME[J]):
    fg.plot(TIME[J], PHAS[J], 'ko',
        markersize = 8, markerfacecolor = 'white')

sgf_PHAS = savgol_filter(PHAS, 25, 1)
fg.plot(TIME, sgf_PHAS, "--",
    linewidth = 1.5, color = fg.color["purple"])

D.exportfigure(f"FIG_PHAS")

# compute the signal frequency from the reference
# frequency and the phase derivative in Hz
FREQ = REF_FREQ + diff(PHAS) / (TIME[1]-TIME[0])

sgf_FREQ = REF_FREQ + diff(sgf_PHAS) / (TIME[1]-TIME[0])

# plot phase in units of cycles: 1 cycle <=> 360 <=> 2 pi
fg = stdFigure(f"FIG_FREQ", "Time", "S", "Frequency", "Hz")
# fg.plot(TIME[1:], FREQ, color = fg.color["grey"])
fg.plot(TIME[1:], sgf_FREQ, "-", linewidth = 1.0, color = fg.color["purple"])

D.exportfigure(f"FIG_FREQ")

#############################
### SAVE RESULTS ###
#############################

# export results
savez(f"{fp}/{fn}.psd.npz",
    DATA_TIME = TIME,
    DATA_AMPLITUDE = AMPL,
    DATA_PERIOD = FREQ,    
    )

############
### DONE ###
############

# export document
D.closedocument()

# display duration and basic infos
time_end = time()
print(f"process duration = {fEng(time_end-time_start)}S")
print(f"frames processed = {getFrameCount()}")
print(f"blocks loaded = {getBlockCount()}")
