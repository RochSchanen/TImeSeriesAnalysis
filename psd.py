#!/usr/bin/python3
# file: psd.py
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
from numpy import pi
from numpy import sin
from numpy import sqrt
from numpy import square
from numpy import arctan2
from numpy import savez
from numpy import array

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
debug_frame_list = [1, 10, 30, 100]
# debug_frame_list = [0, 1, 2, 3, 4, 5]

#############
### USAGE ###
#############

# all arguments have fixed positions
# all arguments are required
if len(argv)<2:
    print("""    --- usage ---
    > python3 psd.py
        1) "sour.csv" source file path
        2) "dest.dat" destination file path
        3) "frequency" reference frequency [Hz]
        4) "time constant" of the PSD filter in [S]
        5) "window size" of the data frames [S]
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
t_constant = float(argv[4])
# intervals between the frames
f_width = float(argv[5])

############################################
### OPEN DOCUMENT FOR DISPLAYING RESULTS ###
############################################

D = Document()
D.opendocument(f"{fp}/{fn}.psd.pdf")

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
fg = stdFigure(f"FIG_PRE_", "Time", "S", "Signal", "V")
fg.plot(X, Y, fg.color["grey"])

# filter data using a Savitzky-Golay filter
YF = savgol_filter(Y,
    # filter length (in points)
    int(1.0 / s_frequency * FILTER_LENGTH / t_delta),
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

#####################################################################
#####################################################################
#####################################################################

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
    # with the last low pass filter levels
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
REF_FREQ = 1/period
# PSD SLOPE (-6dB times the number of low pass filters)
PSD_NUM = 2
# PSD TIME CONSTANT in seconds, adjusted.
PSD_TAU = t_constant / {1:3,2:5,3:7,4:10}[PSD_NUM]
# INITIAL PSD LEVELS
LPFS = ([amplitude/2.0]*PSD_NUM, [0.0]*PSD_NUM)

while True:

    #######################
    ### LOAD NEXT FRAME ###
    #######################

    data = loadFrame(f_start, f_stop)
    if data is None: break
    T, V = data

    ########################
    ### APPLY PSD FILTER ###
    ########################

    # compute PSD outputs (p is the previous fitted phase)
    X, Y, LPFS = PSD_RMS(
        REF_FREQ, PSD_TAU, PSD_NUM,
        T-phase, V, LPFS)

    # get low pass filter channels
    LPFX, LPFY = LPFS

    ############################
    ### EXTRA FRAMES DISPLAY ###
    ############################

    nf = getFrameCount()

    print(nf)

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

    # debug break
    if DEBUG_BREAK:
        if nf > DEBUG_BREAK: break

    # shift frame by an integer number of cycles
    # the fix is important as it allows to keep
    # the initial guess parameter "phase" valid
    f_start += f_width
    f_stop  += f_width

#####################################################################
#####################################################################
#####################################################################

TIME = array(TIME)
PSDX = array(PSDX)
PSDY = array(PSDY)

# compute PSD amplitude
AMPL = sqrt(square(PSDX)+square(PSDY))

# plot amplitudes in [V]
fg = stdFigure(f"FIG_AMPL", "Time", "S", "Amplitude", "V")
fg.plot(TIME, PSDX, fg.color["blue"])
fg.plot(TIME, PSDY, fg.color["red"])
fg.plot(TIME, AMPL, fg.color["black"])
D.exportfigure(f"FIG_AMPL")

# compute PSD phase
PHAS = arctan2(PSDY, PSDX)
# find the down zero crossings
I = PHAS > 0
J = (I[:-1] & ~I[1:]).nonzero()[0]
# unwrap phase (no modulo)
for j in J: PHAS[j+1:] += 2.0*pi
# plot phase in units of cycles: 1 cycle <=> 360
fg = stdFigure(f"FIG_PHAS", "Time", "S", "Phase", "Cycles")
fg.plot(TIME, PHAS/2.0/pi, fg.color["grey"])
D.exportfigure(f"FIG_PHAS")

# #####################################################################
# #                                                   PLOT FREQ SHIFT #
# #####################################################################

# skip = 1

# # compute phase shift in cycles per second, i.e. Hertz
# dPdt = diff(P/2.0/pi) / (T[1]-T[0])

# # plot phase in units of cycles: 1 cycle <=> 360 <=> 2 pi
# fg, ax = stdfig(f"FREQ SHIFT",
#     "Time", "S", T[skip+1:],
#     "Freq. shift", "Hz", dPdt[skip:],
#     )

# stdplot("FREQ SHIFT",
#     T[skip+1:],
#     dPdt[skip:],
#     "-k",
#     )

# D.exportfigure(f"FREQ SHIFT")

# #####################################################################
# #                                                       SIGNAL FREQ #
# #####################################################################

# F = 95.855 + dPdt

# # plot phase in units of cycles: 1 cycle <=> 360 <=> 2 pi
# fg, ax = stdfig(f"FREQ",
#     "Time", "S", T[skip+1:],
#     "Freq. shift", "Hz", F[skip:],
#     )

# stdplot("FREQ",
#     T[skip+1:], 
#     F[skip:], 
#     ".", linewidth = 0.5, color = fg.colors["grey"] 
#     )

# from scipy.signal import savgol_filter
# G = savgol_filter(F, 50, 2)
# stdplot("FREQ",
#     T[skip+1:], 
#     G[skip:], 
#     "--", linewidth = 1.5, color = fg.colors["black"] 
#     )

# D.exportfigure(f"FREQ")

#############################
### SAVE RESULTS ###
#############################

# # convert list to numpy arrays
# TIME = array(TIME)
# PERIOD = array(PERIOD)
# AMPLITUDE = array(AMPLITUDE)

# # export results
# savez(f"{fp}/{fn}.fit.npz",
#     DATA_TIME = TIME,
#     DATA_AMPLITUDE = AMPLITUDE,
#     DATA_PERIOD = PERIOD,    
#     )

#############################
### QUICK RESULTS DISPLAY ###
#############################

# stdfig(f"FIG_A", "Time", "S", TIME, "Amplitude", "V", AMPLITUDE)
# stdplot(f"FIG_A", TIME, AMPLITUDE, fg.color["grey"])

# stdfig(f"FIG_F", "Time", "S", TIME, "Frequency", "Hz", 1.0/PERIOD)
# stdplot(f"FIG_F", TIME, 1/PERIOD, fg.color["grey"])

# D.exportfigure(f"FIG_A")
# D.exportfigure(f"FIG_F")

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
