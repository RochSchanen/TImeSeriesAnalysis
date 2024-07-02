#!/usr/bin/python3
# file: script.py
# author: Roch Schanen
# date: 2024 06 21
# content:
# repository:

# from package: "https://scipy.org/"
from scipy.optimize import curve_fit as fit
# file imports
from ielib import displayFileSize, importCSV
# figure exports
from figlib import Document, stdfig, stdplot, headerText, fEng
# phase sensitive detection
from psdlib import PSD_RMS
# results storage
from numpy import zeros, array, square, sqrt, linspace, sin, pi

#####################################################################
#                                                     OPEN DOCUMENT #
#####################################################################

D = Document()

D.opendocument("./psd_results.pdf")

#####################################################################
#                                              IMPORT FULL DATA SET #
#####################################################################

debug = False

if debug:

    # define test signal
    TEST_FREQ = 96.0            # FREQUENCY             Hertz
    TEST_PHAS = 30.0            # PHASE                 Degrees
    TEST_AMPL = 0.1*sqrt(2)     # AMPLITUDE (0.1Vrms)   Volts 
    SAMPLING  = 20E-6           # SAMPLING INTERVAL     Seconds
    TEST_LEN  = 30.0            # SIGNAL LENGTH         Seconds
    DATA_PTS  = int(TEST_LEN / SAMPLING) + 1
    # export values (debug)
    print(f"Frequency = {fEng(TEST_FREQ)}Hz")
    print(f"Phase = {fEng(TEST_PHAS)}°")
    print(f"Amplitude = {fEng(TEST_AMPL)}V")
    print(f"Sampling intervals = {fEng(SAMPLING)}S")
    print(f"Number of points{fEng(DATA_PTS)}")
    # time vector
    T = linspace(0.0, TEST_LEN, DATA_PTS)
    # test signal (includes a phase shift)
    V = TEST_AMPL*sin(2.0*pi*TEST_FREQ*T - pi/180.0*TEST_PHAS)
    DATA = zeros((DATA_PTS, 2))
    DATA[:, 0], DATA[:, 1] = T, V

else:

    # SELECT FROM FILE LIST
    fp = {
        "verylarge" : "./Recording 3.csv",
        "small"     : "./Recording 0.csv",
        "medium"    : "./Recording 1.csv",
        "large"     : "./Recording 2.csv",
        }["verylarge"]

    # DISPLAY INFOS
    displayFileSize(fp)

    # LOAD DATA FROM EXCEL FORMAT
    DATA = importCSV(fp)

    # MEASURE VECTOR LENGTHS
    DATA_PTS = DATA[:,0].size

#####################################################################
#                                                GENERAL PARAMETERS #
#####################################################################

# frame counter
FRAME_NUMBER = 0

# the time interval must be fixed
TIME_INTERVAL = DATA[1, 0] - DATA[0, 0]

# estimated signal frequency in Hertz
ESTIMATED_FREQ   = 96.0

# compute two and a half time the period length
FRAME_LENGTH = int(2.50 / ESTIMATED_FREQ / TIME_INTERVAL)

# display infos
print(f"""
interval = {fEng(TIME_INTERVAL)}S 
approximate frequency = {fEng(ESTIMATED_FREQ)}Hz
frame index = {FRAME_NUMBER}
frame {FRAME_NUMBER} length = {fEng(FRAME_LENGTH)}
""")

#####################################################################
#                                       FRAME 0 PARAMETERS AND DATA #
#####################################################################

wp, ww = 0, FRAME_LENGTH

# select data
T, X = DATA[wp:wp+ww, 0] - DATA[wp, 0], DATA[wp:wp+ww, 1]

#####################################################################
#                                                      PLOT FRAME 0 #
#####################################################################

# plot raw data
fg, ax = stdfig(
    f"FRAME{FRAME_NUMBER}",
    "Time", "S", T, 
    "Signal", "V", X,
    )

stdplot("FRAME0", T, X, color = fg.colors["grey"])

#####################################################################
#                                                       FIT FRAME 0 #
#####################################################################

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

# trial and error filters:
Y = savgol_filter(X, int(FRAME_LENGTH / 2.5 / 4), 2)
stdplot("FRAME0", T, Y, ":k")

# find the down zero crossings
I = Y > 0
J = (I[:-1] & ~I[1:]).nonzero()[0]
stdplot("FRAME0", T[J], Y[J], 'ko',
    markerfacecolor = 'white',
    markersize = 8,
    )

# collect the first two zeros only
i1, i2 = J[0], J[1]

# estimate the numbers of cycles in the data set
DATA_CYCLES = DATA_PTS // (i2-i1+1)

# compute pre-fitting parameters:
t1, t2 = T[i1], T[i2]       # time positions [S]
m1, m2 = min(Y), max(Y)     # signal extrema [V]

# centre position, period, and amplitude
P = [(t1+t2)/2.0, t2-t1, (m2-m1)/2]

# fit the two and a half cycle of data
P, C = fit(ff, T, X, p0 = P)
stdplot("FRAME0", T, ff(T, *P), "-k")

#####################################################################
#                                         EXPORT FITTING PARAMETERS #
#####################################################################

# extract parameters: center position p, width w (period), height h
p, w, h = P

# define header
def header():
    from numpy import sqrt
    return headerText(fg, f"""
    FRAME {FRAME_NUMBER}:
    PHASE     = {fEng(p)}S
    PERIOD    = {fEng(w)}S
     -> FREQ  = {fEng(1/w)}Hz
    AMPLITUDE = {fEng(h)}V
     -> pp    = {fEng(2*h)}V,
     -> rms   = {fEng(h/sqrt(2))}V
    """)

# export header
header()

#####################################################################
#                                   FINALISE PLOT AND EXPORT FIGURE #
#####################################################################

# plot legends
ax.legend(["data", "filtered", "down zero crossings", "fit"])

# export result figure
D.exportfigure(f"FRAME{FRAME_NUMBER}")

#####################################################################
#                                                   LOOP PARAMETERS #
#####################################################################

# THE PSD FILTER IS APPLIED ON DATA CHUNKS
FRAME_LENGTH = int(20 * w / TIME_INTERVAL)
# NUMBER OF FRAMES TO PROCESS
FRAME_NUMBERS = int(DATA_PTS/FRAME_LENGTH)
# ESTIMATED REFERENCE FREQUENCY
REF_FREQ = 1 / w
# PSD SLOPE (-6db x number of low pass filters)
PSD_NUM = 2
# PSD time constant in seconds
PSD_TAU = 0.100

# DISPLAY INFO
print(f"""
frame length        = {fEng(FRAME_LENGTH)}
number of frames    = {fEng(FRAME_NUMBERS)}
reference frequency = {fEng(REF_FREQ)}Hz
psd time constant   = {fEng(PSD_TAU)}S
""")

# RESERVE MEMORY
DATA_TIME = zeros(FRAME_NUMBERS)
DATA_X    = zeros(FRAME_NUMBERS)
DATA_Y    = zeros(FRAME_NUMBERS)

# INITIAL PSD LEVELS
# (the intial filter values must be computable from p, w, h)
LPFS = None

# SETUP FIRST FRAME
wp, ww = 0, FRAME_LENGTH

# LOOP THROUGH FRAMES
for FRAME_NUMBER in range(FRAME_NUMBERS):
    # SELECT DATA
    T, V = DATA[wp:wp+ww, 0], DATA[wp:wp+ww, 1]
    # psd filter
    X, Y, LPFS = PSD_RMS(REF_FREQ, PSD_TAU, PSD_NUM, T, V, LPFS)
    # record time stamp
    DATA_TIME[FRAME_NUMBER] = T[-1]
    # split channels 
    LPFX, LPFY = LPFS
    # record last stage
    DATA_X[FRAME_NUMBER], DATA_Y[FRAME_NUMBER] = LPFX[-1], LPFY[-1]
    # UPDATE FRAME PARAMTERS
    wp += FRAME_LENGTH
    # safe break
    if wp+ww > DATA_PTS : break
    # test
    if FRAME_NUMBER > 100: break

#####################################################################
#                                                    EXPORT RESULTS #
#####################################################################

T = DATA_TIME[:FRAME_NUMBER]
X = DATA_X[:FRAME_NUMBER]
Y = DATA_Y[:FRAME_NUMBER]

A = sqrt(square(X)+square(Y))

fg, ax = stdfig(
    f"PSD_RESULTS",
    "Time", "S", T,
    "PSD", "V", X, Y, A
    )

stdplot("PSD_RESULTS", T, X, "-b")
stdplot("PSD_RESULTS", T, Y, "-r")
stdplot("PSD_RESULTS", T, A, "-k")

D.exportfigure(f"PSD_RESULTS")

#####################################################################
#                                                              DONE #
#####################################################################

D.closedocument()

# export results
from numpy import savez
savez("./psd_results.npz",
    DATA_TIME = DATA_TIME[:FRAME_NUMBER],
    DATA_X = DATA_X[:FRAME_NUMBER],
    DATA_Y = DATA_Y[:FRAME_NUMBER],
    )
