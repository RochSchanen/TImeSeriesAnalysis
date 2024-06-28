#!/usr/bin/python3
# file: script.py
# author: Roch Schanen
# date: 2024 06 21
# content:
# repository:

# from package: "https://scipy.org/"
from scipy.optimize import curve_fit as fit

# from local modules
from figlib import Document, stdFig, headerText, fEng
from ielib import displayFileSize, importCSV
from numpy import zeros, mean

#####################################################################
#                                                     OPEN DOCUMENT #
#####################################################################

D = Document()
D.opendocument("./script.pdf")

#####################################################################
#                                              IMPORT FULL DATA SET #
#####################################################################

fp = {
    "large"     : "./Recording 3.csv",
    "small"     : "./Recording 0.csv",
    "medium"    : "./Recording 1.csv",
    }["large"]
displayFileSize(fp)
DATA = importCSV(fp)
# maximum vector length
DATA_PTS = DATA[:,0].size

#####################################################################
#                                                GENERAL PARAMETERS #
#####################################################################

# the time interval must be fixed
TIME_INTERVAL = DATA[1, 0] - DATA[0, 0]
# estimated signal frequency in Hertz
APPROX_FREQ   = 96
# compute an approximate two and a half time the period length
FRAME_LENGTH = int(2.50 / APPROX_FREQ / TIME_INTERVAL)
# display
print(f"""
interval = {fEng(TIME_INTERVAL)}S 
approximate frequency = {fEng(APPROX_FREQ)}Hz
frame 0 length = {fEng(FRAME_LENGTH)}
""")

# initialise frame index
FRAME_NUMBER = 0

#####################################################################
#                                                FRAME 0 PARAMETERS #
#####################################################################

wp, ww = 0, FRAME_LENGTH

# select data
T, X = DATA[wp:wp+ww, 0] - DATA[wp, 0], DATA[wp:wp+ww, 1]

#####################################################################
#                                                      PLOT FRAME 0 #
#####################################################################

# plot raw data
fg, ax = stdFig(
    f"FRAME{FRAME_NUMBER}",
    T, "S", "Time", 
    X, "V", "Signal")

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
ax.plot(T*fg.scx, Y*fg.scy, ":k")
# find the down zero crossings
I = Y > 0
J = (I[:-1] & ~I[1:]).nonzero()[0]
ax.plot(T[J]*fg.scx, Y[J]*fg.scy,
    'ko', markerfacecolor = 'white', markersize = 8)
# collect the first two zeros only
i1, i2 = J[0], J[1]
# estimate the numbers of cylces in the data set
DATA_CYCLES = DATA_PTS // (i2-i1)
# compute pre-fitting parameters:
t1, t2 = T[i1], T[i2]       # time positions [S]
m1, m2 = min(Y), max(Y)     # signal extrema [V]
# centre position, period, and amplitude
P = [(t1+t2)/2.0, t2-t1, (m2-m1)/2]
# fit one single cycle
P, C = fit(ff, T, X, p0 = P)
ax.plot(T*fg.scx, ff(T, *P)*fg.scy, "-k")

#####################################################################
#                                                    EXPORT FRAME 0 #
#####################################################################

# extract parameters
p, w, h = P
# print header
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
header()
# plot legends
ax.legend(["data", "filtered", "down zero crossings", "fit"])
# export results
D.exportfigure(f"FRAME{FRAME_NUMBER}")

#####################################################################
#                                                   LOOP PARAMETERS #
#####################################################################

# use 30 cycles to fit the signal (tweak-able)
FRAME_LENGTH = int(250 * w / TIME_INTERVAL)

# shift frame by half its length (tweak-able)
# the shift MUST BE a integer multiple of cycles (see comment below)
FRAME_SHIFT = 50

# approximate number of frames
FRAME_NUMBERS = int((DATA_PTS - FRAME_LENGTH)/FRAME_SHIFT)
print(f"number of frames = {FRAME_NUMBERS}")

# reserve memory
DATA_TIME = zeros(FRAME_NUMBERS)
DATA_PERIOD = zeros(FRAME_NUMBERS)
DATA_AMPLITUDE = zeros(FRAME_NUMBERS)

# pages to display
PAGES = [1, 10, 100]

while True:

    #################################################################
    #                                       UPDATE FRAME PARAMETERS #
    #################################################################

    # fixed frame length
    ww = FRAME_LENGTH
    # shift by five periods
    wp += int(FRAME_SHIFT * w / TIME_INTERVAL)
    # a motivation to shift by an integer multiple of cycles is that
    # it allows to use the previous fit results as the starting
    # parameters for the next one.

    # safe break (no enough data)
    if wp+ww > DATA_PTS : break

    # increment frame number
    FRAME_NUMBER += 1

    # select data
    T, X = DATA[wp:wp+ww, 0] - DATA[wp, 0], DATA[wp:wp+ww, 1]

    #################################################################
    #                                                    PLOT FRAME #
    #################################################################

    if FRAME_NUMBER in PAGES:
        # plot raw data
        fg, ax = stdFig(
            f"FRAME{FRAME_NUMBER}",
            T, "S", "Time", 
            X, "V", "Signal")

    #################################################################
    #                                                     FIT FRAME #
    #################################################################

    # fit frame
    P, C = fit(ff, T, X, p0 = P)
    if FRAME_NUMBER in PAGES:
        ax.plot(T*fg.scx, ff(T, *P)*fg.scy, "-k")

    #################################################################
    #                                                  EXPORT FRAME #
    #################################################################

    # extract parameters and export
    p, w, h = P
    if FRAME_NUMBER in PAGES:
        header()
        # plot legends
        ax.legend(["data", "fit"])
        # export results
        D.exportfigure(f"FRAME{FRAME_NUMBER}")

    #################################################################
    #                                                RECORD RESULTS #
    #################################################################

    # average time of the data fitted in seconds
    DATA_TIME[FRAME_NUMBER-1] = DATA[wp, 0] + mean(T)
    # period in seconds
    DATA_PERIOD[FRAME_NUMBER-1] = w
    # signal in volts
    DATA_AMPLITUDE[FRAME_NUMBER-1] = h

    #################################################################
    #                                                      END LOOP #
    #################################################################

    if FRAME_NUMBER == FRAME_NUMBERS: break

# export results
from numpy import savez
savez("./results.npz",
    DATA_TIME = DATA_TIME[:FRAME_NUMBER],
    DATA_AMPLITUDE = DATA_AMPLITUDE[:FRAME_NUMBER],
    DATA_PERIOD = DATA_PERIOD[:FRAME_NUMBER],    
    )

#####################################################################
#                                                              DONE #
#####################################################################

D.closedocument()
