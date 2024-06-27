#!/usr/bin/python3
# file: script.py
# author: Roch Schanen
# date: 2024 06 21
# content:
# repository:

# from package: "https://scipy.org/"
from scipy.optimize import curve_fit as fit

# from local modules
from figlib     import Document, stdFig, headerText, fEng
from ielib      import displayFileSize, importCSV


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
    }["small"]

displayFileSize(fp)
DATA = importCSV(fp)

#####################################################################
#                                                GENERAL PARAMETERS #
#####################################################################

# the time interval must be fixed
TIME_INTERVAL = DATA[1, 0] - DATA[0, 0]

# estimated signal frequency in Hertz
APPROX_FREQ   = 96

# compute an appoximate two and a half time the period length
FRAME_LENGTH = int(2.50 / APPROX_FREQ / TIME_INTERVAL)

print(f"""
interval = {fEng(TIME_INTERVAL)}S 
approximate frequency = {fEng(APPROX_FREQ)}Hz
frame 0 length = {fEng(FRAME_LENGTH)}
""")

#####################################################################
#                                                 PLOT FRAME 0 DATA #
#####################################################################

FRAME_NUMBER = 0

wp, ww = 0, FRAME_LENGTH

# select data
T, X = DATA[wp:wp+ww, 0] - DATA[wp, 0], DATA[wp:wp+ww, 1]

# plot raw data
fg, ax = stdFig(f"FRAME{FRAME_NUMBER}", T, "S", "Time", X, "V", "Signal")

#####################################################################
#                                              FIT FIRST FEW PERIOD #
#####################################################################

# define fitting function
def ff(t, p, w, h):
    from numpy import pi, sin
    x = 2*pi*(t-p)/w
    y = h*sin(x)
    return y

# compute running average to find zeros
from scipy.signal import savgol_filter
# trial and error filters:
Y = savgol_filter(X, 50, 1)
Y = savgol_filter(Y, 50, 1)
Y = savgol_filter(Y, 50, 1)
ax.plot(T*fg.scx, Y*fg.scy, "-.k")

# find down zero crossings
from numpy import where, diff, signbit
I = Y > 0
J = (I[:-1] & ~I[1:]).nonzero()[0]
ax.plot(T[J]*fg.scx, Y[J]*fg.scy,
    'ko', markerfacecolor = 'white', markersize = 8)

# get first two zeros indices only
i1, i2 = J[0], J[1]
print(f"approximate number of points per cycles = {i2-i1+1}")

# compute pre-fitting paramneters:
t1, t2 = T[i1], T[i2]       # time position [S]
m1, m2 = min(Y), max(Y)     # signal extrema [V]

# compute center position (positive zero crossing), period, amplitude
P = [(t1+t2)/2.0, t2-t1, (m2-m1)/2]

# fit first single cycle
P, C = fit(ff, T, X, p0 = P)
ax.plot(T*fg.scx, ff(T, *P)*fg.scy, "-k")

# extract parameters
p, w, h = P

ax.legend(["data", "filtered", "down zero crossings", "fit"])
D.exportfigure(f"FRAME{FRAME_NUMBER}")











D.closedocument()
