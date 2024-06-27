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

# compute an appoximate period length
FRAME_LENGTH = int(1.25/APPROX_FREQ/TIME_INTERVAL)

print(f"""
interval = {fEng(TIME_INTERVAL)}S 
approximate frequency = {fEng(APPROX_FREQ)}Hz
frame 0 length = {fEng(FRAME_LENGTH)}
""")

#####################################################################
#                                                      PLOT FRAME 0 #
#####################################################################

FRAME_NUMBER = 0

wp, ww = 0, FRAME_LENGTH

# select data
T, X = DATA[wp:wp+ww, 0] - DATA[wp, 0], DATA[wp:wp+ww, 1]


fg, ax = stdFig(f"FRAME{FRAME_NUMBER}", T, "S", "Time", X, "V", "Signal")

from scipy.signal import savgol_filter
# trial and error filters:
Y = savgol_filter(X, 50, 1)
Y = savgol_filter(Y, 50, 1)
Y = savgol_filter(Y, 50, 1)

# find down zero crossings
from numpy import where, diff, signbit
I = Y > 0
J = (I[:-1] & ~I[1:]).nonzero()[0]

# exit if un-expected number zeros (there shoud be two)
if not len(J) == 2:
    print(f"number of down zero crossing {len(J)} ", end = "")
    if len(J) > nc: print(" larger ", end = "")
    if len(J) < nc: print(" smaller ", end = "")
    print(f"than expected.")
    print(f"exiting.")
    exit()

# get zeros indices
i1, i2 = J

print(f"approximate number of points per cycles = {i2-i1+1}")

# compute pre-fitting paramneters
# parameters are center position (positive zero crossing), period, amplitude
t1, t2 = T[i1], T[i2]
m1, m2 = min(Y), max(Y)
P = [(t1+t2)/2.0, t2-t1, (m2-m1)/2]

#####################################################################
#                                                                   #
#                         FITTING FUNCTION                          #
#                                                                   #
#####################################################################

# in-phase fitting function
def ff(t, p, w, h):
    from numpy import pi, sin
    x = 2*pi*(t-p)/w
    y = h*sin(x)
    return y

# fit first single cycle
P, C = fit(ff, T, X, p0 = P)

# extract parameters
p, w, h = P











makeHeader(fg)
D.exportfigure(f"FRAME{FRAME_NUMBER}")
D.closedocument()






# ax.plot(T, Y, '-.', color = "r")
# ax.plot(T, ff(T, *P1), "-", color = fit_color)
# ax.plot(T[J], Y[J], 'ko', markerfacecolor = 'white', markersize = 8)
# marker='o', markeredgewidth=1.5, markersize=16
# ax.legend(["data", "filtered", "fitted", "down zero crossings"])
