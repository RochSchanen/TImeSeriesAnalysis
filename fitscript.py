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
# from "https://numpy.org/"
from numpy import array, log, exp, pi

#####################################################################
#                                                     OPEN DOCUMENT #
#####################################################################

D = Document()
D.opendocument("./results.pdf")
# import results
from numpy import load
DATA = load({
    500 : "./fit_500_100.npz",
    250 : "./fit_250_50.npz",
    100 : "./fit_100_25.npz",
    50  : "./fit_50_20.npz",
    10  : "./fit_10_3.npz",
    }[500])
DATA_TIME = DATA["DATA_TIME"]
DATA_AMPLITUDE = DATA["DATA_AMPLITUDE"]
DATA_PERIOD = DATA["DATA_PERIOD"]

#####################################################################
#                                                PLOT LOG AMPLITUDE #
#####################################################################

X, Y = X, Y = DATA_TIME, log(DATA_AMPLITUDE/DATA_AMPLITUDE[0])

fg, ax = stdFig(
    f"LOG_AMPLITUDE",
    X, "S", "Time",
    Y, "n.u.", r"$\log(A/A_0$)"
    )
ax.set_ylabel(r"$\log\left(\frac{A(t)}{A(0)}\right)$")

# in-phase fitting function
def ff(t, a, w):
    from numpy import exp
    x = t/w
    y = a*exp(-x)
    return y


# find time at half maximmum
from numpy import searchsorted
i = searchsorted(X, 160)

# plot
ax.plot(
    X[i]*fg.scx,
    Y[i]*fg.scy,
    'ko',
    markerfacecolor = 'white',
    markersize = 8
    )

from scipy.stats import linregress
r = linregress(X[i:], Y[i:])
A, B = r.slope, r.intercept
ax.plot(X*fg.scx, (X*A+B)*fg.scy, "-.r")

tmp = r"S$^{-1}$"

headerText(fg, f"""
slope    = {A:.6f}{tmp}
-> tau   = {fEng(-1/A)}S
-> width = {fEng(-A/pi)}Hz
""")

D.exportfigure(f"LOG_AMPLITUDE")

#####################################################################
#                                                    PLOT AMPLITUDE #
#####################################################################

X, Y = DATA_TIME, DATA_AMPLITUDE

fg, ax = stdFig(
    f"AMPLITUDE",
    X, "S", "Time",
    Y, "V", "Amplitude"
    )

# in-phase fitting function
def ff(t, a, w):
    from numpy import exp
    x = t/w
    y = a*exp(-x)
    return y

# setup pre-fit parameters
ax.plot(X*fg.scx, ff(X, Y[i]/exp(A*X[i]), -1.0/A)*fg.scy, "-.r")

headerText(fg, f"""
tau = {fEng(-1/A)}S
amp = {fEng(Y[i]/exp(A*X[i]))}V
""")

D.exportfigure(f"AMPLITUDE")

#####################################################################
#                                             PLOT SMOOTH FREQUENCY #
#####################################################################

X, Y = DATA_TIME, 1.0 / DATA_PERIOD

fg, ax = stdFig(
    f"SMOOTH_FREQUENCY",
    X, "S", "Time", 
    Y, "Hz", "Frequency")

from scipy.signal import savgol_filter
YS = savgol_filter(Y, 100, 2)
ax.plot(X*fg.scx, YS*fg.scy, "-.r")
# YS2 = savgol_filter(Y, 20, 2)
# ax.plot(X*fg.scx, YS2*fg.scy, "--k")
    
D.exportfigure(f"SMOOTH_FREQUENCY")

######################################################################
#                                             FREQUENCY VS AMPLITUDE #
######################################################################

X, Y = DATA_AMPLITUDE, 1.0 / DATA_PERIOD

fg, ax = stdFig(
    f"FREQUENCY_VS_AMPLITUDE",
    X, "V", "Amplitude", 
    Y, "Hz", "Frequency")

from scipy.signal import savgol_filter
YS = savgol_filter(Y, 100, 2)
ax.plot(X*fg.scx, YS*fg.scy, "-.r")
# YS2 = savgol_filter(Y, 20, 2)
# ax.plot(X*fg.scx, YS2*fg.scy, "--k")
    
D.exportfigure(f"FREQUENCY_VS_AMPLITUDE")

#####################################################################
#                                                              DONE #
#####################################################################

D.closedocument()
