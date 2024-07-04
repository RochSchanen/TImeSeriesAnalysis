#!/usr/bin/python3
# file: lockinfit.py
# author: Roch Schanen
# date: 2024 06 21
# content:
# repository:

# results storage
from numpy import square, sqrt, log, pi, arctan2, diff
# figure exports
from figlib import Document, stdfig, stdplot, headerText, fEng

# file imports
# from ielib import displayFileSize, importCSV

# phase sensitive detection
# from psdlib import PSD_RMS

# from package: "https://scipy.org/"
# from scipy.optimize import curve_fit as fit

#####################################################################
#                                                     OPEN DOCUMENT #
#####################################################################

D = Document()
D.opendocument("./psdfit_results.pdf")

#####################################################################
#                                              IMPORT FULL DATA SET #
#####################################################################

from numpy import load

data = load("./psd_results.npz")
T = data["DATA_TIME"]
X = data["DATA_X"]
Y = data["DATA_Y"]

A = sqrt(square(X)+square(Y))

logA = log(A/A[0])

#####################################################################
#                                                    PLOT AMPLITUDE #
#####################################################################

fg, ax = stdfig(
    f"AMPL",
    "Time", "S", T,
    "PSD", "V", X, Y, A
    )

stdplot("AMPL", T, X, "-b", linewidth = 0.5)
stdplot("AMPL", T, Y, "-r", linewidth = 0.5)
stdplot("AMPL", T, A, "--k")

D.exportfigure(f"AMPL")

#####################################################################
#                                                 LOGPLOT AMPLITUDE #
#####################################################################

def ff(t, a, w):
    from numpy import exp
    x = t/w
    y = a*exp(-x)
    return y

# select time index from seconds
from numpy import searchsorted
i = searchsorted(T, 160.0)

from scipy.stats import linregress
r = linregress(T[i:], logA[i:])
A, B = r.slope, r.intercept
fitA = (T*A+B)

fg, ax = stdfig(
    f"LOG",
    "Time", "S", T,
    r"$\log\left(\frac{A(t)}{A(0)}\right)$", "",
    logA, fitA
    )

stdplot("LOG", T, logA, "-k", linewidth = 0.5)
stdplot("LOG", T, fitA, "--r")
stdplot("LOG", T[i], logA[i], 'ko',
    markerfacecolor = 'white', markersize = 6)

tmp = r"S$^{-1}$"
headerText(fg, f"""
slope    = {A:.6f}{tmp}
-> tau   = {fEng(-1/A)}S
-> width = {fEng(-A/pi)}Hz
""")

D.exportfigure(f"LOG")

#####################################################################
#                                                        PLOT PHASE #
#####################################################################

# compute PSD phase
P = arctan2(Y, X)

# find the down zero crossings and unwrap phase
I = P > 0
J = (I[:-1] & ~I[1:]).nonzero()[0]
for j in J: P[j+1:] += 2.0*pi

# plot phase in units of cycles: 1 cycle <=> 360 <=> 2 pi
fg, ax = stdfig(f"PHAS", "Time", "S", T, "Phase", "cycles", P/2.0/pi)
stdplot("PHAS", T, P/2.0/pi, "-k")
D.exportfigure(f"PHAS")

#####################################################################
#                                                   PLOT FREQ SHIFT #
#####################################################################

skip = 5

# compute phase shift in cycles per second, i.e. Hertz
dPdt = diff(P/2.0/pi) / (T[1]-T[0])

# plot phase in units of cycles: 1 cycle <=> 360 <=> 2 pi
fg, ax = stdfig(f"FREQ SHIFT",
    "Time", "S", T[skip+1:],
    "Freq. shift", "Hz", dPdt[skip:],
    )

stdplot("FREQ SHIFT",
    T[skip+1:],
    dPdt[skip:],
    "-k",
    )

D.exportfigure(f"FREQ SHIFT")

#####################################################################
#                                                       SIGNAL FREQ #
#####################################################################

# compute phase shift in cycles per second, i.e. Hertz
F = 95.855 + dPdt

# plot phase in units of cycles: 1 cycle <=> 360 <=> 2 pi
fg, ax = stdfig(f"FREQ",
    "Time", "S", T[skip+1:],
    "Freq. shift", "Hz", F[skip:],
    )

stdplot("FREQ",
    T[skip+1:], 
    F[skip:], 
    "-k",
    )

D.exportfigure(f"FREQ")

#####################################################################
#                                                              DONE #
#####################################################################

D.closedocument()
