#!/usr/bin/python3
# file: decayfit.py
# author: Roch Schanen
# date: 2024 06 18
# content: analyse decay signal fitting through window
# repository:

# from standard libraries:
from sys import exit

#####################################################################
#                                                                   #
#                      EXTERNAL PACKAGES                            #
#                                                                   #
#####################################################################

# from package: "https://numpy.org/"
from numpy import mean
from numpy import array
from numpy import zeros
from numpy import linspace
from numpy import ceil, floor
from numpy import log10, exp, log

# from package: "https://matplotlib.org/"
from matplotlib.pyplot import figure
from matplotlib.pyplot import fignum_exists
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import CSS4_COLORS

# from package: "https://scipy.org/"
from scipy.optimize import curve_fit as fit

#####################################################################
#                                                                   #
#              FIGURE, HEADER, FOOTER, TICKS, DOCUMENT              #
#                                                                   #
#####################################################################

def selectfigure(name):
    if not fignum_exists(name):
        # create figure
        fg = figure(name)
        # set A4 paper dimensions
        fg.set_size_inches(8.2677, 11.6929)
        # create square axis
        w, h = array([1, 1 / 1.4143])*0.7
        x, y = (1-w)/2, (1-h)/2
        ax = fg.add_axes([x, y, w, h])
    else:
        # select figure
        # (here the figure can be of any type)
        fg = figure(name)
        # get axes
        ax = fg.get_axes()[0]
    # done
    return fg, ax

def headerText(text, fg):
    w, h = array([1, 1 / 1.4143])*0.7
    x, y = (1-w)/2, (1-h)/2
    # tx = fg.text(x+w/2, 3*y/2+h, text)
    tx = fg.text(x, 3*y/2+h, text)
    tx.set_fontfamily('monospace')
    tx.set_horizontalalignment('left')
    tx.set_verticalalignment('center')
    tx.set_fontsize("small")
    return tx

def footerText(text, fg):
    w, h = array([1, 1 / 1.4143])*0.7
    x, y = (1-w)/2, (1-h)/2
    # tx = fg.text(x+w/2, y/2, text)
    tx = fg.text(x, y/2, text)
    tx.set_fontfamily('monospace')
    tx.set_horizontalalignment('left')
    tx.set_verticalalignment('center')
    tx.set_fontsize("small")
    return tx

def _getTickIntervals(start, stop, ticks):

    ln10 = 2.3025850929940459

    # trial table
    T = [0.010, 0.020, 0.025, 0.050,
         0.100, 0.200, 0.250, 0.500,
         1.000, 2.000, 2.500, 5.000]

    # corresponding tick sub division intervals
    S = [5.0,   4.0,   5.0,   5.0,
         5.0,   4.0,   5.0,   5.0,
         5.0,   4.0,   5.0,   5.0]

    span = stop - start                         # get span
    d = exp(ln10 * floor(log10(span)))          # find decade
    span /= d                                   # re-scale

    # find number of ticks below and closest to n
    i, m = 0, floor(span / T[0])                # start up
    while m > ticks:                            # next try?
        i, m = i + 1, floor(span / T[i + 1])    # try again 

    # re-scale
    mi =  d * T[i]   # main tick intervals
    si = mi / S[i]   # sub tick intervals

    # done
    return mi, si

def _getTickPositions(start, stop, ticks):

    # get intervals
    mi, si = _getTickIntervals(start, stop, ticks)

    # main ticks (round is the built-in python version)
    ns = ceil(start / mi - 0.001) * mi  # start
    ne = floor(stop / mi + 0.001) * mi  # end
    p  = round((ne - ns) / mi) + 1      # fail safe
    M  = linspace(ns, ne, p)            # main positions

    # sub ticks (round is the built-in python version)
    ns = ceil(start / si + 0.001) * si  # start
    ne = floor(stop / si - 0.001) * si  # end
    p  = round((ne - ns) / si) + 1      # fail safe
    S  = linspace(ns, ne, p)            # sub positions

    # done
    return M, S

class Document():

    def __init__(self, pathname = None):
        if pathname is not None:
            self._DOC = self.opendocument(pathname)
        return

    def opendocument(self, pathname):
        self._DOC = PdfPages(pathname)
        return self._DOC

    def exportfigure(self, name):
        args = selectfigure(name)
        self._DOC.savefig(args[0])
        return

    def closedocument(self):
        self._DOC.close()
        return

#####################################################################
#                                                                   #
#                      SETUP ENGINEERING UNITS                      #
#                                                                   #
#####################################################################

def engineering(data):
    return {
         0: (1E+00,  ""),
        -1: (1E+03, "m"),
        -2: (1E+06, "µ"),
        -3: (1E+09, "n"),
        -4: (1E+12, "p"),
        +1: (1E-03, "K"),
        +2: (1E-06, "M"),
        +3: (1E-09, "G"),
        +4: (1E-12, "T"),
    }[int(floor(log10(max(data)-min(data))/3))]

#####################################################################
#                                                                   #
#                       create pdf document                         #
#                                                                   #
#####################################################################

D = Document()
D.opendocument("./display.pdf")

#####################################################################
#                                                                   #
#                             COLORS                                #
#                                                                   #
#####################################################################

data_color = CSS4_COLORS["grey"]
fit_color  = CSS4_COLORS["black"]

#####################################################################
#                                                                   #
#                         IMPORT RESULTS                            #
#                                                                   #
#####################################################################

from numpy import load

DATA = load(f"./decay.npz")
DATA_T = DATA["time"]
DATA_A = DATA["amplitude"]
DATA_P = DATA["period"]

T, X = DATA_T, DATA_A

# get engineering units
scaler_t, suffix_t = engineering(T)
scaler_x, suffix_x = engineering(X)

# re-scale data to engineering units
T, X = T*scaler_t, X*scaler_x

#####################################################################
#                                                                   #
#                         FITTING FUNCTION                          #
#                                                                   #
#####################################################################

# in-phase fitting function
def ff(t, w):
    from numpy import exp
    x = t/w
    y = amplitude*exp(-x)
    return y

amplitude = X[0]

# find half height position
from numpy import searchsorted
i = searchsorted(-X, -X[0]/2)
th, xh = T[i], X[i]
# compute pre-fit parameters
P1 = [th/log(2)]
# fit exponential decay time
# P1, C = fit(ff, T, X, p0 = P1)
# extract parameters
w1 = P1

amplitude = X[i]

# find second half maximum
j = searchsorted(-X, -X[i]/2)
th, xh = T[j], X[j]
P2 = [th/log(2)]
# fit exponential decay time
# P2, C = fit(ff, T, X, p0 = P2)
# extract parameters
w2 = P2

#######################

fg, ax = selectfigure(f"AMPLITUDE")

headerText(f"""
    RESULTS: AMPLITUDE
    tau1 = {w1}
    tau2 = {w2}
    """, fg)

# fix T labels and ticks
ts, te = min(T), max(T)
dt = 0.1*(te-ts)
ax.set_ylim(ts-dt, te+dt)
MT, ST = _getTickPositions(ts-dt, te+dt, 7)
ax.set_xticks(MT)
ax.set_xticks(ST, minor = True)

# fix X labels and ticks
xs, xe = min(X), max(X)
dx =  0.1*(xe-xs)
ax.set_ylim(xs-dx, xe+dx)
MX, SX = _getTickPositions(xs-dx, xe+dx, 7)
ax.set_yticks(MX)
ax.set_yticks(SX, minor = True)

# fix grid style
ax.tick_params(axis = "both", which = "both", direction = "in")
ax.grid("on", which = "minor", linewidth = 0.5)
ax.grid("on", which = "major", linewidth = 1.0)

# set axes labels
ax.set_xlabel(f"Time / {suffix_t}S")
ax.set_ylabel(f"Amplitude / {suffix_x}V")

# plot data and fit
ax.plot(T, X, '.', color = data_color)

amplitude = X[0]
ax.plot(T, ff(T, *P1), "-k", linewidth = 1.0)

amplitude = X[i]
ax.plot(T, ff(T, *P2), "-k", linewidth = 1.0)

D.exportfigure(f"AMPLITUDE")

#####################################################################

# load data
T = DATA_T
X = 1.0/DATA_P

# get engineering units
scaler_t, suffix_t = engineering(T)
scaler_x, suffix_x = engineering(X)

# re-scale data to engineering units
T, X = T*scaler_t, X*scaler_x

fg, ax = selectfigure(f"FREQUENCY")

headerText(f"""
    RESULTS: FREQUENCY
    """, fg)

# fix T labels and ticks
ts, te = min(T), max(T)
dt = 0.1*(te-ts)
ax.set_ylim(ts-dt, te+dt)
MT, ST = _getTickPositions(ts-dt, te+dt, 7)
ax.set_xticks(MT)
ax.set_xticks(ST, minor = True)

# fix X labels and ticks
xs, xe = min(X), max(X)
xs, xe = 95.8, 96.0
xs, xe = 95.8-0.5, 96.0+0.5
dx =  0.1*(xe-xs)
ax.set_ylim(xs-dx, xe+dx)
MX, SX = _getTickPositions(xs-dx, xe+dx, 7)
ax.set_yticks(MX)
ax.set_yticks(SX, minor = True)

# fix grid style
ax.tick_params(axis = "both", which = "both", direction = "in")
ax.grid("on", which = "minor", linewidth = 0.5)
ax.grid("on", which = "major", linewidth = 1.0)

# set axes labels
ax.set_xlabel(f"Time / {suffix_t}S")
ax.set_ylabel(f"Period / {suffix_x}S")

# plot data and fit
ax.plot(T, X, '.', color = data_color)

from scipy.signal import savgol_filter
Y = savgol_filter(X, 50, 1)
Y = savgol_filter(Y, 50, 1)
ax.plot(T, Y, '-', color = fit_color)

D.exportfigure(f"FREQUENCY")

#####################################################################
# finalize document
#####################################################################

D.closedocument()
