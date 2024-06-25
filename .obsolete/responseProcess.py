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
#                       create pdf document                         #
#                                                                   #
#####################################################################

D = Document()
D.opendocument("./display.pdf")

#####################################################################
#                                                                   #
#                            IMPORT DATA                            #
#                                                                   #
#####################################################################

# just for info
def displayfilesize(fp):
    from os import stat
    fs = stat(fp).st_size
    print(f'File Size is {fs:.0f} Bytes.')
    print(f'File Size is {fs/1024:.0f} KB.')
    print(f'File Size is {fs/1024/1024:.0f} MB.')
    return fs

# select file (debug)
fp = {
    "verylarge"     : "./Recording 3.csv",
    "medium"        : "./partial.csv",
    "small"         : "./partialsmall.csv",
    "large"         : "./partiallarge.csv",
    }["verylarge"]

# display file size (debug)
displayfilesize(fp)

# from numpy import genfromtxt
from numpy import loadtxt

# data = genfromtxt(fp, delimiter=',', skip_header = 4)
# # the above takes 36s for 15055399 data points

data = loadtxt(fp, delimiter=',', skiprows = 4)
# the above takes 12.7s for 15055399 data points

print(f"number of points = {data[:,0].size}.")

from sys import getsizeof
print(f"memory usage = {getsizeof(data)} Bytes.")

# split data array, select window zero and shift time (zero shift)
wp, ww = 0, 700
T, X = data[wp:wp+ww, 0] - data[wp, 0], data[wp:wp+ww, 1]

# get time interval
TIME_INTERVAL = T[1]-T[0]

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

# get engineering units
scaler_t, suffix_t = engineering(T)
scaler_x, suffix_x = engineering(X)

# re-scale data
T, X = T*scaler_t, X*scaler_x

#####################################################################
#                                                                   #
#                            PRE-FIT DATA                           #
#                                                                   #
#####################################################################

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

#####################################################################
#                                                                   #
#                             COLORS                                #
#                                                                   #
#####################################################################

data_color = CSS4_COLORS["grey"]
fit_color  = CSS4_COLORS["black"]

#####################################################################
#                                                                   #
#                             FRAME 0                               #
#                                                                   #
#####################################################################

FRAME_NUMBER = 0

fg, ax = selectfigure("FRAME0")

headerText(f"""
    FRAME 0: determination of the pre-fitting values
    {ww} points from {wp} to {wp+ww-1}
    time offset      = {data[wp, 0]:08.4f}S
    time frame width = {T[-1]:08.4f}{suffix_t}S
    first center     = {p:.3f}{suffix_t}S
    period           = {w:.4f}{suffix_t}S
    frequency        = {scaler_t/w:.4f}Hz
    amplitude        = {h:.1f}{suffix_x}V
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
ax.set_ylabel(f"Signal / {suffix_x}V")

# plot data
ax.plot(T, X, '.', color = data_color)
ax.plot(T, Y, '-.', color = "r")
ax.plot(T, ff(T, *P), "-", color = fit_color)
ax.plot(T[J], Y[J], 'ko', markerfacecolor = 'white', markersize = 8)
ax.legend(["data", "filtered", "fitted", "down zero crossings"])

D.exportfigure("FRAME0")

#####################################################################
#                                                                   #
#                     DECLARE STORAGE LISTS                         #
#                                                                   #
#####################################################################

# estimate array size
pts = data[:,0].size
cycle_pts = i2-i1+1
cycles = pts / cycle_pts
FRAME_NUMBERS = int(ceil(cycles / 5))

print(f"max number of frames = {FRAME_NUMBERS}")

DATA_TIME = zeros(FRAME_NUMBERS)
DATA_PERIOD = zeros(FRAME_NUMBERS)
DATA_AMPLITUDE = zeros(FRAME_NUMBERS)

#####################################################################
#                                                                   #
#                             FRAME 1                               #
#                                                                   #
#####################################################################

# first frame
wp, ww = 0, 4400 # approximately 8 cycles

# increment frame number
FRAME_NUMBER += 1

# select window and shift time origin
T, X = data[wp:wp+ww, 0] - data[wp, 0], data[wp:wp+ww, 1]

# get engineering units
scaler_t, suffix_t = engineering(T)
scaler_x, suffix_x = engineering(X)

# re-scale data to engineering units
T, X = T*scaler_t, X*scaler_x

# re-fit using previous results
P, C = fit(ff, T, X, p0 = P)

# extract parameters
p, w, h = P

# FIGURE "FRAME1"
fg, ax = selectfigure(f"FRAME{FRAME_NUMBER}")

headerText(f"""
    FRAME {FRAME_NUMBER}:
    {ww} points from {wp} to {wp+ww-1}
    time offset      = {data[wp, 0]:08.4f}S
    time frame width = {T[-1]:08.4f}{suffix_t}S
    first center     = {p:.3f}{suffix_t}S
    period           = {w:.4f}{suffix_t}S
    frequency        = {scaler_t/w:.4f}Hz
    amplitude        = {h:.1f}{suffix_x}V
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
ax.set_ylabel(f"Signal / {suffix_x}V")

# plot data and fit
ax.plot(T, X, '.', color = data_color)
ax.plot(T, ff(T, *P), "-k", linewidth = 1.0)
ax.legend(["data", "fitted"])

D.exportfigure(f"FRAME{FRAME_NUMBER}")

# COLLECT DATA
DATA_TIME[FRAME_NUMBER-1] = mean(T)/scaler_t
DATA_PERIOD[FRAME_NUMBER-1] = w/scaler_t
DATA_AMPLITUDE[FRAME_NUMBER-1] = h/scaler_x

#####################################################################
#                                                                   #
#                         START LOOP                                #
#                                                                   #
#####################################################################

# manual breaks

# FRAME_NUMBERS = 50
# displaypages = [10, 20, 30, 40, FRAME_NUMBERS-1]

# FRAME_NUMBERS = 1000
# displaypages = [100, 500, FRAME_NUMBERS-1]

# FRAME_NUMBERS = 3000
# displaypages = [100, 500, 1000, 2000, FRAME_NUMBERS-1]

# FRAME_NUMBERS = 5000
# displaypages = [100, 1000, 3500, FRAME_NUMBERS-1]

FRAME_NUMBERS = 4800
displaypages = [2, 10, 100, 500, 1000, 2000, 3000, 4000, FRAME_NUMBERS-1]

# FRAME_NUMBERS = 4800
# displaypages = [FRAME_NUMBERS-1]

while True:

    #####################################################################
    #                                                                   #
    #                         FRAME SHIFT                               #
    #                                                                   #
    #####################################################################

    # shift by five periods
    wp += int(5*w/scaler_t/TIME_INTERVAL)
    # a motivation to shift by an integer multiple of cycles is that
    # it leaves the initial parameters unchanged. We can re-use the
    # previous results as the initial parameters.
    if wp+ww>pts : break # safe break (not enougth data available)

    # increment frame number
    FRAME_NUMBER += 1

    #####################################################################
    #                                                                   #
    #                     LOAD NEXT FRAME DATA                          #
    #                                                                   #
    #####################################################################

    # select window and shift time origin
    T, X = data[wp:wp+ww, 0] - data[wp, 0], data[wp:wp+ww, 1]

    # get engineering units
    scaler_t, suffix_t = engineering(T)
    scaler_x, suffix_x = engineering(X)

    # re-scale data to engineering units
    T, X = T*scaler_t, X*scaler_x

    #####################################################################
    #                                                                   #
    #                     COMPUTE FRAME RESULTS                         #
    #                                                                   #
    #####################################################################

    # fit using previous results
    P, C = fit(ff, T, X, p0 = P)

    # extract parameters
    p, w, h = P

    #####################################################################
    #                                                                   #
    #                        ADD TO DISPLAY                             #
    #                                                                   #
    #####################################################################

    if FRAME_NUMBER in displaypages:

        fg, ax = selectfigure(f"FRAME{FRAME_NUMBER}")

        headerText(f"""
            FRAME {FRAME_NUMBER}:
            {ww} points from {wp} to {wp+ww-1}
            time offset      = {data[wp, 0]:08.4f}S
            time frame width = {T[-1]:08.4f}{suffix_t}S
            first center     = {p:.3f}{suffix_t}S
            period           = {w:.4f}{suffix_t}S
            frequency        = {scaler_t/w:.4f}Hz
            amplitude        = {h:.1f}{suffix_x}V
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
        ax.set_ylabel(f"Signal / {suffix_x}V")

        # plot data and fit
        ax.plot(T, X, '.', color = data_color)
        ax.plot(T, ff(T, *P), "-k", linewidth = 1.0)
        ax.legend(["data", "fitted"])

        D.exportfigure(f"FRAME{FRAME_NUMBER}")

    #####################################################################
    #                                                                   #
    #                        COLLECT DATA                               #
    #                                                                   #
    #####################################################################

    # average time of the data fitted in seconds
    DATA_TIME[FRAME_NUMBER-1] = data[wp, 0] + mean(T)/scaler_t
    # period in seconds
    DATA_PERIOD[FRAME_NUMBER-1] = w/scaler_t
    # signal in volts
    DATA_AMPLITUDE[FRAME_NUMBER-1] = h/scaler_x

    #####################################################################
    #                                                                   #
    #                         END LOOP                                  #
    #                                                                   #
    #####################################################################

    if FRAME_NUMBER == FRAME_NUMBERS: break


#####################################################################
#                                                                   #
#                         EXPORT RESULTS                            #
#                                                                   #
#####################################################################

from numpy import savez, savetxt

savez(f"./decay.npz",
    time = DATA_TIME[:FRAME_NUMBER],
    amplitude = DATA_AMPLITUDE[:FRAME_NUMBER],
    period = DATA_PERIOD[:FRAME_NUMBER],
    )

# savetxt()

#####################################################################
# finalize document
#####################################################################

D.closedocument()
