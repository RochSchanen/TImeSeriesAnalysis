# file = psdf.py
# created = 2024 06 25
# author = Roch Schanen
# content = phase sensitive digital filter
# comment =

#####################################################################
#                                                                   #
#              FIGURE, HEADER, FOOTER, TICKS, DOCUMENT              #
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

data_color = CSS4_COLORS["grey"]
fit_color  = CSS4_COLORS["black"]

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

def E(value, prec = 3):
    value = abs(value)
    C, S = (1E+00, "") if float(value) == 0.0 else {
         0: (1E+00,  ""),
        -1: (1E+03, "m"),
        -2: (1E+06, "µ"),
        -3: (1E+09, "n"),
        -4: (1E+12, "p"),
        +1: (1E-03, "K"),
        +2: (1E-06, "M"),
        +3: (1E-09, "G"),
        +4: (1E-12, "T"),
    }[int(floor(log10(value)/3))]
    fs = f".{prec}f" # format string
    return f"{value*C:{fs}}{S}"

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

from numpy import linspace, empty, zeros
from numpy import sin, cos, sqrt
from numpy import pi

# define test signal
TEST_FREQ = 96.0        # FREQUENCY             Hertz
TEST_PHAS = 0.0         # PHASE                 Degrees
TEST_AMPL = 1.0         # AMPLITUDE             Volts 
TEST_SAMP = 100E-6      # SAMPLING INTERVAL     Seconds
TEST_LEN  = 350E-3      # SIGNAL LENGTH         Seconds

TEST_PTS  = int(TEST_LEN / TEST_SAMP) + 1

print(f"Frequency = {E(TEST_FREQ)}Hz")
print(f"Phase = {E(TEST_PHAS)}°")
print(f"Amplitude = {E(TEST_AMPL)}V")
print(f"Sampling intervals = {E(TEST_SAMP)}S")
print(f"Number of points{E(TEST_PTS)}")

# time vector
T = linspace(0.0, TEST_LEN, TEST_PTS)
# phase vector
P = 2.0*pi*TEST_FREQ*T
# test signal (includes a phase shift)
V = TEST_AMPL*sin(P-pi/180.0*TEST_PHAS)

########################################
# PHASE SENSITIVE DETECTION PARAMETERS #
########################################

# time constant (same as the lockin time time constant)
PSDF_TIMC = 0.050 # Seconds

# normally we should have TEST_SAMP << PSDF_TIMC and
# PSDF_ALPH is approximately TEST_SAMP / PSDF_TIMC
PSDF_ALPH = TEST_SAMP / (TEST_SAMP + PSDF_TIMC);

# the slope depends on the number of low-pass filters
PSDF_NUM = 3

# compute reference vectors
SIN, COS = sin(P), cos(P)

# reserve memory for the low pass filter outputs
VX = zeros((TEST_PTS, PSDF_NUM))
VY = zeros((TEST_PTS, PSDF_NUM))

# sweep through buffered signal: reference and signal vectors
for i, (s, c, v) in enumerate(zip(SIN, COS, V)):
    # compute and record first stage detection output
    VX[i, 0] = VX[i-1, 0] + (s*v - VX[i-1, 0])*PSDF_ALPH
    VY[i, 0] = VY[i-1, 0] + (c*v - VY[i-1, 0])*PSDF_ALPH
    # compute and record upper low pass filter ouputs
    for j in range(1, PSDF_NUM):
        VX[i, j] = VX[i-1, j] + (VX[i, j-1] - VX[i-1, j])*PSDF_ALPH
        VY[i, j] = VY[i-1, j] + (VY[i, j-1] - VY[i-1, j])*PSDF_ALPH

# copy last filter data sets
X, Y = 2*VX[:, -1], 2*VY[:, -1]

#####################################################################

scaler_t, suffix_t = engineering(T)
scaler_v, suffix_v = engineering(V)

T, V = T*scaler_t, V*scaler_v
X, Y = X*scaler_v, Y*scaler_v

p, w, h = 1, 1, 1

fg, ax = selectfigure("FRAME0")

# headerText(f"""
#     FRAME 0: determination of the pre-fitting values
#     {ww} points from {wp} to {wp+ww-1}
#     time offset      = {T[wp]:08.4f}S
#     time frame width = {T[-1]:08.4f}{suffix_t}S
#     first center     = {p:.3f}{suffix_t}S
#     period           = {w:.4f}{suffix_t}S
#     frequency        = {scaler_t/w:.4f}Hz
#     amplitude        = {h:.1f}{suffix_x}V
#     """, fg)

# fix T labels and ticks
ts, te = min(T), max(T)
dt = 0.1*(te-ts)
ax.set_ylim(ts-dt, te+dt)
MT, ST = _getTickPositions(ts-dt, te+dt, 7)
ax.set_xticks(MT)
ax.set_xticks(ST, minor = True)

# fix X labels and ticks
vs, ve = min(V), max(V)
dv =  0.1*(ve-vs)
ax.set_ylim(vs-dv, ve+dv)
MV, SV = _getTickPositions(vs-dv, ve+dv, 7)
ax.set_yticks(MV)
ax.set_yticks(SV, minor = True)

# fix grid style
ax.tick_params(axis = "both", which = "both", direction = "in")
ax.grid("on", which = "minor", linewidth = 0.5)
ax.grid("on", which = "major", linewidth = 1.0)

# set axes labels
ax.set_xlabel(f"Time / {suffix_t}S")
ax.set_ylabel(f"Signal / {suffix_v}V")

# plot data
ax.plot(T, V, '-', color = data_color)
ax.plot(T, X, '-', color = "r")
ax.plot(T, Y, '-', color = "b")
# ax.plot(T, Y, '-.', color = "r")
# ax.plot(T, ff(T, *P), "-", color = fit_color)
# ax.plot(T[J], Y[J],
#     'ko', 
#     markerfacecolor = 'white',
#     markersize = 8)
ax.legend(["data"])

D.exportfigure("FRAME0")

D.closedocument()
