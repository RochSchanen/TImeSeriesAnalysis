#!/usr/bin/python3
# file: figlib.py
# author: Roch Schanen
# date: 2024 06 21
# content:
# repository:

#####################################################################
#                                                          PACKAGES #
#####################################################################

# matplotlib    from "https://matplotlib.org/"
# numpy         from "https://numpy.org/"

from numpy import array
from numpy import linspace
from numpy import ceil, floor
from numpy import log10, exp

#####################################################################
#                                                            FIGURE #
#####################################################################

# add other page formats

def selectfigure(name):
    from matplotlib.pyplot import figure
    from matplotlib.pyplot import fignum_exists
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

#####################################################################
#                                                            HEADER #
#####################################################################

def headerText(fg, text):
    w, h = array([1, 1 / 1.4143])*0.7
    x, y = (1-w)/2, (1-h)/2
    tx = fg.text(x, 3*y/2+h, text)
    tx.set_fontfamily('monospace')
    tx.set_horizontalalignment('left')
    tx.set_verticalalignment('center')
    tx.set_fontsize("small")
    return tx

#####################################################################
#                                                            FOOTER #
#####################################################################

def footerText(fg, text):
    w, h = array([1, 1 / 1.4143])*0.7
    x, y = (1-w)/2, (1-h)/2
    tx = fg.text(x, y/2, text)
    tx.set_fontfamily('monospace')
    tx.set_horizontalalignment('left')
    tx.set_verticalalignment('center')
    tx.set_fontsize("small")
    return tx

#####################################################################
#                                                             TICKS #
#####################################################################

def getTickIntervals(start, stop, ticks):

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

def getTickPositions(start, stop, ticks):

    # get intervals
    mi, si = getTickIntervals(start, stop, ticks)

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

#####################################################################
#                                                    ENGINEER UNITS #
#####################################################################

def engineerFormat(value):
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
    }[int(floor(log10(value)/3))]

#####################################################################
#                                                          DOCUMENT #
#####################################################################

class Document():

    def __init__(self, pathname = None):
        if pathname is not None:
            self._DOC = self.opendocument(pathname)
        return

    def opendocument(self, pathname):
        from matplotlib.backends.backend_pdf import PdfPages
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
#                                                   STANDARD FIGURE #
#####################################################################

def stdFig(name, X, xu, xl, Y, yu, yl):

    from matplotlib.colors import CSS4_COLORS
    data_color = CSS4_COLORS["grey"]

    fg, ax = selectfigure(name)

    scaler_x, suffix_x = engineerFormat(max(X)-min(X))
    scaler_y, suffix_y = engineerFormat(max(Y)-min(Y))

    # re-scale data
    X, Y = X*scaler_x, Y*scaler_y

    # fix X labels and ticks
    s, e = min(X), max(X)
    d =  0.1*(e-s)
    ax.set_xlim(s-d, e+d)
    M, S = getTickPositions(s-d, e+d, 7)
    ax.set_xticks(M)
    ax.set_xticks(S, minor = True)

    # fix Y labels and ticks
    s, e = min(Y), max(Y)
    d =  0.1*(e-s)
    ax.set_ylim(s-d, e+d)
    M, S = getTickPositions(s-d, e+d, 7)
    ax.set_yticks(M)
    ax.set_yticks(S, minor = True)

    # fix grid style
    ax.tick_params(axis = "both", which = "both", direction = "in")
    ax.grid("on", which = "minor", linewidth = 0.5)
    ax.grid("on", which = "major", linewidth = 1.0)

    # set axes labels
    ax.set_xlabel(f"{xl} / {suffix_x}{xu}")
    ax.set_ylabel(f"{yl} / {suffix_y}{yu}")

    # plot data
    ax.plot(X, Y, '.', color = data_color)

    return fg, ax
