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

# compute 3-decades 
def vEng(value):
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
    return C, S

def fEng(value, format = ".3f"):
    C, S = vEng(value)
    return f"{value*C:{format}}{S}"

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
    # find color keys at
    # matplotlib.org/3.1.0/gallery/color/named_colors.html
    data_color = CSS4_COLORS["grey"]

    fg, ax = selectfigure(name)

    fg.scx, fg.suffix_x = vEng(max(X)-min(X))
    fg.scy, fg.suffix_y = vEng(max(Y)-min(Y))

    # re-scale data
    X, Y = X*fg.scx, Y*fg.scy

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
    ax.grid("on", which = "minor", linewidth = 0.5, linestyle = "--")
    ax.grid("on", which = "major", linewidth = 1.0)

    # set axes labels
    XU = f" / {fg.suffix_x}{xu}" if xu else ""
    ax.set_xlabel(f"{xl}{XU}")
    YU = f" / {fg.suffix_y}{yu}" if yu else ""
    ax.set_ylabel(f"{yl}{YU}")

    # plot data
    ax.plot(X, Y, '-', color = data_color)

    return fg, ax

#####################################################################
#                                                 STANDARD FIGURE 2 #
#####################################################################

def scaleParameters_0(*args):
    vmins, vmaxs = [], []
    for V in args:
        vmins.append(min(V))
        vmaxs.append(max(V))
    vspan = max(vmaxs) - min(vmins)
    va, vf = vEng(vspan)
    vs, ve = min(vmins)*va, max(vmaxs)*va
    vd = 0.1*(ve - vs)
    vs, ve = vs-vd, ve+vd
    M, S = getTickPositions(vs, ve, 7)
    return va, vf, vs, ve, M, S

def setXTicks(ax, *args):
    scale, suffix, s, e, M, S = scaleParameters_0(*args)
    ax.set_xlim(s, e)
    ax.set_xticks(M)
    ax.set_xticks(S, minor = True)
    return scale, suffix

def setYTicks(ax, *args):
    scale, suffix, s, e, M, S = scaleParameters_0(*args)
    ax.set_ylim(s, e)
    ax.set_yticks(M)
    ax.set_yticks(S, minor = True)
    return scale, suffix

def stdfig(name, xl, xu, X, yl, yu, *Yargs, **kwargs):
    fg, ax = selectfigure(name)
    # setup X scale
    fg.scx, fg.suffix_x = setXTicks(ax, X)
    XU = f" / {fg.suffix_x}{xu}" if xu else ""
    ax.set_xlabel(f"{xl}{XU}")
    # setup Y scale
    fg.scy, fg.suffix_y = setYTicks(ax, *Yargs)
    YU = f" / {fg.suffix_y}{yu}" if yu else ""
    ax.set_ylabel(f"{yl}{YU}")
    # setup grid
    ax.tick_params(axis = "both", which = "both", direction = "in")
    ax.grid("on", which = "minor", linewidth = 0.5, linestyle = "--")
    ax.grid("on", which = "major", linewidth = 1.0)
    # matplotlib.org/3.1.0/gallery/color/named_colors.html
    from matplotlib.colors import CSS4_COLORS
    fg.colors = CSS4_COLORS
    return fg, ax

def stdplot(name, X, Y, *args, **kwargs):
    fg, ax = selectfigure(name)
    ax.plot(X*fg.scx, Y*fg.scy, *args, **kwargs)
    return

#####################################################################
#                                                 STANDARD FIGURE 3 #
#####################################################################

class stdFigure():

    def axisScale(self, s, e):
        # get engineer values
        scale, suffix = vEng(e-s)
        # rescale
        s, e = s*scale, e*scale
        # expand
        d = 0.1*(e-s)
        s, e = s-d, e+d
        # done
        return scale, suffix, s, e

    def setXTicks(self, s, e, tiks = 7):
        self.scx, self.sfx, s, e = self.axisScale(s, e)
        M, S = getTickPositions(s, e, tiks)
        self.ax.set_xlim(s, e)
        self.ax.set_xticks(M)
        self.ax.set_xticks(S, minor = True)
        return

    def setYTicks(self, s, e, tiks = 7):
        self.scy, self.sfy, s, e = self.axisScale(s, e)
        M, S = getTickPositions(s, e, tiks)
        self.ax.set_ylim(s, e)
        self.ax.set_yticks(M)
        self.ax.set_yticks(S, minor = True)
        return

    def __init__(self, name,
        xl = "X", xu = "", # x labels and units
        yl = "Y", yu = "", # y labels and units
        ):
        # figure name
        self.fn = name
        # instanciate matplotlib figure
        self.fg, self.ax = selectfigure(name)
        # set default style
        self.ax.tick_params(axis = "both",  which = "both",  direction = "in")
        self.ax.grid("on", which = "minor", linewidth = 0.5, linestyle = "--")
        self.ax.grid("on", which = "major", linewidth = 1.0)
        # record other parameters
        self.labels = xl, yl
        self.units  = xu, yu
        self.data = []
        # copy reference to matplotlib colors
        from matplotlib.colors import CSS4_COLORS
        self.color = CSS4_COLORS
        # link back self to matplotlib figure
        self.fg.stdFigureReference = self
        # try first plot
        self.rescale()
        # done
        return

    def plot(self, X, Y, *args, **kwargs):
        # add plot with styles
        line, = self.ax.plot(
            X*self.scx,
            Y*self.scy,
            *args, **kwargs)
        # record original set
        self.data.append((X, Y, line))
        # rescale axes
        self.rescale()
        # done 
        return

    def rescale(self):
        # defaults range is -1, +1
        xs, xe, ys, ye = -1.0, +1.0, -1.0, +1.0
        if self.data:
            # recompute mins and maxs for all axes
            xmins, xmaxs, ymins, ymaxs = [], [], [], []
            for X, Y, line in self.data:
                xmins.append(min(X)); xmaxs.append(max(X))
                ymins.append(min(Y)); ymaxs.append(max(Y))
            xs, xe = min(xmins), max(xmaxs)
            ys, ye = min(ymins), max(ymaxs)
        # get axis labels and units
        xl, yl = self.labels
        xu, yu = self.units
        # setup X scale
        self.setXTicks(xs, xe)
        XU = f" / {self.sfx}{xu}" if xu else ""
        self.ax.set_xlabel(f"{xl}{XU}")
        # setup Y scale
        self.setYTicks(ys, ye)
        YU = f" / {self.sfy}{yu}" if yu else ""
        self.ax.set_ylabel(f"{yl}{YU}")
        # re-scale all
        if self.data:
            for X, Y, line in self.data:
                line.set_xdata(X*self.scx)
                line.set_ydata(Y*self.scy)
        # done
        return

    def header(self, text):
        w, h = array([1, 1 / 1.4143])*0.7
        x, y = (1-w)/2, (1-h)/2
        tx = self.fg.text(x, 3*y/2+h, text)
        tx.set_fontfamily('monospace')
        tx.set_horizontalalignment('left')
        tx.set_verticalalignment('center')
        tx.set_fontsize("small")
        return tx
