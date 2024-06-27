#!/usr/bin/python3
# file: processrawfits.py
# author: Roch Schanen
# date: 2024 06 21
# content:
# repository:

# from package: "https://scipy.org/"
from scipy.optimize import curve_fit as fit

# from local modules
from figlib     import Document, stdFig, headerText
from ielib      import displayFileSize, importCSV

#####################################################################
#                                                            IMPORT #
#####################################################################

fp = {
    "large"         : "./Recording 3.csv",
    "medium"        : "./Recording 2.csv",
    }["medium"]

displayFileSize(fp)
DATA = importCSV(fp)

#####################################################################
#                                                            HEADER #
#####################################################################

def makeHeader(fg):
    headerText(fg, f"""
    FRAME {FRAME_NUMBER}:
    {ww} points from {wp} to {wp+ww-1}
    time offset      = {DATA[wp, 0]:08.4f}S
    time frame width =
    first center     = 
    period           = 
    frequency        = 
    amplitude        = 
    """)

#####################################################################
#                                                      QUICK FORMAT #
#####################################################################

def ENG(value, units = ""):
    from figlib import engineerFormat
    scaler, prefix = engineerFormat(value)
    return f"{value*scaler:.3f}{prefix}{units}"

#####################################################################
#                                                        PARAMETERS #
#####################################################################

TIME_INTERVAL = DATA[1, 0] - DATA[0, 0]
APPROX_FREQ   = 96  # Hz
FRAME0_LENGTH = int(1.25/APPROX_FREQ/TIME_INTERVAL)

print(f"""
interval = {ENG(TIME_INTERVAL, "S")} 
approximate frequency = {ENG(APPROX_FREQ, "Hz")} 
frame 0 length = {ENG(FRAME0_LENGTH)} 
""")

#####################################################################
#                                                           FRAME 0 #
#####################################################################

wp, ww, FRAME_NUMBER = 0, FRAME0_LENGTH, 0
T, X = DATA[wp:wp+ww, 0] - DATA[wp, 0], DATA[wp:wp+ww, 1]
D = Document()
D.opendocument("./display.pdf")
fg, ax = stdFig(f"FRAME{FRAME_NUMBER}", T, "S", "Time", X, "V", "Signal")

# ax.plot(T, Y, '-.', color = "r")
# ax.plot(T, ff(T, *P1), "-", color = fit_color)
# ax.plot(T[J], Y[J], 'ko', markerfacecolor = 'white', markersize = 8)
# marker='o', markeredgewidth=1.5, markersize=16
# ax.legend(["data", "filtered", "fitted", "down zero crossings"])

makeHeader(fg)
D.exportfigure(f"FRAME{FRAME_NUMBER}")
D.closedocument()
