#!/usr/bin/python3
# file: display.py
# author: Roch Schanen
# created: 2024 06 21
# content:
# repository:

# import results
from numpy import load
from ielib import fBin
from sys import getsizeof

files = {
    "000mK" : "./.outputs/Recording 3.fit.npz",
    "600mK" : "./.outputs/Ringdown_1325_600mK.fit.npz",
    "750mK" : "./.outputs/Ringdown_1325mV_750mK.fit.npz",
    }

DATA = {k:load(files[k]) for k in files.keys()}

PLOTS = [
    ("600mK",  0.0, ["r"], {}),
    ("000mK", 11.0, ["b"], {}),
    ("750mK",  0.0, ["k"], {}),
]

###################
# plot amplitudes #
###################

from numpy import log
from figlib import Document
from figlib import stdFigure

D = Document()
D.opendocument("./.outputs/display.pdf")

fa = stdFigure(f"FIG_A", "Time", "S", r"Amplitude", "V")
fl = stdFigure(f"FIG_L", "Time", "S", r"ln(a/a$_0$)", "")
ff = stdFigure(f"FIG_F", "Time", "S", r"Frequency", "Hz")

for k, t, args, kwargs in PLOTS:

    T = DATA[k]["DATA_TIME"]
    A = DATA[k]["DATA_AMPLITUDE"]
    P = DATA[k]["DATA_PERIOD"]

    from scipy.signal import savgol_filter
    F = savgol_filter(1.0/P, 250, 2)

    fa.plot(T+t, A, *args, **kwargs)
    fl.plot(T+t, log(A/A[0]), *args, **kwargs)
    # ff.plot(T, 1.0/P, *args, **kwargs)
    ff.plot(T+t, F, *args, **kwargs)

D.exportfigure(f"FIG_A")
D.exportfigure(f"FIG_L")
D.exportfigure(f"FIG_F")
D.closedocument()
