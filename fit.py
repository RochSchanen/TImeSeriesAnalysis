#!/usr/bin/python3
# file: fit.py
# author: Roch Schanen
# created: 2024 07 10
# content:

BLOCK_SIZE = 1<<20 # 1MB

#############
### USAGE ###
#############

from sys import argv, exit

if len(argv)<2:
    print("""    --- usage ---
    > python3 fit.py "source.csv" "destination.dat" "frequency" "span" "resolution"
    - frequency: estimated frequency [Hz] (used to compute the pre-fit parameters).
    - span: the time window used to fit the harmonic signal [S].
    - resolution: the time by which the window is shifted between harmonic fits [S].
    The resolution is recommanded to be similar or less than the value of the span.
    """)
    exit()

###################
### SOURCE FILE ###
###################

from sys import argv, exit
from os.path import exists as validPath

if not validPath(argv[1]):
    print(f"file '{argv[1]}' not found.")
    print(f"exiting...")
    exit()

fpi = argv[1]

###################
### OUTPUT FILE ###
###################

from os.path import split as splitPath
p, n = splitPath(argv[2])

if p.strip() == "":
    p = f".outputs"

from os.path import isdir as validDir
if not validDir(p):
    print(f"directory '{p}' not found.")
    print(f"exiting...")
    exit()

if n.strip() == "":
    n = splitPath(argv[1])[1]

print(f"{p}/{n}")











# from os.path import splitext
# n, e = splitext(n)
