# file = psdf.py
# created = 2024 06 25
# author = Roch Schanen
# content = phase sensitive digital filter
# comment =

from numpy import linspace
from numpy import sin
from numpy import pi

# define test signal
TEST_FREQ = 96.0 	# FREQUENCY 	Hertz
TEST_PHAS = 10.0 	# PHASE 		Degrees
TEST_AMP  = 0.2  	# AMPLITUDE 	Volts 
TEST_INT  = 20E-6 	# INTERVALS		Seconds
TEST_LEN  = 1.0		# LENGTH		Seconds
TEST_PTS  = int(TEST_LEN / TEST_INT) + 1

# compute test signal
T = linspace(0.0, TEST_LEN, TEST_PTS)

print(T)
