# file = psdf.py
# created = 2024 06 25
# author = Roch Schanen
# content = phase sensitive digital filter
# comment =

def PSD_RAW(PSD_FRQ, PSD_TMC, PSD_NUM, T, V):
    # - PSD_FRQ is frequency
    # - PSD_TMC is timer constant
    # - PSD_NUM is number of low pass filters:
    # 1 ->  -6dB (wait  5 PSD_TMC)
    # 2 -> -12dB (wait  7 PSD_TMC)
    # 3 -> -18dB (wait  9 PSD_TMC)
    # 4 -> -24dB (wait 10 PSD_TMC)
    # - T is the time vector in units of Seconds
    # - V is the signal vector in units of Volts
    
    # the time sampling intervals must be fixed
    SAMPLING = T[1]-T[0]
    # imports
    from numpy import pi, sin, cos, zeros
    # compute phase vector
    P = 2.0*pi*PSD_FRQ*T
    # compute reference vectors
    SIN, COS = sin(P), cos(P)
    # compute digital low pass filter alpha value
    PSD_ALPH = SAMPLING / (SAMPLING + PSD_TMC);
    # reserve memory for the low pass filter outputs
    VX = zeros((len(T), PSD_NUM))
    VY = zeros((len(T), PSD_NUM))
    # sweep through buffered signal: reference and signal vectors
    for i, (s, c, v) in enumerate(zip(SIN, COS, V)):
        # compute and record the first stage detection output
        VX[i, 0] = VX[i-1, 0] + (s*v - VX[i-1, 0])*PSD_ALPH
        VY[i, 0] = VY[i-1, 0] + (c*v - VY[i-1, 0])*PSD_ALPH
        # compute and record the next low pass filters state
        for j in range(1, PSD_NUM):
            VX[i, j] = VX[i-1, j] + (VX[i, j-1]-VX[i-1, j])*PSD_ALPH
            VY[i, j] = VY[i-1, j] + (VY[i, j-1]-VY[i-1, j])*PSD_ALPH
    # done, return the last low pass filters data sets
    return VX, VY

def PSD_RMS(PSD_FRQ, PSD_TMC, PSD_NUM, T, V, LPFS = None):
    # - PSD_FRQ is frequency
    # - PSD_TMC is timer constant
    # - PSD_NUM is number of low pass filters:
    # 1 ->  -6dB (wait  5 PSD_TMC)
    # 2 -> -12dB (wait  7 PSD_TMC)
    # 3 -> -18dB (wait  9 PSD_TMC)
    # 4 -> -24dB (wait 10 PSD_TMC)
    # - T is the time vector in units of Seconds
    # - X is the signal vector in units of Volts
    # - use LPFS to initialise the low pass filter levels

    # the time sampling intervals must be fixed
    SAMPLING = T[1]-T[0]
    # imports
    from numpy import pi, sin, cos, zeros, sqrt
    # compute phase vector
    P = 2.0*pi*PSD_FRQ*T
    # compute reference vectors
    SIN, COS = sin(P), cos(P)
    # compute digital low pass filter alpha value
    PSD_ALPH = SAMPLING / (SAMPLING + PSD_TMC)
    # reserve memory for the low pass filter outputs
    VX = zeros((len(T), PSD_NUM))
    VY = zeros((len(T), PSD_NUM))
    # setup initial low pass filter levels (defaults to zeros)
    if LPFS is not None: VX[0, :], VY[0, :] = LPFS
    # sweep through references and the signal vector V
    for i, (s, c, v) in enumerate(zip(SIN, COS, V)):
        # use start up value on the first iteration
        k = 0 if i == 0 else i-1
        # compute and record the first stage detection output
        VX[i, 0] = VX[k, 0] + (s*v - VX[k, 0])*PSD_ALPH
        VY[i, 0] = VY[k, 0] + (c*v - VY[k, 0])*PSD_ALPH
        # compute and record the next low pass filters state
        for j in range(1, PSD_NUM):
            VX[i, j] = VX[k, j] + (VX[i, j-1]-VX[k, j])*PSD_ALPH
            VY[i, j] = VY[k, j] + (VY[i, j-1]-VY[k, j])*PSD_ALPH
    # record low pass filters final levels
    LPFS = VX[-1, :], VY[-1, :]
    # keep last low pass filter data set
    VXN, VYN = sqrt(2)*VX[:,-1], sqrt(2)*VY[:,-1]
    # done: return the last low pass filters data sets
    # with the last low pass filter levels
    return VXN, VYN, LPFS

#####################################################################

if __name__ == "__main__":

    # required package: "https://numpy.org/"
    from figlib import Document
    from figlib import stdFig
    from figlib import fEng
    
    from numpy import linspace
    from numpy import square
    from numpy import sqrt
    from numpy import sin
    from numpy import pi

    #################################################################

    D = Document()
    D.opendocument("./display.pdf")

    #################################################################

    ##########################
    # TEST SIGNAL PARAMETERS #
    ##########################

    # define test signal
    TEST_FREQ = 96.0            # FREQUENCY             Hertz
    TEST_PHAS = 45.0            # PHASE                 Degrees
    TEST_AMPL = 0.1*sqrt(2)     # AMPLITUDE (0.1Vrms)   Volts 
    SAMPLING = 100E-6          # SAMPLING INTERVAL     Seconds
    TEST_LEN  = 350E-3          # SIGNAL LENGTH         Seconds
    TEST_PTS  = int(TEST_LEN / SAMPLING) + 1
    # export values (debug)
    print(f"Frequency = {fEng(TEST_FREQ)}Hz")
    print(f"Phase = {fEng(TEST_PHAS)}°")
    print(f"Amplitude = {fEng(TEST_AMPL)}V")
    print(f"Sampling intervals = {fEng(SAMPLING)}S")
    print(f"Number of points{fEng(TEST_PTS)}")
    # time vector
    T = linspace(0.0, TEST_LEN, TEST_PTS)
    # test signal (includes a phase shift)
    V = TEST_AMPL*sin(2.0*pi*TEST_FREQ*T - pi/180.0*TEST_PHAS)

    #################################################################

    ########################
    # APPLY PSD_RAW FILTER #
    ########################
     
    VX, VY = PSD_RAW(TEST_FREQ, 0.050, 2, T, V)

    ################
    # PLOT RESULTS #
    ################

    # plot data
    fg, ax = stdFig("FIG0", T, "S", "Time", V, "V", "Signal")
    # add plots to figure
    X, Y = sqrt(2)*VX[:, 0], sqrt(2)*VY[:, 0]
    A = sqrt(square(X)+square(Y))
    ax.plot(T*fg.scx, X*fg.scy, '-.', color = "b")
    ax.plot(T*fg.scx, Y*fg.scy, '-.', color = "b")
    ax.plot(T*fg.scx, A*fg.scy, '-', color = "r")
    # add plots to figure
    X, Y = sqrt(2)*VX[:, 1], sqrt(2)*VY[:, 1]
    A = sqrt(square(X)+square(Y))
    ax.plot(T*fg.scx, X*fg.scy, '-.', color = "b")
    ax.plot(T*fg.scx, Y*fg.scy, '-.', color = "b")
    ax.plot(T*fg.scx, A*fg.scy, '-', color = "r")
    # setup legends
    ax.legend(["data", "VX1", "VY1", "VA1", "VX2", "VY2", "VA2"])
    # export figure
    D.exportfigure("FIG0")

    #################################################################

    ########################
    # APPLY PSD_RMS FILTER #
    ########################

    X, Y, LPFS = PSD_RMS(TEST_FREQ, 0.050, 2, T, V)
    A = sqrt(square(X)+square(Y))
    
    ################
    # PLOT RESULTS #
    ################

    # plot data
    fg, ax = stdFig("FIG1", T, "S", "Time", V, "V", "Signal")
    # add plot to figure
    ax.plot(T*fg.scx, X*fg.scy, '-.', color = "b")
    ax.plot(T*fg.scx, Y*fg.scy, '-.', color = "b")
    ax.plot(T*fg.scx, A*fg.scy, '-', color = "r")
    # setup legends
    ax.legend(["data", "VX", "VY", "V"])
    # export figure
    D.exportfigure("FIG1")

    #################################################################

    # when files are large and memory scarce,
    # the digital filtering can be done in parts.
    # The LPFS values allows to link the filtered sections.

    #################################
    # APPLY PSD_RMS FILTER IN PARTS #
    #################################

    T1, V1 = T[0:1000], V[0:1000]
    X1, Y1, LPFS = PSD_RMS(TEST_FREQ, 0.050, 2, T1, V1)
    A1 = sqrt(square(X1)+square(Y1))

    T2, V2 = T[1001:2000], V[1001:2000]
    X2, Y2, LPFS = PSD_RMS(TEST_FREQ, 0.050, 2, T2, V2, LPFS)
    A2 = sqrt(square(X2)+square(Y2))

    T3, V3 = T[2001:], V[2001:]
    X3, Y3, LPFS = PSD_RMS(TEST_FREQ, 0.050, 2, T3, V3, LPFS)
    A3 = sqrt(square(X3)+square(Y3))
    
    ################
    # PLOT RESULTS #
    ################

    # plot data
    fg, ax = stdFig("FIG2", T, "S", "Time", V, "V", "Signal")
    # add plot to figure
    ax.plot(T1*fg.scx, X1*fg.scy, '-.', color = "b")
    ax.plot(T1*fg.scx, Y1*fg.scy, '-.', color = "b")
    ax.plot(T1*fg.scx, A1*fg.scy, '-',  color = "r")
    # add plot to figure
    ax.plot(T2*fg.scx, X2*fg.scy, '-.', color = "g")
    ax.plot(T2*fg.scx, Y2*fg.scy, '-.', color = "g")
    ax.plot(T2*fg.scx, A2*fg.scy, '-',  color = "k")
    # add plot to figure
    ax.plot(T3*fg.scx, X3*fg.scy, '-.', color = "b")
    ax.plot(T3*fg.scx, Y3*fg.scy, '-.', color = "b")
    ax.plot(T3*fg.scx, A3*fg.scy, '-',  color = "r")
    # setup legends
    ax.legend(["data", "VX", "VY", "V"])
    # export figure
    D.exportfigure("FIG2")

    #################################################################

    # done
    D.closedocument()
