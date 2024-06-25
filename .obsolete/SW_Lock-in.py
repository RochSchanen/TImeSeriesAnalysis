# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 09:20:53 2024

@author: David
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


def downsample(x, n):
    # downsamples data by taking averages
    return np.mean(np.reshape(x, (x.shape[0]//n, n)), 1)



def get_rFFT_PSD(t,x):
    # computes the discrete Fourier transform of a real signal and its power spectral density
    
    dt = t[1]-t[0]                          # time step
    myfft = np.fft.rfft(x)
    myfreqs = np.fft.rfftfreq(x.size, dt)
    PSD = 2*myfft * np.conjugate(myfft) / len(x) * dt
    PSD[0] = myfft[0] * np.conjugate(myfft)[0] / len(x) * dt
    return dt, myfreqs, myfft, PSD



def bandpass_filt(times, signal, f1, f2, order):
    # applies a bandpass filter
    
    dt, myfreqs, myfft, PSD = get_rFFT_PSD(times, signal)
    LPfilter = (f2 / (f2 + 1.0j * myfreqs))**order
    HPfilter = (1.0j * myfreqs / (1.0j * myfreqs + f1))**order
    BPfilter = LPfilter*HPfilter
    filtered_fft = myfft * BPfilter
    filtered_volts = np.fft.irfft(filtered_fft, len(signal))
    filtered_PSD = filtered_fft * np.conjugate(filtered_fft) / len(signal) * dt
    return filtered_volts, filtered_fft, filtered_PSD



def moving_average(x, w):
    # computes moving averages, discarding border values as needed
    return np.convolve(x, np.ones(w), 'valid') / w



def LockinXY(t, sig, f_ref, window_pts):
    # computes software lockin response

    # prepare synthetic signals
    phases = 2*np.pi*t*f_ref
    sx = np.sqrt(2) * np.cos(phases)       # synthetic cosine signal 
    sy = np.sqrt(2) * np.sin(phases)       # synthetic sine signal 

    # computes lockin values in rms
    X = np.sum(sig*sx) / len(sx)
    Y = np.sum(sig*sy) / len(sy)
    return X, Y



def fit_exp(t,A,t0,y0):
    # exponential dacay fitting function
    return A*np.exp(-t/t0) + y0


# read data file

print('Loading data ...')
filename = r'd:\Work\Lancaster\Pillbox\Recording 3.csv'
data = np.loadtxt(filename, delimiter = ',', skiprows = 4)
print('... done.')


times = data[:,0]
volts = data[:,1]

# optional downsampling (uncomment below)
print('Downsampling ...')
lines = data.shape[0]
n = 10                                          # decimation factor
times = downsample(data[:(lines//n)*n,0],n)
volts = downsample(data[:(lines//n)*n,1],n)
print('... done.')


# optional truncation (uncomment below)
#start = 0*2500
#end = 20*2500
#times = times[start:end]
#volts = volts[start:end]

# compute Fourier transform and power spectral density
print('Calculating FFT ...')
dt, myfreqs, myfft, PSD = get_rFFT_PSD(times, volts)
print('... done.')


#optional filtering

use_filter = False
f1 = 10.0                                       # low cutoff, Hz
f2 = 1000.0                                     # high cutoff, Hz
order = 1
if use_filter:
    print('Filtering ...')
    filtered_volts, filtered_fft, filtered_PSD = bandpass_filt(times, volts, f1, f2, order)
    print('... done.')

# signal graphs in time and frequency domain
plt.figure()
plt.title('Decay signal time series')
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.plot(times,volts, label = 'Decay signal')
if use_filter:
    plt.plot(times, filtered_volts, label = 'Filtered decay signal')
plt.legend()

plt.figure()
plt.title('Decay signal power spectral density')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (V$^2$ / Hz)')
plt.xscale('log')
plt.yscale('log')
plt.plot(myfreqs, PSD, 'k-', label = 'Power spectral density')
if use_filter:
    plt.plot(myfreqs, filtered_PSD, 'r-', label = 'Power spectral density after filter')
plt.legend()
plt.show()

if use_filter:
    print('Using bandpass filter')
    volts = filtered_volts


# The code below allows plotting the integral of power spectral density - mainly for checking normalization
# =============================================================================
# df = myfreqs[1]-myfreqs[0]
# PSI = np.cumsum(PSD * df)
# if use_filter:
#     filtered_PSI = np.cumsum(filtered_PSD*myfreqs)
# 
# plt.figure()
# plt.title('Decay signal power spectral density integral')
# plt.xlabel('Frequency (Hz)')
# plt.ylabel('PSI (V$^2$)')
# plt.xscale('log')
# plt.yscale('log')
# plt.plot(myfreqs, PSI, 'k-', label = 'Power spectral integral')
# if use_filter:
#     plt.plot(myfreqs, filtered_PSI, 'r-', label = 'Power spectral integral after filter')
# plt.legend()
# plt.show()
# =============================================================================


# prepare software lock-in computation
f_ref = 95.93                               # reference frequency, Hz (picked by hand as an estimate)
TC = 5.0                                    # time constant, s
window_pts = int(TC/dt)
Ts = []
Xs = []
Ys = []

non_overlap_windows = len(volts)//window_pts
multiplier = 10
N_eval = non_overlap_windows * multiplier
evaluation_pts = np.linspace(window_pts//2+1, len(volts)-(window_pts//2+1), N_eval, dtype = 'int32')

print('Calculating Lock-in response ...')

# compute lock-in readings
for i in range(N_eval):
    vs = volts[evaluation_pts[i] - window_pts//2:evaluation_pts[i] + window_pts//2]
    ts = times[evaluation_pts[i] - window_pts//2:evaluation_pts[i] + window_pts//2]
    X,Y = LockinXY(ts, vs, f_ref, window_pts)
    T = np.mean(ts)
    Xs.append(X)
    Ys.append(Y)
    Ts.append(T)


Ts = np.array(Ts)
Xs = np.array(Xs)
Ys = np.array(Ys)
print('... done.')


# skip some points at beginning
skip_points = 0
Ts = np.array(Ts[skip_points:])
Xs = np.array(Xs[skip_points:])
Ys = np.array(Ys[skip_points:])

Rs = np.sqrt(Xs**2 + Ys**2)
Phases = np.arctan2(Ys,Xs)


# unwrap phase
for i in range(len(Phases)-1):
    if np.abs(Phases[i+1]-Phases[i]) > np.pi:
        direction = np.sign(Phases[i] - Phases[i+1])
        Phases[i+1:] += direction * 2 * np.pi


# determine true signal frequencies from derivative of lockin phase
T_step = Ts[1]-Ts[0]
Freqs = f_ref - np.diff(Phases) / (2 * np.pi * T_step)
Ts_plots = moving_average(Ts,2)                              # for purposes of graphs only (otherwise different number of points)
Rs_plots = moving_average(Rs,2)                              # for purposes of graphs only (otherwise different number of points)


# additional smoothing (also removes points at borders)
print('Additional smoothing ...')
n = multiplier
Freqs = moving_average(Freqs,n)
Ts_plots = moving_average(Ts_plots,n)
Rs_plots = moving_average(Rs_plots,n)
print('... done.')



# exponential fit of amplitude
initial_guess = [np.amax(Rs), np.amax(Ts)/4.0 , 0.0]
params,_ = curve_fit(fit_exp, Ts, Rs, p0 = initial_guess)
print(f'Amplitude fit parameters: A0 = {params[0]:.3e} V, offset = {params[2]:.3e} V, tau = {params[1]:.2f} s')
linewidth = 1 / (2*np.pi * params[1])
print(f'Estimated linewidth: {linewidth*1000:.2f} mHz')


# final graphs
plt.figure()
plt.title(f'Ringdown demodulated @ $f_r$ = {f_ref:.6f} Hz, TC = {TC:.2f} s')
plt.xlabel('Time (s)')
plt.ylabel('Signals (V rms)')
plt.plot(Ts, Xs, 'r-', label = 'Absorption')
plt.plot(Ts, Ys, 'b-', label = 'Dispersion')
plt.plot(Ts, Rs, 'k-', label = 'Amplitude')
plt.plot(Ts, fit_exp(Ts,*params), 'g--', label = 'Exponential fit')
plt.legend()
plt.show()

plt.figure()
plt.title(f'Ringdown demodulated @ $f_r$ = {f_ref:.6f} Hz, TC = {TC:.2f} s')
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')
plt.plot(Ts_plots, Freqs, 'k-', label = 'Frequency')
plt.legend()
plt.show()

plt.figure()
plt.title(f'Ringdown demodulated @ $f_r$ = {f_ref:.6f} Hz, TC = {TC:.2f} s')
plt.xlabel('Amplitude (V rms)')
plt.ylabel('Frequency (Hz)')
plt.plot(Rs_plots, Freqs, 'k-', label = 'Frequency')
plt.legend()
plt.show()
