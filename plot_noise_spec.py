import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import matplotlib.mlab as mlab

## noise file from Aldo
dat = np.loadtxt("Gain_6X_Ready/outnoise_tps.csv", skiprows=2, delimiter=",")

plt.figure()  ## figure showing Aldo's noise spectrum
plt.loglog( dat[:,0], dat[:,2], lw=2)

## Function for 1/f noise
def ffn(x,A):
    return A/x

## Function for f^n noise
def ffn2(x,A,n):
    return A * x**n

## Function for total noise spectrum
def ffntot(x, A1, A2, A3,n):
    return np.abs(A1) + np.abs(A2)/x + np.abs(A3) * x**n

fpts = dat[:,0] < 2e3

p = [0.3e-12,5e-9,6e-17,0.9]  #best fit parameters to describe Aldo's spectrum

xx2 = np.linspace(1e3, 1e5)  ## f values for white
xx4 = np.linspace(5e3, 2e5) ## f values for white + f^n
xx = np.linspace( dat[fpts,0][0], 10*dat[fpts,0][-1], 1e3) ## f values for 1/f

## plots of spectral components
plt.plot( xx2, p[0]*np.ones_like(xx2), 'c--', label="white")
plt.plot( xx4, p[0]*np.ones_like(xx4) + ffn2(xx4, p[2], p[3]), 'g--', label="white + f$^{0.9}$")
plt.plot( xx, ffn(xx, p[1]), 'm--', label="1/f")

xx3 = np.linspace(dat[fpts,0][0], 2e5, 1e2) ## f values for full spectrum
plt.plot(xx3, ffntot(xx3, *p), 'r--', label="sum" )

plt.xlabel("Freq [Hz]")
plt.ylabel("Vnoise [V$^2$/Hz]")
plt.legend()
plt.xlim(1e2, 1e9)
plt.ylim(1e-14, 1e-10)
plt.title(r"Noise plot from Aldo, 20191206, $\tau_p = 1.2\mu$s")
plt.savefig("noise_spec_1_2us.pdf")
plt.close('all')

## now simulate some time streams to check the noise

xf = np.linspace( dat[0,0], 1e6, 2048)  ## frequencies to simulate noise for, up to nyquist at 2 MHz
cd = np.interp( xf, dat[:,0], dat[:,2] ) ## interpolate Aldo's spectrum at these frequencies

tint = [0.5, 1,5,10, 20, 40, 100, 200] #Integration times for analysis [us]  

nmc = 10000 # number of noise traces to generate
out_dat = []
charge_spec = np.zeros( 1025 )
curr_spec = np.zeros( 1025 )
for i in range(nmc):

    cdat = np.sqrt( cd )*6e4 ## amplitude spectrum with proper norm
    cph = np.random.rand( len( cd ) )*2*np.pi ## random phase

    ts = np.fft.irfft( cdat*np.exp(1j*cph) )  ## current waveform
    tsc = np.cumsum( ts ) ## charge waveform found by integrating current waveform
                          ## note that the 1/f noise really shows up here -- claim it doesn't matter
                          ## whether this is integrated with software or an analog capacitor after
                          ## band-limiting the noise with the anti-aliasing filter (and neglecting
                          ## digitization noise

    p,f = mlab.psd( tsc, NFFT=2048, Fs = 2e6 )                      
    charge_spec += p
    p,f = mlab.psd( ts, NFFT=2048, Fs = 2e6 )                      
    curr_spec += p 
                          
    ## enable to double check the waveforms, and also that noise spectrum matches what was intended
    #plt.figure()
    #plt.plot(tsc)
    #p,f = mlab.psd( ts, NFFT=2048, Fs = 2e6 )
    #plt.loglog(f, p)
    #plt.loglog(dat[:,0], dat[:,2])    
    #plt.show()

    ## measure noise versus integration time
    cn = []
    for ti in tint:
        ## factor of 2 to get samples from us to index (at 2 MHz sampling)
        ## the following gives the noise on the baseline -- i.e. the noise on a charge
        ## measurement integrated over a given time window
        ca = np.mean( tsc[int(2*ti):int(4*ti)] ) - np.mean( tsc[:int(2*ti)] ) 
        cn.append(ca)
    out_dat.append(cn)

out_dat = np.array(out_dat) ## array holding noise realizations vs integration time

## plot the power spectra of charge and current to compare
charge_spec /= nmc
curr_spec /= nmc

AD8655 = np.loadtxt("ralph_AD8655.txt", skiprows=1, delimiter=",")
LTC6268 = np.loadtxt("ralph_LTC6268.txt", skiprows=1, delimiter=",")
plt.figure()
plt.loglog(dat[:,0], dat[:,2], label="Aldo's spec.")
plt.loglog(AD8655[:,0]*1e3, (10**(AD8655[:,1]/20))**2, label="Ralph, AD8655")
plt.loglog(LTC6268[:,0]*1e3, (10**(LTC6268[:,1]/20))**2, label="Ralph, LTC6268")  
plt.loglog( f, curr_spec, label="Sim. curr.")
plt.loglog( f, charge_spec, label="Sim. charge")
plt.legend(loc="upper right")
plt.xlabel("Freq. (Hz)")
plt.ylabel("Noise PSD (V$^2$/Hz)")
plt.xlim([1e3, 2e6])
plt.ylim([1e-15, 1e-5])



## histogram up the noise realizations to get the distribution
## In the end we record the sigma of these distributions vs integration time
plt.figure()
noise_dat = []
for i,ti in enumerate(tint):

    hh, be = np.histogram( out_dat[:,i], bins=100, range=(-0.2, 0.2) )

    plt.step(be[:-1], hh, where='post', label=r"$\tau_{int} = %d \mu s$"%ti)

    noise_dat.append( [ti, np.std(out_dat[:,i])] )

plt.legend()


## now make a plot of the noise expected for charge measurement vs integration time,
## and compare to the expectation for white noise, integral grows like ~sqrt(tau)
noise_dat = np.array(noise_dat)

plt.figure()
plt.plot(noise_dat[:,0], noise_dat[:,1], 'ko-', label="Sim. noise")
xx = np.linspace(0, 200, 1e3)
scale_fac = noise_dat[0,1]/np.sqrt(0.5)
plt.plot( xx, np.sqrt(xx) *scale_fac, 'r', label=r"$1/\sqrt{\tau}$")

plt.xlabel(r"Integration time, $\tau$ [$\mu$s]")
plt.ylabel("Baseline noise [arb. units]")
plt.xlim([0,200])
plt.ylim([0, 0.07])
plt.legend(loc="upper left")
plt.savefig("noise_vs_int_time.pdf")

plt.show()
