import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sp
import scipy.optimize as opt
import matplotlib.mlab as mlab
from matplotlib.backends.backend_pdf import PdfPages

nsamp = 2048
Fs = 2e6 #Hz
navg = 700
nelec = 1.3e5
nquant = 1.9e5

noise_elec = 200

nmc=1000

trace = np.loadtxt("test_trace.txt", delimiter=",", skiprows=1)

trace = trace[::5,:]

xvec = np.arange(-1024,1024)
tr = np.interp(xvec, trace[:,0], trace[:,1])
tr[xvec<trace[0,0]] = 0
tr[xvec>trace[-1,0]] = 0

tr = tr * nelec/np.sum(tr)

print(np.sum(tr))

b,a = sp.butter(3,0.01)

def ffn(x, A):
    return A/nelec * tr


tvec = xvec/Fs

# plt.figure()
# plt.plot(tvec*1e6,tr)
# plt.figure()
# plt.plot(tr)
# plt.show()

## noise file from Aldo
dat = np.loadtxt("Gain_6X_Ready/outnoise_tps.csv", skiprows=2, delimiter=",")

plt.figure()  ## figure showing Aldo's noise spectrum
plt.loglog( dat[:,0], dat[:,2], lw=2)

fvec = np.linspace(0, Fs/2, int(nsamp/2)+1)
camps = np.interp( fvec, dat[:,0], np.sqrt(dat[:,2]) )
#plt.loglog(fvec, camps**2)
#plt.show()

def get_aldo_noise(nlev):

    rphase = np.random.rand(int(nsamp/2)+1)*2*np.pi
    rvec = camps * rphase

    noise_dat = np.fft.irfft(rvec)
    noise_dat *= nlev/np.std(noise_dat)
    
    if(False):
        p,f = mlab.psd( noise_dat, Fs = Fs, NFFT=nsamp)
        plt.figure()
        plt.loglog(f,p)

        plt.figure()
        plt.plot( noise_dat)
        
        plt.show()
        
    return noise_dat
        
outvec = []
outvec2 = []
for n in range(nmc):

    tf = 1.0*tr

    noise = get_aldo_noise(noise_elec)
    tf += noise

    ##trapezoidal filter
    outvec.append( np.sum(tf[1186:1226]) )

    ## time domain fit
    bp, bcov = opt.curve_fit(ffn, xvec, tf, p0 = [125000,])
    outvec2.append( bp[0] )

outvec, outvec2 = np.array(outvec), np.array(outvec2)

## now fix the ballistic deficit
outvec *= nelec/np.median(outvec)
outvec2 *= nelec/np.median(outvec2)

h, b = np.histogram(outvec)
h2, b2 = np.histogram(outvec2)
bc = b[:-1] + np.diff(b)/2
bc2 = b2[:-1] + np.diff(b2)/2

def gfit(x, A, m, s):
    return A*np.exp(-(x-m)**2/(2*s**2))

sig=np.sqrt(h)
sig[sig==0]=1
fpts = h>0
bf, bcov = opt.curve_fit(gfit, bc[fpts], h[fpts], sigma=sig[fpts], p0=[nmc, nelec, 5*noise_elec])
sig=np.sqrt(h2)
sig[sig==0]=1
fpts = h2>0
bf2, bcov2 = opt.curve_fit(gfit, bc2[fpts], h2[fpts], sigma=sig[fpts], p0=[nmc, nelec, 5*noise_elec])

xx = np.linspace(nelec-10*noise_elec, nelec+10*noise_elec, int(1e3))

plt.figure()
plt.errorbar(bc, h, yerr=np.sqrt(h), fmt='k.')
plt.errorbar(bc2, h2, yerr=np.sqrt(h2), fmt='r.')
plt.plot(xx, gfit(xx, *bf), color='k', alpha=0.5, label="Integ. (20 $\mu$s), $\sigma = %d\ e$"%bf[2])
plt.plot(xx, gfit(xx, *bf2), color='r', alpha=0.5, label="Fit, $\sigma = %d\ e$"%bf2[2])

plt.legend()

print(np.mean(outvec), np.std( outvec ))

print( np.std( outvec )/np.mean( outvec ) )
print( np.std( outvec )/nquant )

print( np.std( outvec2 )/np.mean( outvec2 ) )
print( np.std( outvec2 )/nquant )




plt.show()
