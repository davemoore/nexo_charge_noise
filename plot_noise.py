import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sp
import scipy.optimize as opt
import matplotlib.mlab as mlab
from matplotlib.backends.backend_pdf import PdfPages

pdf = PdfPages('noise_sim_plots.pdf')

nsamp = 2048
Fs = 2e6 #Hz
navg = 700
nelec = 1.3e5
nquant = 1.9e5


nmc=2000

trace = np.loadtxt("test_trace.txt", delimiter=",", skiprows=1)

trace = trace[::5,:]

xvec = np.arange(-1024,1024)
tr = np.interp(xvec, trace[:,0], trace[:,1])
tr[xvec<trace[0,0]] = 0
tr[xvec>trace[-1,0]] = 0

tr = tr * nelec/np.sum(tr)

print(np.sum(tr))

b,a = sp.butter(3, 3e5/(Fs/2))

fitpts = [1100,1250]

def ffn(x, A):
    return A/nelec * tr[fitpts[0]:fitpts[1]]


tvec = xvec/Fs

## noise file from Aldo
dat = np.loadtxt("Gain_6X_Ready/outnoise_tps.csv", skiprows=2, delimiter=",")

plt.figure()  ## figure showing Aldo's noise spectrum
plt.loglog( dat[:,0], dat[:,2], lw=2)
plt.xlabel("Freq [Hz]")
plt.ylabel("Vnoise [V$^2$/Hz]")
plt.xlim(1e2, 5e6)
plt.ylim(1e-14, 1e-10)
plt.title(r"Noise plot from Aldo, 20191206, $\tau_p = 1.2\mu$s")
pdf.savefig()

fvec = np.linspace(0, Fs/2, int(nsamp/2)+1)
camps = np.interp( fvec, dat[:,0], np.sqrt(dat[:,2]) )
#plt.loglog(fvec, camps**2)
#plt.show()

def get_aldo_noise(nlev):

    rphase = np.random.rand(int(nsamp/2)+1)*2*np.pi
    rvec = camps * rphase

    noise_dat = np.fft.irfft(rvec)
    noise_dat *= (nlev/np.sqrt(2))/np.std(noise_dat[100:-100])  #cut out edge effects in ifft, sqrt(2) samples for ENC
    
    if(False):
        p,f = mlab.psd( noise_dat, Fs = Fs, NFFT=nsamp)
        plt.figure()
        plt.loglog(f,p)

        plt.figure()
        plt.plot( noise_dat)
        
        plt.show()
        
    return noise_dat


noise_rms = [100, 150, 200, 250, 300]

noise_type = ["aldo", "white"]

white_vec = []
aldo_vec = []

for noise_val in noise_rms:

    for nt in noise_type:

        outvec = []
        outvec2 = []
        for n in range(nmc):

            tf = 1.0*tr

            if(nt == "aldo"):
                noise = get_aldo_noise(noise_val)
            else:
                ## AA filter
                noise = sp.filtfilt(b,a,noise)
                noise = np.random.randn(nsamp) * noise_val/np.sqrt(2)

                
            tf += noise

            ##trapezoidal filter
            outvec.append( np.sum(tf[1206:1226]) )

            ## time domain fit
            bp, bcov = opt.curve_fit(ffn, xvec[fitpts[0]:fitpts[1]], tf[fitpts[0]:fitpts[1]], p0 = [125000,])
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
        bf, bcov = opt.curve_fit(gfit, bc[fpts], h[fpts], sigma=sig[fpts], p0=[nmc, np.median(outvec), np.std(outvec)])
        sig=np.sqrt(h2)
        sig[sig==0]=1
        fpts = h2>0
        bf2, bcov2 = opt.curve_fit(gfit, bc2[fpts], h2[fpts], sigma=sig[fpts], p0=[nmc, np.median(outvec2), np.std(outvec2)])

        plt.figure()
        plt.errorbar(bc, h, yerr=np.sqrt(h), fmt='k.')
        plt.errorbar(bc2, h2, yerr=np.sqrt(h2), fmt='r.')
        xx = np.linspace(bf[1]-3*np.abs(bf[2]), bf[1]+3*np.abs(bf[2]), int(1e3))
        plt.plot(xx, gfit(xx, *bf), color='k', alpha=0.5, label="Integ. (20 $\mu$s), $\sigma = %d\ e$"%np.abs(bf[2]))
        xx = np.linspace(bf2[1]-3*np.abs(bf2[2]), bf2[1]+3*np.abs(bf2[2]), int(1e3))
        plt.plot(xx, gfit(xx, *bf2), color='r', alpha=0.5, label="Fit, $\sigma = %d\ e$"%np.abs(bf2[2]))
        plt.legend()
        if(nt == "aldo"):
            plt.title("ASIC spec, RMS$=%d\ e$, charge noise (fit), $\sigma=%.2f$%%"%(noise_val, 100*np.abs(bf[2])/nquant))

            aldo_vec.append([noise_val, np.abs(bf[2]), np.abs(bf2[2])])
        else:
            plt.title("White spec, RMS$=%d\ e$, charge noise (fit), $\sigma=%.2f$%%"%(noise_val, 100*np.abs(bf[2])/nquant))
            white_vec.append([noise_val, np.abs(bf[2]), np.abs(bf2[2])])
            
        pdf.savefig()
            
        print(np.mean(outvec), np.std( outvec ))

        print( np.std( outvec )/np.mean( outvec ) )
        print( np.std( outvec )/nquant )

        print( np.std( outvec2 )/np.mean( outvec2 ) )
        print( np.std( outvec2 )/nquant )

aldo_vec = np.array(aldo_vec)
white_vec = np.array(white_vec)

plt.figure()
plt.plot( aldo_vec[:,0], aldo_vec[:,1], 'ko-', label="Trap. filt")
plt.plot( aldo_vec[:,0], aldo_vec[:,2], 'ro-', label='Fit')
plt.xlabel("Noise RMS [$e$]")
plt.ylabel("Single channel charge resolution (noise only), [$e$]")
plt.title("Recon. noise vs charge noise RMS, ASIC spectrum")
plt.legend()
pdf.savefig()

plt.figure()
plt.plot( white_vec[:,0], white_vec[:,1], 'ko-', label="Trap. filt")
plt.plot( white_vec[:,0], white_vec[:,2], 'ro-', label='Fit')
plt.xlabel("Noise RMS [$e$]")
plt.ylabel("Single channel charge resolution (noise only), [$e$]")
plt.title("Recon. noise vs charge noise RMS, white spectrum")
plt.legend()
pdf.savefig()


pdf.close()


