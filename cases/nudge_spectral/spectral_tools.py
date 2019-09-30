import matplotlib.pyplot as pl
import numpy as np
from numba import jit

pl.close('all')
pl.ion()

np.random.seed(3)

def calc_radial_distance(itot, jtot):

    j0 = jtot//2
    jj,ii = np.indices((jtot, itot//2+1))
    jj -= j0
    return np.sqrt(ii**2 + jj**2)


def calc_spectral_mask(itot, jtot, k0):

    dist = calc_radial_distance(itot, jtot)
    mask = dist <= k0

    return mask.astype(int)


def spectral_filter_2d(field, k0, method='low'):

    jtot = field.shape[0]
    itot = field.shape[1]

    # 2D forward FFT
    f = np.fft.rfft2(field)
    f = np.fft.fftshift(f, axes=0)

    # Create/define/apply wavenumber mask
    mask = calc_spectral_mask(itot, jtot, k0)
    
    if method == 'high':
        mask = 1.-mask

    f *= mask

    # Shift back, and calculate/return 2D backward FFT
    return np.fft.irfft2(np.fft.ifftshift(f, axes=0)), mask


def spectral_blend_2d(field1, field2, k0):

    jtot = field1.shape[0]
    itot = field1.shape[1]

    # 2D forward FFT
    f1 = np.fft.rfft2(field1)
    f2 = np.fft.rfft2(field2)

    f1 = np.fft.fftshift(f1, axes=0)
    f2 = np.fft.fftshift(f2, axes=0)

    # Create/define wavenumber mask
    mask1 = calc_spectral_mask(itot, jtot, k0)
    mask2 = 1-mask1

    # Blend
    f = mask1*f1 + mask2*f2

    # Shift back, and calculate/return 2D backward FFT
    return np.fft.irfft2(np.fft.ifftshift(f, axes=0))


def spectral_noise(itot, jtot, k, fac, ks, ampl):

    # Generate 2D spectrum with random noise
    rand = 2*np.pi*np.random.rand(jtot, itot//2+1)
    fft_rand = np.cos(rand) + 1j*np.sin(rand)

    # Shift small wavenumbers to center
    fft_rand = np.fft.fftshift(fft_rand, axes=0)

    # Calculate the radial wave numbers
    l = calc_radial_distance(itot, jtot)

    # Filter on radial wave number using a Gaussian function
    factor = np.zeros_like(rand)
    for i in range(len(k)):
        factor += fac[i]*np.exp(-(l-k[i])**2. / (2.*ks[i]))

    # Create the filtered field
    afft = factor*fft_rand

    # Inverse shift & 2D FFT
    a = np.fft.irfft2(np.fft.ifftshift(afft, axes=0))

    # Scale to requested amplitude
    a *= ampl / (a.max() - a.min())

    return a


def radial_spectrum(field, kmax=None):

    jtot = field.shape[0]
    itot = field.shape[1]

    # Calculate 2D FFT, remove mean, and shift output
    fft = np.fft.rfft2(field) / (itot*jtot)
    fft[0,0] = 0.
    fft = np.fft.fftshift(fft, axes=0)

    # Radial distance
    l = calc_radial_distance(itot, jtot)

    # Calculate spectral power
    power = fft.real**2 + fft.imag**2
    power[:,1:-1] *= 2

    # Bin on wavenumber
    if kmax is None:
        kmax = int(np.sqrt((itot//2+1)**2 + (jtot//2+1)**2))
    bin_edge = np.arange(0, kmax+1, 1)

    n, bin_edges = np.histogram(l.flatten(), bin_edge, weights=power.flatten())
    bins = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    return bins, n


if __name__ == '__main__':

    itot = 512
    jtot = 512

    xsize = 2*np.pi
    ysize = 2*np.pi

    dx = xsize/itot
    dy = ysize/jtot

    x  = np.linspace(dx/2, xsize, itot)
    y  = np.linspace(dy/2, ysize, jtot)

    kx = np.arange(itot/2+1)
    ky = np.arange(jtot/2+1)


    if True:
        z1 = spectral_noise(itot, jtot, k=[2],  fac=[1], ks=[3], ampl=3)
        z2 = spectral_noise(itot, jtot, k=[20], fac=[1], ks=[10], ampl=1)
        z3 = spectral_blend_2d(z1, z2, 10)

        if True:
            f = 200 / xsize


            pl.figure(figsize=(11,3.8))
            pl.subplot(131, aspect='equal')
            pl.title('Large scale', loc='left')
            pl.pcolormesh(f*x, f*y, z1, cmap=pl.cm.RdBu_r)
            pl.xlabel(r'$x$ (km)')
            pl.ylabel(r'$y$ (km)')

            pl.subplot(132, aspect='equal')
            pl.title('Small scale', loc='left')
            pl.pcolormesh(f*x, f*y, z2, cmap=pl.cm.RdBu_r)
            pl.xlabel(r'$x$ (km)')
            pl.ylabel(r'$y$ (km)')

            pl.subplot(133, aspect='equal')
            pl.title('Blended', loc='left')
            pl.pcolormesh(f*x, f*y, z3, cmap=pl.cm.RdBu_r)
            pl.xlabel(r'$x$ (km)')
            pl.ylabel(r'$y$ (km)')

            pl.tight_layout()

            pl.figure()
            pl.plot(x, z1[jtot//2,:], 'k-', label='Large scale')
            pl.plot(x, z2[jtot//2,:], 'k-', label='Small scale', dashes=[5,2])
            pl.plot(x, z3[jtot//2,:], 'r-', label='Blended')
            pl.legend()



        b1, f1 = radial_spectrum(z1, kmax=30)
        b2, f2 = radial_spectrum(z2, kmax=30)
        b3, f3 = radial_spectrum(z3, kmax=30)


        pl.figure()
        pl.semilogy(b1, f1, 'k--x', label='large scale')
        pl.semilogy(b2, f2, 'k--o', label='small scale')
        pl.semilogy(b3, f3, 'r-', linewidth=1.5, label='blended')
        pl.legend()
        pl.xlabel('K')
        pl.ylabel('power')




#for n in range(time.size):
#    qtpath = data.variables['qtpath'][n,:,:]
#    qtpath_fft = np.fft.rfft2(qtpath) / (x.size*y.size)
#    qtpath_fft[0,0] = 0.
#    qtpath_fft_power = abs(qtpath_fft)**2
#    qtpath_fft_power[:,1:-1] *= 2
#    qtpath_spectra[n,:,:] = qtpath_fft_power

#bin_seq = (2.*np.pi/L)*np.arange(.5, k.size+10, 2)
#n, bin_edges = np.histogram(K.flatten(), bin_seq, weights=(qtpath_spectra[0,:,:]).flatten())
#pdfy = np.zeros((time.size, n.size))
#pdfx = 0.5*(bin_edges[0:-1] + bin_edges[1:])
#for i in range(time.size):
#    n, bins = np.histogram(K.flatten(), bin_seq, weights=(qtpath_spectra[i,:,:]).flatten())
#    pdfy[i,:] = n[:]





    if False:
        z = np.zeros((jtot, itot))

        for j in range(jtot):
            z[j,:] += np.cos(x*1)
            z[j,:] += np.cos(x*5)
            z[j,:] += np.cos(x*10)

        zz1, mask1 = spectral_filter_2d(z, 5, 'low')
        zz2, mask2 = spectral_filter_2d(z, 5, 'high')

        pl.figure()
        pl.subplot(231)
        pl.pcolormesh(x, y, z)
        pl.colorbar()

        pl.subplot(232)
        pl.pcolormesh(x, y, zz1)
        pl.colorbar()

        pl.subplot(233)
        pl.pcolormesh(kx-0.5, np.arange(jtot)-0.5, mask1)
        pl.colorbar()

        pl.subplot(235)
        pl.pcolormesh(x, y, zz2)
        pl.colorbar()

        pl.subplot(236)
        pl.pcolormesh(kx-0.5, np.arange(jtot)-0.5, mask2)
        pl.colorbar()


        pl.figure()
        pl.plot(x, z  [jtot//2, :], label='orig')
        pl.plot(x, zz1[jtot//2, :], label='low-pass')
        pl.plot(x, zz2[jtot//2, :], label='high-pass')
        pl.legend()


