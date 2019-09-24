import matplotlib.pyplot as pl
import numpy as np
from numba import jit

pl.close('all')
pl.ion()

np.random.seed(3)

def radial_distance(itot, jtot):

    j0 = jtot//2
    jj,ii = np.indices((jtot, itot//2+1))
    jj -= j0    # Shift center
    return np.sqrt(ii**2 + jj**2)


def spectral_mask(itot, jtot, k0):

    dist = radial_distance(itot, jtot)
    mask = dist <= k0

    return mask.astype(int)


def spectral_filt_2d(field, k0, method='low'):

    jtot = field.shape[0]
    itot = field.shape[1]

    # 2D forward FFT
    f = np.fft.rfft2(field)

    # Result of 2D FFT = real fft over last dimension, complex over 
    # first dimension. Shift the complex part for convenience
    f = np.fft.fftshift(f, axes=0)

    # Create/define wavenumber mask
    mask = spectral_mask(itot, jtot, k0)
    
    if method == 'high':
        mask = 1.-mask

    # Apply mask
    f *= mask

    # Shift back, and calculate/return 2D backward FFT
    return np.fft.irfft2(np.fft.ifftshift(f, axes=0)), mask


def spectral_blend_2d(field1, field2, k0):

    jtot = field1.shape[0]
    itot = field1.shape[1]

    # 2D forward FFT
    f1 = np.fft.rfft2(field1)
    f2 = np.fft.rfft2(field2)

    # Result of 2D FFT = real fft over last dimension, complex over 
    # first dimension. Shift the complex part for convenience
    f1 = np.fft.fftshift(f1, axes=0)
    f2 = np.fft.fftshift(f2, axes=0)

    # Create/define wavenumber mask
    mask1 = spectral_mask(itot, jtot, k0)
    mask2 = 1-mask1

    # Blend
    f = mask1*f1 + mask2*f2

    # Shift back, and calculate/return 2D backward FFT
    return np.fft.irfft2(np.fft.ifftshift(f, axes=0))


def spectral_noise(itot, jtot, k, fac, ks, ampl):

    # Generate 2D spectrum with random noise
    arnd = 2*np.pi*np.random.rand(jtot, itot//2+1)
    afftrnd = np.cos(arnd) + 1j*np.sin(arnd)

    # Shift small wavenumbers to center
    afftrnd = np.fft.fftshift(afftrnd, axes=0)

    # Calculate the radial wave numbers
    l = radial_distance(itot, jtot)

    # Filter on radial wave number using a Gaussian function
    if isinstance(k, list):
        factor = np.zeros_like(arnd)
        for i in range(len(k)):
            factor += fac[i]*np.exp(-(l-k[i])**2. / (2.*ks[i]))
    else:
        factor  = fac*np.exp(-(l-k)**2. / (2.*ks))

    # Create the filtered field
    afft = factor*afftrnd

    # Inverse shift & 2D FFT
    a = np.fft.irfft2(np.fft.ifftshift(afft, axes=0))

    # Scale to requested amplitude
    a *= ampl / (a.max() - a.min())

    return a





if True:

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
        z1 = spectral_noise(itot, jtot, k=2, fac=1, ks=1, ampl=2)
        z2 = spectral_noise(itot, jtot, k=[4,20], fac=[2,1], ks=[3,1], ampl=1)
        z3 = spectral_blend_2d(z1, z2, 5)

        if True:
            pl.figure(figsize=(10,4))
            pl.subplot(131, aspect='equal')
            pl.pcolormesh(x, y, z1, cmap=pl.cm.RdBu_r)

            pl.subplot(132, aspect='equal')
            pl.pcolormesh(x, y, z2, cmap=pl.cm.RdBu_r)

            pl.subplot(133, aspect='equal')
            pl.pcolormesh(x, y, z3, cmap=pl.cm.RdBu_r)

            pl.tight_layout()

            pl.figure()
            pl.plot(x, z1[jtot//2,:], 'k-', label='Large scale')
            pl.plot(x, z2[jtot//2,:], 'r-', label='Small scale')
            pl.plot(x, z3[jtot//2,:], 'k--', label='Blended')






    if False:
        z = np.zeros((jtot, itot))

        for j in range(jtot):
            z[j,:] += np.cos(x*1)
            z[j,:] += np.cos(x*5)
            z[j,:] += np.cos(x*10)

        zz1, mask1 = spectral_filt_2d(z, 5, 'low')
        zz2, mask2 = spectral_filt_2d(z, 5, 'high')

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


