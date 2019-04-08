import numpy as np


class Soil_type:
    def __init__(self, phi_sat, phi_fc, phi_wp, name):
        self.phi_sat = phi_sat
        self.phi_fc  = phi_fc
        self.phi_wp  = phi_wp
        self.name = name

    def calc_f2(self, phi_in):

        f2 = (phi_in - self.phi_wp) / (self.phi_fc - self.phi_wp)

        if np.any(f2 > 1):
            print('Warning: phi > phi_fc')
        if np.any(f2 < 0):
            print('Warning: phi < phi_wp')

        f2 = np.minimum(1, np.maximum(0, f2))

        return f2

    def rescale(self, phi_in, new_type):

        f2 = self.calc_f2(phi_in)
        phi_out = new_type.phi_wp + f2 * (new_type.phi_fc - new_type.phi_wp)

        return phi_out


med_fine = Soil_type(0.430, 0.383, 0.133, 'ECMWF medium fine')
fine     = Soil_type(0.520, 0.448, 0.279, 'ECMWF fine')
wosten   = Soil_type(0.590, 0.528, 0.320, 'Wosten')


if __name__ == '__main__':

    import matplotlib.pyplot as pl
    pl.close('all')


    phi  = np.array([0.3106892,  0.2866719,  0.29377648, 0.35078132])
    phi2 = med_fine.rescale(phi, fine)
    phi3 = med_fine.rescale(phi, wosten)
    
    f21 = med_fine.calc_f2(phi)
    f22 = wosten.  calc_f2(phi2)
    
    cc = pl.cm.bwr(np.linspace(0,1,4))
    
    pl.figure()
    ax=pl.subplot(111)
    pl.scatter(np.ones(3)*1, [med_fine.phi_wp, med_fine.phi_fc, med_fine.phi_sat], marker='x', color='k')
    pl.scatter(np.ones(3)*2, [fine.    phi_wp, fine.    phi_fc, fine.    phi_sat], marker='x', color='k')
    pl.scatter(np.ones(3)*3, [wosten.  phi_wp, wosten.  phi_fc, wosten.  phi_sat], marker='x', color='k')
    for i in range(4):
        pl.plot([1,3], [phi[i], phi[i]], color=cc[i], label='L{} original'.format(i+1))
    for i in range(4):
        pl.plot([1,2,3], [phi[i], phi2[i], phi3[i]], '--', color=cc[i], label='L{} scaled'.format(i+1))
    ax.set_xticks(np.arange(1,4))
    ax.set_xticklabels(['ECMWF medium_fine', 'ECMWF fine', 'Wosten'])
    pl.ylabel('phi (m3/m3)')
    pl.grid()
    pl.legend(ncol=2)
