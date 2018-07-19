import matplotlib.pyplot as plt
from copy import deepcopy
import numpy as np

if __package__ == '':
    __package__ = 'dendrocat'
from .aperture import ellipse, annulus

def specindex(nu1, nu2, f1, alpha):
    return f1*(nu2/nu1)**(alpha)
    
def ffplot(self, rsobj1, rsobj2, apertures=[], specs=[], peak=False,
           label=False, log=True):
    
    if rsobj1.nu > rsobj2.nu:
        rsobj1, rsobj2 = rsobj2, rsobj1
    
    catalog = deepcopy(self.catalog)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    apertures = list(set(apertures)|{ellipse, annulus})
    if specs is None:
        specs = [1, 2, 3]
    
    cols = []
    for aperture in apertures:
        
        cols.append(rsobj1.freq_id+'_'+aperture.__name__+'_rms')
        cols.append(rsobj2.freq_id+'_'+aperture.__name__+'_rms')
        
        if peak:
            cols.append(rsobj1.freq_id+'_'+aperture.__name__+'_peak')
            cols.append(rsobj2.freq_id+'_'+aperture.__name__+'_peak')
        else:
            cols.append(rsobj1.freq_id+'_'+aperture.__name__+'_sum')
            cols.append(rsobj2.freq_id+'_'+aperture.__name__+'_sum')
    
    try:
        index = list(set(range(len(catalog)))^
                     set(np.nonzero(catalog.mask[cols])[0])
                     .union(set(np.where(catalog['rejected']==1)[0])))
    except KeyError:
        for aperture in apertures:
            self.photometer(aperture)
        catalog = deepcopy(self.catalog)
        index = list(set(range(len(catalog)))^
                     set(np.nonzero(catalog.mask[cols])[0])
                     .union(set(np.where(catalog['rejected']==1)[0])))
        
    catalog = catalog[index]
    
    flux1 = []
    flux2 = []
    npix1 = []
    npix2 = []
    err1 = []
    err2 = []
    
    for aperture in list(set(apertures)^set([annulus])):
    
        if peak:
            flux1.append(catalog[rsobj1.freq_id+'_'+aperture.__name__+'_peak'])
            flux2.append(catalog[rsobj2.freq_id+'_'+aperture.__name__+'_peak'])
        else:
            flux1.append(catalog[rsobj1.freq_id+'_'+aperture.__name__+'_sum'])
            flux2.append(catalog[rsobj2.freq_id+'_'+aperture.__name__+'_sum'])
                                 
        npix1.append(catalog[rsobj1.freq_id+'_'+aperture.__name__+'_npix'])
        npix2.append(catalog[rsobj2.freq_id+'_'+aperture.__name__+'_npix'])
    
    err1.append(catalog[rsobj1.freq_id+'_annulus_rms'])
    err2.append(catalog[rsobj2.freq_id+'_annulus_rms'])
    
    marker_labels = catalog['_idx']
    
    xflux = np.linspace(np.min(flux1), np.max(flux1), 10)
    yfluxes = []
    
    for spec in specs:
        yfluxes.append(specindex(rsobj1.nu, rsobj2.nu, xflux, spec))
    
    n_images = len(apertures)
    xplots = int(np.around(np.sqrt(n_images)))
    yplots = xplots
    fig, axes = plt.subplots(ncols=yplots, nrows=xplots, figsize=(12, 12))
    
    for i in range(len(apertures)-1):    
        ax = np.ndarray.flatten(np.array(axes))[i]
        ax.errorbar(flux1[i], flux2[i], xerr=err1[i], yerr=err2[i], fmt='o', 
                    ms=2, alpha=0.75, elinewidth=0.5, color=colors[i], 
                    label='{} Aperture Sums'.format(apertures[i].__name__))
        for j, yflux in enumerate(yfluxes):
            ax.plot(xflux, yflux, '--', color=colors[np.negative(j)], 
                    label='Spectral Index = {}'.format(specs[j]))
        
        ax.set_xticks([])
        ax.set_yticks([])
        
        ax.set_xlim([.6*np.min(flux1), 1.4*np.max(flux1)])
        ax.set_ylim([.1*np.min(flux2), 1.9*np.max(flux2)])
        
        if label:
            for j, label in enumerate(marker_labels):
                ax.annotate(label, 
                            (flux1[i][j], flux2[i][j]), 
                            size=8)
        if peak:
            if log:
                plt.xlabel('Log Peak Flux {}'.format(rsobj1.freq_id))
                plt.ylabel('Log Peak Flux {}'.format(rsobj2.freq_id))
                plt.xscale('log')
                plt.yscale('log')
            else:
                plt.xlabel('Peak Flux {}'.format(rsobj1.freq_id))
                plt.ylabel('Peak Flux {}'.format(rsobj2.freq_id))
        else:
            if log:
                plt.xlabel('Log Flux {}'.format(rsobj1.freq_id))
                plt.ylabel('Log Flux {}'.format(rsobj2.freq_id))
                plt.xscale('log')
                plt.yscale('log')
            else:
                plt.xlabel('Flux {}'.format(rsobj1.freq_id))
                plt.ylabel('Flux {}'.format(rsobj2.freq_id))
        plt.suptitle('{} Flux v. {} Flux'.format(rsobj2.freq_id, 
                                                 rsobj1.freq_id))
        plt.legend()
