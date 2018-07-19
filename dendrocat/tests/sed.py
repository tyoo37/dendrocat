import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np


def findrow(idx, catalog):
    return catalog[np.where(catalog['_idx'] == idx)]

def plot_sed(self, idx=None, aperture=None, alphas=None, peak=False, log=True, 
             outfile=None):
    '''
    Plot a spectral energy distribution for a specific source in the catalog.
    
    Parameters
    ----------
    idx : str
        The identifier used to specify a row in the MasterCatalog.
    tags : list of str, optional
        Keywords to search for in catalog column names, to get information 
        about which aperture shape and flux extraction method (i.e., sum, peak) 
        to use.
    alphas : list of float, optional
        Spectral indices to plot under flux data.
    log : bool, optional
        If enabled, spectral energy distribution will be plotted on a log scale
    '''
    

        
    
    if callable(aperture) is True:
        aperture = aperture.__name__
    if aperture is None:
        aperture = 'ellipse'
    method = ['peak' if peak else 'sum'][0]
    
    row = findrow(idx, self.catalog)
    cols = [s for s in row.colnames if aperture in s and method in s]
    err_cols = [e for e in row.colnames if 'annulus' in e and 'rms' in e]
    freq_ids = [col.split('_')[0] for col in cols]
    nus = []
    for obj in self.__dict__.values():
        try:
            for freq_id in freq_ids:
                if obj.freq_id == freq_id:
                    nus.append(obj.nu.to(u.GHz).value)
        except AttributeError:
            pass
    fluxes = [row[0][col] for col in cols]
    errs = [row[0][err_col] for err_col in err_cols]
    x = np.linspace(0.8*np.min(nus), 1.1*np.max(nus), 100)
    ys = []
    
    if alphas:
        k = int(np.floor(len(fluxes)/2))
        for a in alphas:
            constant = fluxes[k]/(nus[k]**a)
            ys.append(constant*(x**a))
    
    fig, ax = plt.subplots()              
    ax.errorbar(nus, fluxes, yerr=errs, fmt='o', ms=2, elinewidth=0.75, 
                color='k', label='{} aperture {}'.format(aperture, method), 
                zorder=3)
    if ys:
        for i, y in enumerate(ys):                     
            ax.plot(x, y, label=r'$\alpha$ = {}'.format(alphas[i], zorder=2))
     
    if log is True:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Log Frequency (GHz)')
        if peak is True:
            ax.set_ylabel('Log Peak Flux (Jy)')
        else:
            ax.set_ylabel('Log Flux (Jy)')
    else:
        ax.set_xlabel('Frequency (GHz)')
        if peak is True:
            ax.set_ylabel('Peak Flux (Jy)')
        else:
            ax.set_ylabel('Flux (Jy)')
    
    ax.set_xticks(nus, ['{} GHz'.format(nu) for nu in nus])
    ax.set_title('Spectral Energy Distribution for {}'.format(idx))
    ax.legend()
    
    if outfile is not None:
        ax.savefig(outfile, dpi=300, bbox_inches='tight')
