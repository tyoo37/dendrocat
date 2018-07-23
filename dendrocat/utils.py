import numpy as np
from radio_beam import Beams
import astropy.units as u
from astropy.table import MaskedColumn, Column, vstack
from astropy.utils.console import ProgressBar
from copy import deepcopy
import warnings

def mask(reg, cutout):
    n = cutout.shape[0]
    mask = reg.to_mask(mode='center')
    return np.array(mask.to_image((n, n)), dtype='int')
    
def specindex(nu1, nu2, f1, alpha):
    return f1*(nu2/nu1)**(alpha) 
    
def findrow(idx, catalog):
    idx = int(idx)
    return catalog[np.where(catalog['_idx'] == idx)]

def rms(x):
    return (np.absolute(np.mean(x**2) - (np.mean(x))**2))**0.5

def check_units(quantity, unit=u.deg):
    if unit.is_equivalent(quantity):
        return quantity.to(unit)
    else:
        return quantity * u.Unit(unit)
        warnings.warn("Assuming quantity is in {}".format(unit))

def commonbeam(major1, minor1, pa1, major2, minor2, pa2):
    """
    Create a smallest bounding ellipse around two other ellipses. 
    Give ellipse dimensions as astropy units quantities.
    """
    major1 = check_units(major1, unit=u.deg)
    minor1 = check_units(minor1, unit=u.deg)
    pa1 = check_units(pa1, unit=u.deg)
    major2 = check_units(major2, unit=u.deg)
    minor2 = check_units(minor2, unit=u.deg)
    pa2 = check_units(pa2, unit=u.deg)
    
    somebeams = Beams([major1.to(u.arcsec), major2.to(u.arcsec)]*u.arcsec, 
                      [minor1.to(u.arcsec), minor2.to(u.arcsec)]*u.arcsec, 
                      [pa1, pa2]*u.deg)
                      
    common = somebeams.common_beam()
    new_major = common._major
    new_minor = common._minor
    new_pa = common._pa
    
    return new_major.to(u.deg), new_minor.to(u.deg), new_pa

def save_regions(catalog, outfile, hide_rejects=True):
    """
    Save a catalog as a a DS9 region file.
    
    Parameters
    ----------
    catalog : astropy.table.Table, RadioSource, or MasterCatalog object
        The catalog or catalog-containing object from which to extract source
        coordinates and ellipse properties.
    outfile : str
        Path to save the region file.
    hide_rejects : bool, optional
        If enabled, rejected sources will not be saved. Default is True
    """
    
    if outfile.split('.')[-1] != 'reg':
        warnings.warn('Invalid or missing file extension. Self-correcting.')
        outfile = outfile.split('.')[0]+'.reg'
    
    if hide_rejects:
        catalog = catalog[np.where(catalog['rejected'] == 0)]
            
    with open(outfile, 'w') as fh:
        fh.write("icrs\n")
        for row in catalog:
            fh.write("ellipse({x_cen}, {y_cen}, {major_fwhm}, " \
                     "{minor_fwhm}, {position_angle}) # text={{{_idx}}}\n"
                     .format(**dict(zip(row.colnames, row))))



def _matcher(obj1, obj2, verbose=True):
    """
    Find sources that match up between two radio objects. 
    
    Parameters
    ----------
    obj1 : rsprocess.RadioSource object or rsprocess.MasterCatalog object
        A catalog with which to compare radio sources.
    obj2 : rsprocess.RadioSource object or rsprocess.MasterCatalog object
        A catalog with which to compare radio sources.
        
    Returns
    ----------
    astropy.table.Table object
    """
  
    all_colnames = set(obj1.catalog.colnames + obj2.catalog.colnames)
    stack = vstack([obj1.catalog, obj2.catalog])
    
    all_colnames.add('_index')
    try:
        stack.add_column(Column(range(len(stack)), name='_index'))
    except ValueError:
        stack['_index'] = range(len(stack))
    stack = stack[sorted(list(all_colnames))]
    
    rejected = np.where(stack['rejected'] == 1)[0]
    
    if verbose:
        print('Combining matches')
        pb = ProgressBar(len(stack) - len(rejected))
    
    i = 0
    while True:
        
        if i == len(stack) - 1:
            break
        
        if i in rejected:
            i += 1
            continue
        
        teststar = stack[i]
        delta_p = vstack([stack[:i], stack[i+1:]]
                         )['_idx', '_index', 'x_cen', 'y_cen']
        delta_p['x_cen'] = np.abs(delta_p['x_cen'] - teststar['x_cen'])                
        delta_p['y_cen'] = np.abs(delta_p['y_cen'] - teststar['y_cen'])
        delta_p.sort('x_cen')
        
        threshold = 1e-5
        found_match = False
        
        dist_col = MaskedColumn(length=len(delta_p), name='dist', 
                                mask=True)
        
        for j in range(10):
            dist_col[j] = np.sqrt(delta_p[j]['x_cen']**2. 
                                  + delta_p[j]['y_cen']**2)            
            if dist_col[j] <= threshold:
                found_match = True
                
        delta_p.add_column(dist_col)
        delta_p.sort('dist')
        
        if found_match:
            match_index = np.where(stack['_index'] == delta_p[0]['_index'])
            match = deepcopy(stack[match_index])
            stack.remove_row(match_index[0][0])
            
            # Find the common bounding ellipse
            new_x_cen = np.average([match['x_cen'], teststar['x_cen']])
            new_y_cen = np.average([match['y_cen'], teststar['y_cen']])
            
            # Find new ellipse properties
            new_maj, new_min, new_pa = commonbeam(
                                         float(match['major_fwhm']), 
                                         float(match['minor_fwhm']), 
                                         float(match['position_angle']),
                                         float(teststar['major_fwhm']),
                                         float(teststar['minor_fwhm']),
                                         float(teststar['position_angle'])
                                         )
            
            # Replace properties of test star
            stack[i]['x_cen'] = new_x_cen       
            stack[i]['y_cen'] = new_y_cen
            stack[i]['major_fwhm'] = new_maj.value
            stack[i]['minor_fwhm'] = new_min.value
            stack[i]['position_angle'] = new_pa.value
    
            # Replace masked data with available values from the match
            for k, masked in enumerate(stack.mask[i]):
                colname = stack.colnames[k]
                if masked:
                    stack[i][colname] = match[colname]
        i += 1
        if verbose:
            pb.update()
    
    # Fill masked detection column fields with 'False'
    for colname in stack.colnames:
        if colname.split('_')[0] == 'detected':
            stack[colname].fill_value = 0
    
    return stack

