

def pixel_scale(imagefile):
    from astropy.io import fits
    from astropy.wcs import WCS

    fitsdata = fits.open(imagefile)
    hdr = fits.getheader(imagefile)
    wcsHR = WCS(hdr,naxis=2)
    scale = wcsHR.proj_plane_pixel_scales()[0]

    return scale
    

