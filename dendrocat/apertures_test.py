def get_pixels(self, aperture, cutouts, cutout_data):
    # aperture is a function that makes a mask from a source and a cutout.
    # I have several different functions that can be used as 'aperture' here,
    # each of which returns a differently shaped mask.
    
    # Each aperture function should ideally be able to take a few different
    # keyword arguments, like 'radius' for circles and 'padding' and 'width'
    # for annuli.
    
    pix_arrays = []
    
    for i in range(len(cutouts)):
        this_mask = aperture(self.catalog[i], cutouts[i], self)
                            # ^ need some aperture-specific kwargs to go here          
        pix_arrays.append(cutouts[i].data[this_mask.astype('bool')])

    return pix_arrays


def annulus(source, cutout, obj, **kwargs):
    # ^ Need kwargs for padding and width, but need to preserve generality so
    # other aperture functions can still be passed to 'get_pixels'
    
    center = regions.PixCoord(cutout.center_cutout[0], cutout.center_cutout[1])
    inner_r = source['major_fwhm']*u.deg # + padding
    outer_r = inner_r # + width
    
    innerann_reg = regions.CirclePixelRegion(center, inner_r/obj.pixel_scale)
    outerann_reg = regions.CirclePixelRegion(center, outer_r/obj.pixel_scale)
    
    annulus_mask = (mask(outerann_reg, cutout) - mask(innerann_reg, cutout))
    return annulus_mask


# what I want:
get_pixels(annulus, padding=1e-5*u.deg, width=1e-5*u.deg)

# which then calls
this_mask = annulus(self.catalog[i], cutouts[i], self.pixel_scale, 
                    padding=1e-5*u.deg, width=1e-5*u.deg)
# within get_pixels, in order to return a mask with the specified padding and 
# width. 


# I think decorators and wrapping are related to what I'm trying to do, but
# I'm finding it hard to get help online. I basically just need to pass keyword
# arguments from a function into one of that function's arguments, which is 
# also a function.
