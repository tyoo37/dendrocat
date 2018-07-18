def plot_grid(self, catalog=None, data=None, cutouts=None, cutout_data=None,
              apertures=None, skip=True):
    
    if catalog is None:
        try:
            catalog = self.catalog
        except AttributeError:
            catalog = self.to_catalog()
    
    if data is None:
        data = self.data
    
    # Make sure ellipse and annulus apertures are used
    if apertures is None:
        apertures = [ellipse, annulus]
    else:
        apertures = list(set(apertures+[ellipse, annulus]))
    
    # Get cutouts
    if cutouts is None or cutout_data is None:
        cutouts, cutout_data = self._make_cutouts(catalog=catalog, data=data)
    
    # Get pixels and masks in each aperture
    ap_names = []
    pixels = []
    masks = []
    
    for aperture in apertures:
        some_pixels, a_mask = self.get_pixels(aperture,
                                              catalog=catalog,
                                              cutouts=cutouts,
                                              cutout_data=cutout_data)
        ap_names.append(aperture.__name__)
        pixels.append(some_pixels)
        masks.append(a_mask)
    
    # Find SNR
    ellipse_pix = pixels[ap_names.index('ellipse')]
    annulus_pix = pixels[ap_names.index('annulus')]
    
    snr_vals = self.get_snr(source=ellipse_pix, background=annulus_pix)
    names = np.array(catalog['_idx'])
    rejected = np.array(catalog['rejected'])
    
    if skip:
        accepted_indices = np.where(catalog['rejected'] == 0)[0]
        snr_vals = snr_vals[accepted_indices]
        cutout_data = cutout_data[accepted_indices]
        masks = masks[accepted_indices]
        names = names[accepted_indices]
        rejected = rejected[accepted_indices]
        
    not_nan = ~np.isnan(cutout_data).any(axis=1).any(axis=1)
    snr_vals = snr_vals[not_nan]
    cutout_data = cutout_data[not_nan]
    masks = masks[not_nan]
    names = names[not_nan]
    rejected = rejected[not_nan]
    
    n_images = len(cutout_data)
    xplots = int(np.around(np.sqrt(n_images)))
    yplots = xplots + 1
    gs1 = gs.GridSpec(yplots, xplots, wspace=0.0, hspace=0.0)
    plt.figure(figsize=(9.5, 10))
    
    for i in range(n_images):
        
        image = cutout_data[i]
        plt.subplot(gs1[i])
        
        if rejected[i] == 1:
            plt.imshow(image, origin='lower', cmap='gray')
        else:
            plt.imshow(image, origin='lower')

        for j in range(len(masks[i])):
            plt.imshow(masks[i][j], origin='lower', cmap='gray', 
                       alpha=0.15)
            
        plt.text(0, 0, '{}  SN {:.1f}'.format(names[i], snr_vals[i]), 
                 fontsize=7, color='w')
        plt.xticks([])
        plt.yticks([])
    
    plt.tight_layout()
    plt.show()
