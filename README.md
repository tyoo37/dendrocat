### DendroCat
#### REBUILD BRANCH
`master` made use of some quick hacks that aren't entirely good coding practice. This branch removes those hacks and attempts to rebuild their ends through better means.

#### To Do
 - [X] Port old modules into object-oriented package
    - [X] Detect
    - [X] Reject
    - [X] Match
    - [X] Flux
    - [X] FFPlot
 - [X] Test and debug FFPlot
 - [X] Write new module to plot spectral energy distributions
 - [X] Add upper limits for detections with 2~3 times higher noise than signal
    - [X] Only plot noise level (downward arrow)
    - [X] Plot three arrows for 1, 2, and 3 sigma
 - [X] Preserve original dendrogram _idx, rework to use a different unique identifier for everything else
 - [X] Match sources with external catalogs
 - [X] Add ability to use external catalog fluxes to extend SED plots
 - [X] Support for adding custom sources to any source catalog, by vstacking a source table
    - [X] `add_sources` method for RadioSource and MasterCatalog
 - [ ] Add two new apertures using the radio beam
    - [ ] Beam centered on peak flux (point-source aperture)
    - [ ] Ellipse with dimensions 2 x beam centered on peak flux (point-source background)
 - [ ] Combined aperture grid and SED plots
 - [ ] `match_external` takes an argument for source names -- if a match is made, the source name is replaced
 - [ ] MasterCatalog and RadioSource methods to grab a specific source by name, idx, etc (i.e., `mc.grab('w51e2')`)
 - [ ] Non-rejected catalog is an attribute of the RadioSource or MasterCatalog object (i.e., `rs1.nonrejected`)
 - [ ] Improve rejection algorithm
    - [ ] Reject high eccentricity source ellipses with off-center peak fluxes
    - [ ] Reject large areas where the peak flux isn't much different than the median (indicator of a large noise pocket)
 - [ ] Make `_name` a unique identifier
 - [ ] Finish documenting methods
 
#### Fixed Bugs
 - Issue: Large amount of `NaN` SNR values when autorejecting a catalog made from a dendrogram with low min values
    - Solution: Caused by empty arrays for annulus pixels. Solution is to adjust default annulus width and padding to ensure annulus region isn't infinitely thin.
 - Issue: Elliptical flux sums are sometimes lower than the peak flux
    - Solution The ellipse may contain negative flux values, if it includes a decent portion of the background. In this case, the sum can be less than the peak flux. Temporary solution is to only sum positive values.
 - Issue: Ellipses in plot_grid still appear to be smaller than those in DS9
    - Solution: Made same correction as before -- multiply all astropy ellipse region dimensions by 2
 - Issue: Some sources disappear when matched
    - Solution: Matching algorithm was matching accepted sources with rejected sources, overwriting the accepted sources so they don't appear in the master catalog. Fixed by only matching accepted sources to each other.

#### Unresolved Bugs
 - Issue: Bright sources with contaminated annululi are mistaken for upper limit detections due to their low SN
