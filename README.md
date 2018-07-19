### DendroCat

#### To Do
 - [X] Port old modules into object-oriented package
    - [X] Detect
    - [X] Reject
    - [X] Match
    - [X] Flux
    - [X] FFPlot
 - [X] Test and debug FFPlot
 - [X] Write new module to plot spectral energy distributions
 - [ ] Add upper limits for detections with 2~3 times higher noise than signal
    - [ ] Only plot noise level (downward arrow)
 - [ ] Preserve original dendrogram _idx, rework to use a different unique identifier for everything else
 - [ ] Save and load object information with pickle
 - [ ] Implement circular apertures of different radii
    - [ ] Function name changes with radius
 - [ ] Improve rejection algorithm
    - [ ] Reject high eccentricity source ellipses with off-center peak fluxes
    - [ ] Reject large areas where the peak flux isn't much different than the median (indicator of a large noise pocket)
 - [ ] Make `freq_id` a unique identifier
 - [ ] Find an easier way to grab a specific source by name (indexing astropy tables is clumsy)
 - [ ] Make default dendrogram parameter generation more robust
     - [ ] Generate default npix parameter based on source size, image scale
 - [X] Generate default annulus parameters based on image scale
 - [ ] Sort / remove RadioSource instance attributes when catalog is sorted or rows are removed
 - [ ] Finish documenting methods
 
#### Bugs
 - Large amount of `NaN` SNR values when autorejecting a catalog made from a dendrogram with low min values
    - [SOLVED] Caused by empty arrays for annulus pixels. Solution is to adjust default annulus width and padding to ensure annulus region isn't infinitely thin.
