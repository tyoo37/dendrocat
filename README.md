### DendroCat

#### To Do
 - [X] Port old modules into object-oriented package
    - [X] Detect
    - [X] Reject
    - [X] Match
    - [X] Flux
    - [X] FFPlot
 - [ ] Test and debug FFPlot
 - [ ] Make `_idx` a unique identifier
 - [ ] Make `freq_id` a unique identifier
 - [ ] Write new module to plot spectral energy distributions
 - [ ] Make default dendrogram parameter generation more robust
     - [ ] Generate default npix parameter based on source size, image scale
 - [X] Generate default annulus parameters based on image scale
 - [ ] Enable dendrogram property saving and loading, using pickle
 - [ ] Sort / remove RadioSource instance attributes when catalog is sorted or rows are removed
 - [ ] Finish documenting methods
 
#### Bugs
 - Large amount of `NaN` SNR values when autorejecting a catalog made from a dendrogram with low min values
    - [SOLVED] Caused by empty arrays for annulus pixels. Solution is to adjust default annulus width and padding to ensure annulus region isn't infinitely thin.
