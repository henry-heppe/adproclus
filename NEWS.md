# adproclus (development version)

# adproclus 2.0.0

## Breaking changes
* The default value of the `algorithm` argument in `adproclus()` has been
  changed to `ALS2` from `ALS1`. This is because the `ALS2` algorithm
  often is less time-consuming with similar performance.

## New features
* New model selection tools for full dimensional ADPROCLUS and low 
  dimensional ADPROCLUS. `mselect_adproclus()` 
  and `mselect_adproclus_low_dim()` to estimate a range of models.
  
* New functions to plot the results of the model selection tools. 
  `plot_scree_adpc()` and `plot_scree_adpc_preselected()`.
  
* Automatic model selection with `select_by_CHull()`.

* Compute cluster means on arbitrary data with `cluster_means()`.

## Minor improvements and fixes
* Characteristics of the baseline cluster `Cl0` are now also included
  in the summary output.
  
* Documentation of `plot_vars_by_comp()` now includes the correct
  argument `vars_by_comp` instead of `VarsByComp`.
  
* Give control over the number of decimal places to be displayed in
  the summary output with the `digits` argument in `summary.adproclus()`
  and `cluster_means()`.
  
* Allow partial matches of arguments for plot choice in `plot()`.

* Remove unnecessary arguments passed to the plotting functions
  internally to give user the flexibility to supply them via `...`.
  
* Correct description of the `iterations` argument in `adproclus()`
  and `adproclus_low_dim()`.

* Set default title of plots to `NULL`.

# adproclus 1.0.2
* Using `plot_cluster_network()` now also works when clusters have an 
overlap of 1 observation.

# adproclus 1.0.0

* Initial CRAN submission.
