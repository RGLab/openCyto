# openCyto 1.3.17

## Enhancements

+ more robust csv parsing
+ register wrapper function for quadGate.seq
+ parse 'subSample' argument from 'gating_args' column in csv template to allow subsampling the data prior to gating.
+ pass ... from quantileGate to quantile function
+ add method for multipleFitlerResult
+ support multi-dimensions in framework by deprecating x,y with channels.

## Bug Fixes:
+ fix the bug that stop.at argument to gating doesn't recognize partial or full path

# openCyto 1.7.1

## Enhancements
+ New API: add_pop function allows users to apply a single gating method to GatingSet without writing the compelete csv template

# openCyto 1.23.5
(Bioconductor version 3.10)

## API Changes

### Simple renaming

* `gate_flowClust_1d` -> `gate_flowclust_1d`
* `gate_flowClust_2d` -> `gate_flowclust_2d`
* `tautStringGate` -> `gate_tautstring`
* `templateGen` -> `gh_generate_template`
* `add_pop_init` -> `gs_add_gating_method_init`
* `add_pop` -> `gs_add_gating_method`
* `remove_pop` -> `gs_remove_gating_method`
* `get.helperGates` -> `gs_get_helpergates`
* `toggle.helperGates` -> `gt_toggle_helpergates`
* `delete.helperGates` -> `gs_delete_helpergates`
* `getNodes` -> `gt_get_nodes`
* `getParent` -> `gt_get_parent`
* `getChildren` -> `gt_get_children`
* `getGate` -> `gt_get_gate`
* `gating` -> `gt_gating`
* `registerPlugins` -> `register_plugins`
* 

### Classes and methods no longer exported

* `registerGatingFunction`
* `gtMethod`
* `gtPopulation`
* `polyFunctions`
* `ppMethod`

## Bug Fixes

* Some minor fixes to `gt_toggle_helpergates`
* Fix "positive" argument to `gate_mindensity` to match doc

# openCyto 1.25.x
(Bioconductor version 3.11)

## Enhancements

* Allow control of mininum number of events for each step of `GatingTemplate` #298

## Bug Fixes

* Handle a few more edge cases in `.improvedMinDensity`
* Added a fix to the density estimate used by `gate_tautstring`