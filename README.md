# environment_calibration_common

common files 'under the hood' for environment calibration workflow. 

#### Some helper files:

* clean_all.py 
* get_eradication.py - Downloads new eradication if manifest doesn't dictate to use a local copy
* get_version.py - Gives version of eradication being used if needed
* load_inputs.py - Loads the coordinator and sites requested
* my_func.py
* run_analyzers.py - Dictates which analyzers to run for which sites based on simulation coordinator. Makes finished.txt when all analyzers for a site are done running.
* run_sims.py - Simulation submision script that manages bringing all of the helpers together to create simulations. During simulation submission we also submit a scheduled waiting script and analyzer script. 
* translate_parameters.py - script for converting locations in unit parameter space to emod-compatible values, and vice-versa, according to 'parameter_key.csv'
* utils_slurm.py - slurm-specific helper functions for chain job submission on QUEST
* helpers.py (modified per IIVT) - Contains general helper functions for simulation setup. If trying to change the core workings of sims, config files, etc, this is the place to do it.
  
### analyzers

* analyzer_collection.py
* analyze.py

### compare_to_data

* calculate_all_scores.py
* run_full_comparison.py
