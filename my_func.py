##### Import required packages #####
# standard packages
import os
import shutil
import time
import pandas as pd
import numpy as np
import sys
# from within environment_calibration_common submodule
from helpers import load_coordinator_df
from translate_parameters import translate_parameters, get_initial_samples
from run_sims import submit_sim
from get_eradication import get_eradication
from compare_to_data.run_full_comparison import compute_scores_across_site
# from source 'simulations' directory
sys.path.append("../simulations")
import manifest as manifest


def my_func(X,wdir):
    coord_df=load_coordinator_df()
    site = coord_df.at['site','value']
    sites =[site]
    n_sims = int(coord_df.at['nSims','value'])
    # Supply parameters to X
    get_eradication(manifest.use_local_eradication)
    param_key=pd.read_csv("parameter_key.csv")
    df = pd.DataFrame({'parameter':[], 'unit_value': [], 'emod_value':[], 'type':[], 'param_set':[]})
    i=1
    for x in X:
        a = translate_parameters(param_key,x,i)
        a['param_set'] = np.repeat(i,len(a))
        i=i+1
        df = pd.concat([df,a])
    df.to_csv(f"{wdir}/translated_params.csv")
  
    for i, my_site in enumerate(sites):
        if os.path.exists(os.path.join(manifest.simulation_output_filepath,my_site)):
            shutil.rmtree(os.path.join(manifest.simulation_output_filepath,my_site))
        submit_sim(site=my_site, nSims=n_sims, X=df)
      
  
    while True:
        outputs = []
        for my_site in sites:
            outputs.append(os.path.exists(os.path.join(manifest.simulation_output_filepath,my_site,'finished.txt')))
        if all(outputs):#os.path.exists(manifest.simulation_output_filepath,my_site): 
            for my_site in sites:
                Y = compute_scores_across_site(my_site)
            break
        time.sleep(120)   
    
    return(Y)

if __name__ == '__main__':
    param_key=pd.read_csv("parameter_key.csv")
    X = get_initial_samples(param_key, 100)
    X = X[1:5]
    my_func(X)
