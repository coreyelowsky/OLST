#!/bin/bash

# sleep for a few seconds
sleep 5

# move log files
mv "${cur_dir}"*${JOB_ID}* ${logs_path}

export feature_index=$((SGE_TASK_ID-1))
export out_path=${data_path}coeffs_${feature_index}

# call python code to run logistic regression
python $log_reg_python_script $feature_index $data_path $morphometrics_path $clusters_path $bootstrap_samples $out_path
