#!/bin/bash

# wait a few seconds to make sure job will show up in qstat
sleep 10

# sleep until the job is no longer in qstat output (jobs are complete)
qstat_output=`qstat | grep $job_id`

while [[ "$qstat_output" == *"$job_id"* ]]
do
  	sleep 1
	qstat_output=`qstat | grep $job_id` 
done

# run python script to merge
python $merge_coeffs_python_script $data_path $num_features $bootstrap_samples
