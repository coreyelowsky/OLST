#!/bin/bash

echo ""
echo "####################################"
echo "Logistic Regression - Bootstrapping"
echo "####################################"
echo ""

# path for bootsrapping
# all data should be in this path
export data_path=/grid/osten/data_norepl/elowsky/OLST/logistic_regression_cluster/barrel_cortex_layer_6/
export morphometrics_path=${data_path}normalized_morphometrics.txt
export clusters_path=${data_path}cluster_info.csv
export bootstrap_samples=10000

# get current directory
export cur_dir=`pwd`"/"

# make logs path
export logs_path=${data_path}logs/
mkdir -p $logs_path

# calculate number of features
export num_features=$((`awk '{print NF}' $morphometrics_path | sort -nu | tail -n 1` - 2))

echo "# Bootstrap Samples: ${bootstrap_samples}"
echo "# Features: ${num_features}"
echo ""

# submit jobs to cluster
export job_name=logistic_regression_bootstrap
export threads_per_job=2
export memory_per_job=2
export memory_per_thread=$((memory_per_job/threads_per_job+1)) 
export log_reg_script=/grid/osten/data_norepl/elowsky/OLST/logistic_regression_cluster/logistic_regression.sh
export log_reg_python_script=/grid/osten/data_norepl/elowsky/OLST/logistic_regression_cluster/logistic_regression.py
export merge_coeffs_bash_script=/grid/osten/data_norepl/elowsky/OLST/logistic_regression_cluster/merge_coeffs.sh
export merge_coeffs_python_script=/grid/osten/data_norepl/elowsky/OLST/logistic_regression_cluster/merge_coeffs.py

echo "Load Modules to use correct python version...."
module load EBModules
module load Python/3.8.6-GCCcore-10.2.0

export qsub_output=`qsub -N $job_name -cwd -binding linear_per_task:1 -pe threads $((threads_per_job/2)) -l m_mem_free="$((memory_per_thread*2))"G -t 1-$num_features $log_reg_script`

# parse qsub output to get job id
export job_id=`echo $qsub_output | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`

# call merge coeffs
nohup $merge_coeffs_bash_script > "${logs_path}nohup_merge_coeffs.txt" &









