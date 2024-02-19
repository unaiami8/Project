#!/bin/bash
#SBATCH --job-name=api-poller    # create a short name for your job
#SBATCH --partition=cpu          # use cpu partition
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=1G         # memory per cpu-core
#SBATCH --time=00:05:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email for all job steps
#SBATCH --mail-user=username@domain.com


# Here we choose which Python script we submit to the cpu cluster for 
# checking an API
/usr/bin/python3 /home/username/poller.py
