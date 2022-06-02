#!/bin/bash
# This options tells gridware to change to the current directory before
# executing the job (default is the home of the user)
#$-cwd
#
# This option declares the amount of memory the job will
# use. Specify this value as accurat as possible. If you
# declare too low memory usage, your job may be aborted.
# If you declare too high memory usage, your job might be
# wait for a long time until it is started on a machine
# that has the sufficient amount of free memory.
#$-l vf=2G
#
# Specify this option only for multithreaded jobs that use
# more than one cpu core. The value 1..8 denotes the number
# of requested cpu cores. Please note that multithreaded jobs
# are always calculated on a single machine - for parallel
# jobs you should use MPI instead.
# Another important hint: memory request by -l vf=... are
# multiplied by the number of requested cpu cores. Thus you
# should divide the overall memory consumption of your job by
# the number of parallel threads.
#$-pe serial 10
#$-q bigmem
#
# -- Job name ---
#$ -N cafe
#$ -S /bin/bash
#
# Please insert your mail address!
#$-M kai.wei@tum.de

###split gene with copy number smaller than 100 and larger than 100 to two input 
python /data/proj/chilense/30_genomes_outputs/kai_wei/variant_call/CNVs/metaSV/VST/cafe/python_scripts/cafetutorial_clade_and_size_filter.py -i cafe_input.txt -o filtered_cafe_input.txt

###write parameters to a file named by cafe.sh:
## we run first cafe for gene copy number smaller than 100 to estimate lammda value using -s parameter
#!cafe
load -i smaller_cafe_input.txt -t 1 -l reports/log_run1.txt -p 0.05 
tree ((LA4107:77.22659400,LA2932:77.22659400):19.23579517,(highland:89.82318742,central:89.82318742):6.63920175)
lambda -s -t ((1,1)1,(1,1)1)
report reports/report_run1

###run cafe.sh
 export PATH=/data/home/users/k.wei/.conda/envs/cnven/bin:$PATH
 cafe cafe.sh


## then run cafe for gene copy number larger than 100 using lammda value (-l 0.00206736781311) from last step (smaller gene copy number)
#!cafe
load -i larger_cafe_input.txt -t 1 -l reports/log_run2.txt -p 0.05 
tree ((LA4107:77.22659400,LA2932:77.22659400):19.23579517,(highland:89.82318742,central:89.82318742):6.63920175)
lambda -l 0.00206736781311  -t ((1,1)1,(1,1)1)
report reports/report_run2

###run cafe.sh again
 export PATH=/data/home/users/k.wei/.conda/envs/cnven/bin:$PATH
 cafe cafe.sh

### a python script can be used to summary cafe output
python2 /data/proj/chilense/30_genomes_outputs/kai_wei/variant_call/CNVs/metaSV/VST/cafe/python_scripts/cafetutorial_report_analysis.py -i reports/report_run1.cafe -o reports/summary_run1
python2 /data/proj/chilense/30_genomes_outputs/kai_wei/variant_call/CNVs/metaSV/VST/cafe/python_scripts/cafetutorial_report_analysis.py -i reports/report_run2.cafe -o reports/summary_run2
