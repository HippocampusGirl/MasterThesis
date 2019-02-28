#!/bin/bash

CWD=$(pwd)
cd ~/Dropbox/Univ/masterThesis/data/Ephys

ulimit -n 4096

# autaptomatic -regions -measure AMPLITUDE CHARGE RISE_TIME N_PEAKS -baseline 0 0.1 -average $(cat PP_sort.txt) 
# mv measurements.xlsx ${CWD}/measurements.xlsx

# autaptomatic -regions -measure AMPLITUDE CHARGE RISE_TIME N_PEAKS -baseline 0 0.1 -average $(cat DG_PP_sort.txt) 
# mv measurements.xlsx ${CWD}/hippocampus.xlsx

# autaptomatic -regions -measure AMPLITUDE CHARGE -baseline 0 0.1 $(cat Sucrose_sort.txt) 
# mv measurements.xlsx ${CWD}/sucrose_1ap.xlsx

# autaptomatic -sucrose -measure CHARGE -sucrose $(cat Sucrose_sort.txt)  
# mv measurements.xlsx ${CWD}/sucrose.xlsx

cd ${CWD}