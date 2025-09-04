#!/bin/bash


# Modify the below lines
##########################################################

scripts_dir='/data/Project1/MSHP/Code/EPI/multivariate_analysis/searchlight'
cd ${scripts_dir}

##########################################################

sub_number=($(seq 1 42))

run_matlab() {
    s=${1}
    /usr/local/MATLAB/R2023b/bin/matlab -nodisplay -nosplash -nodesktop -r "Alternative_decoding_accuracy_batch_RD(${s});exit;"
}

export -f run_matlab

# use gnu paralllel to execute function
njobs=10
parallel -j $njobs run_matlab ::: ${sub_number[@]}
