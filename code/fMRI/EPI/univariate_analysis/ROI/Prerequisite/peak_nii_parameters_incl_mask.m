clear all
clc

result_dir = '/data/Project1/MSHP/analysis/results_2nd/parametric_analysis_2_-inf/D_parametric_whole_brain';

mkdir(sprintf('%s/peak_nii_incl_mask',result_dir))
cd(sprintf('%s/peak_nii_incl_mask',result_dir))

addpath(genpath('/data/Project1/MSHP/Code/EPI/univariate_analysis/ROI/Prerequisite_2/MRtools_Hoffman2-master/'))

mapparameters.UID='';
mapparameters.out=['D'];
mapparameters.sign='pos';
mapparameters.thresh=3.301273;
mapparameters.type='T';
mapparameters.voxlimit=20;
mapparameters.separation=12;
mapparameters.SPM=1;
mapparameters.cluster=25;
mapparameters.df1=41;
mapparameters.label='aal_MNI_V4';
%mapparameters.mask='/data/Project1/MSHP/left_PFC.nii';
mapparameters.conn=18;
mapparameters.nearest=1;

peak_nii(strcat(result_dir,'/increasing_complexity/T_cluster_3_D_unc_0_001_incl_mask_D_F_0_05.nii'), mapparameters)

rmpath(genpath('/data/Project1/MSHP/Code/EPI/univariate_analysis/ROI/Prerequisite_2/MRtools_Hoffman2-master/'))
