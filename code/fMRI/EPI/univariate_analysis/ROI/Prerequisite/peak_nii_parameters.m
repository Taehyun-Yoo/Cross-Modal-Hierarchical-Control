clear all
clc

result_dir = '/data/Project1/MSHP/analysis/results_2nd/parametric_analysis/D_parametric_PFC';

mkdir(sprintf('%s/peak_nii',result_dir))
cd(sprintf('%s/peak_nii',result_dir))

addpath(genpath('/data/Project1/MSHP/Code/EPI/univariate_analysis/ROI/Prerequisite/MRtools_Hoffman2-master/'))

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
mapparameters.mask='/data/Project1/MSHP/left_PFC.nii';
mapparameters.conn=18;
mapparameters.nearest=1;

peak_nii(strcat(result_dir,'/increasing_complexity/T_cluster_1.nii'), mapparameters)