clear all
clc

spm('defaults','fmri');
spm_jobman('initcfg');


% Modify the below lines


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


base_dir = '/data/';    % Your project folder
project_dir = strcat(base_dir, 'Project1/MSHP/');

spm_dir = '/usr/local/MATLAB/R2023b/toolbox/spm12/';
addpath(spm_dir)

init = 1;
total_sub_num = 53;    % Number of subjects
total_sub_number = [init:1:init+total_sub_num-1];

exc_sub_num = 11;
exc_sub_number = [8, 11, 13, 20, 29, 30, 34, 38, 39, 47, 50];

%exc_sub_num = 14;
%exc_sub_number = [2, 6, 8, 11, 13, 17, 20, 29, 35, 38, 39, 41, 47, 50];

sub_num = total_sub_num - exc_sub_num;
sub_number = setdiff(total_sub_number, exc_sub_number);

subjects = {};

for i = 1:length(sub_number)
    subjects{i} = sprintf('sub%02d', sub_number(i));
end

% name of contrast
contrast_name = {'R','D'};

% second level contrast results folder
contrast_folder = contrast_name; 

mkdir(sprintf('%sanalysis/mvpa_2_results_2nd',project_dir))
for c = 1:length(contrast_folder)
    mkdir(sprintf('%sanalysis/mvpa_2_results_2nd/searchlight/%s/between_alternative/PFC',project_dir,contrast_folder{c}))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Second level model specification
for i = 1:length(contrast_name)
    for s = 1:length(subjects)
        second_contrasts{i}{s} = strcat(project_dir, 'analysis/', subjects{s}, '/analysis/multivariate_analysis_2/', contrast_folder{i}, '/between_alternative/mswres_accuracy_minus_chance.nii,1');
    end
end

contrast_dir = {};

% Estimation directory
Est_file = {};

for i=1:length(contrast_folder)
    contrast_dir{i} = strcat(project_dir, 'analysis/mvpa_2_results_2nd/searchlight/', contrast_folder{i}, '/between_alternative/PFC');
    Est_file{i} = strcat(project_dir, 'analysis/mvpa_2_results_2nd/searchlight/', contrast_folder{i}, '/between_alternative/PFC/SPM.mat');
end

for i=1:length(contrast_folder)
    matlabbatch = {};
    
    matlabbatch{1}.spm.stats.factorial_design.dir = {contrast_dir{i}};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = second_contrasts{i}';
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/data/Project1/MSHP/mask/left_PFC.nii'};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    matlabbatch{1}.spm.stats.factorial_design.delete = 1;
    
    
    
    %% Second level model estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {Est_file{i}};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    matlabbatch{2}.spm.stats.fmri_est.delete = 1;
    
    
    
    fprintf('2nd-level model specification \n')
    fprintf('Contrast: %d/%d \n', i, length(contrast_name))
    
    
    %% Second level contrast manager
    matlabbatch{3}.spm.stats.con.spmmat = {Est_file{i}};
    
    
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = contrast_name{i};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.delete = 1;
    
    
    spm_jobman('run',matlabbatch);
    
    fprintf('2nd-level Contrast Manager \n')
    fprintf('Contrast: %d/%d \n', c, length(contrast_folder))

end
