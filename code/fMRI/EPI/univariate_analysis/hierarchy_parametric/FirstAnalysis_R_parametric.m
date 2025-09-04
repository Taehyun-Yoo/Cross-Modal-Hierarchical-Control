clear all
clc

base_dir = '/data/';    % Your project folder
project_dir = strcat(base_dir, 'Project1/MSHP/');

spm_dir = '/usr/local/MATLAB/R2023b/toolbox/spm12/';
addpath(spm_dir)

init = 1;
total_sub_num = 53;    % Number of subjects
total_sub_number = [init:1:init+total_sub_num-1];

exc_sub_num = 11;
exc_sub_number = [8, 11, 13, 20, 29, 30, 34, 38, 39, 47, 50];

sub_num = total_sub_num - exc_sub_num;
sub_number = setdiff(total_sub_number, exc_sub_number);

run_num = 2;        % Number of runs

runs = [run_num];

units = 'secs';     % scans or secs
TR = 2;             % Repetition time

%slice_num = 42;     % Number of slices per volume
highpassfilter = 128;       % high pass filter

% In case of slice timing corrected data
% you should modified microtime_resulution and microtime_onset.
microtime_resolution = 16;  % number of slices
microtime_onset = 8;       % half of number of slices

R_duration = 23;    % time window

subjects = {};

for i = 1:length(sub_number)
    subjects{i} = sprintf('sub%02d', sub_number(i));
end

% onset

onsets_para = {};

for s = 1:length(subjects) 
    for r = 1:run_num
        fid_para = fopen(strcat(project_dir, 'RAW/main/behavioral_data/Onset/', subjects{s}, '/', num2str(r), '_task.txt'), 'rt');
        onsets_para{s}{r}{1} = cell2mat(textscan(fid_para, '%f', 'Delimiter', '\n') );
        fclose(fid_para);
    end
end
    
contrast_name = {'increasing complexity'};

cont_w = [zeros(1,1), 1, zeros(1,6), zeros(1,1), 1, zeros(1,6), zeros(1,2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% First level model specification
for s = 1:length(subjects) 
    % subject specific directory
    subj_dir = fullfile(project_dir, 'analysis', sprintf('sub%02d', sub_number(s)));
    
    % initialize matlabbatch for each subject
    
    matlabbatch = {};
    
    
    % Modify the below lines
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    %if isfile(strcat(base_dir, subjects{s}, '/analysis/SPM.mat'))
    %    delete(strcat(base_dir, subjects{s}, '/analysis/*'))
    %end
    
    % For R2017a and previous releases, use the "exist" function.
    if exist(strcat(subj_dir, '/analysis/parametric_analysis/R_parametric/SPM.mat'), 'file') == 2
       delete(strcat(subj_dir, '/analysis/parametric_analysis/R_parametric/*'))
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
        
    matlabbatch{1}.spm.stats.fmri_spec.dir = {fullfile(subj_dir, 'analysis/parametric_analysis/R_parametric')};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = units;
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = microtime_resolution;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = microtime_onset;
    
    for r = 1:run_num
        
        run_dir = fullfile(subj_dir, 'functional', sprintf('run_%d', r));
        
        % input image
        subj_scans =  spm_select('FPList', run_dir, '^ms.*');
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cellstr(subj_scans);
            
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.name = 'task';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.onset = onsets_para{s}{r}{1};
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.duration = R_duration;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.pmod.name = 'complexity';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.pmod.param = [1; 1; 1; 1; 2; 2; 2; 2; 4; 4; 4; 4];
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.pmod.poly = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond.orth = 1;               
        
        % other settings for each session
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {spm_select('FPList', run_dir, '^rp.*')};
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = highpassfilter;
        
    end % !session
    
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    fprintf('1st-level model specification \n')
    fprintf('Subject: %d/%d \n', s, sub_num)
    
    %% First level model estimation
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(subj_dir, 'analysis/parametric_analysis/R_parametric', 'SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    
    fprintf('1st-level model estimation \n')
    fprintf('Subject: %d/%d \n', s, sub_num)
    
    
    
    %% First level contrast manager
    
    matlabbatch{3}.spm.stats.con.spmmat = {fullfile(subj_dir, 'analysis/parametric_analysis/R_parametric', 'SPM.mat')};
    matlabbatch{3}.spm.stats.con.delete = 1;
    
    for c = 1:length(contrast_name)
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = contrast_name{1};
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = cont_w;
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
    end
    
    spm_jobman('run',matlabbatch);   
    
    fprintf('1st-level Contrast Manager \n')
    fprintf('Subject: %d/%d \n', s, sub_num)
end