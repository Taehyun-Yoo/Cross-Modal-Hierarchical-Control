clear all
clc

% MarsBaR batch script to show off MarsBaR batching
% This script replicates most of the work in the MarsBaR tutorial.
% See http://marsbar.sourceforge.net
%
% The script assumes that the current directory is the 'batch' directory
% of the example data set files.
%
% $Id: run_tutorial.m,v 1.3 2004/08/15 01:19:43 matthewbrett Exp $

% Start marsbar to make sure spm_get works
spm_dir = '/usr/local/MATLAB/R2023b/toolbox/spm12/';
addpath(spm_dir)

addpath(strcat(spm_dir, 'toolbox/marsbar-0.45'));
marsbar('on')

% You might want to define the path to the example data here, as in
% subjroot = '/my/path/somewhere';
base_dir = '/data/';    % Your project folder
project_dir = strcat(base_dir, 'Project1/MSHP/');

% Directory to store (and load) ROIs
roi_dir = strcat(project_dir, 'Code/EPI/univariate_analysis/ROI/Prerequisite/'); %  fullfile(subjroot, 'rois');

% MarsBaR version check
v = str2num(marsbar('ver'));
if v < 0.35
  error('Batch script only works for MarsBaR >= 0.35');
end

% SPM version check. We need this to guess which model directory to use and
% to get the SPM configured design name. 
spm_ver = spm('ver');
sdirname = [spm_ver '_ana'];
if strcmp(spm_ver, 'SPM99')
  conf_design_name = 'SPMcfg.mat';
else
  spm_get_defaults;
  conf_design_name = 'SPM.mat';
end

% Set up the SPM defaults, just in case
spm('defaults', 'fmri');

init = 1;
total_sub_num = 53;    % Number of subjects
total_sub_number = [init:1:init+total_sub_num-1];

exc_sub_num = 11;
exc_sub_number = [8, 11, 13, 20, 29, 30, 34, 38, 39, 47, 50];

sub_num = total_sub_num - exc_sub_num;
sub_number = setdiff(total_sub_number, exc_sub_number);

subjects = {};

for i = 1:length(sub_number)
    subjects{i} = sprintf('sub%02d', sub_number(i));
end

% load roi
rois = {'R', 'F', 'D', 'D(2)'};

roi_array = {};

for r = 1:length(rois)
    roi_array{r} = maroi(fullfile(roi_dir, sprintf('%s_2voxels_roi.mat', rois{r})));
end

% Number of runs
run_num = 4;

condition_name = {'increasing_complexity_F'};    % name of conditions

event_name = repmat(condition_name,run_num,1);

% for each subject
for s = 1:length(sub_number)
    % subject specific directory
    sub_dir = fullfile(project_dir, 'analysis', subjects{s}, 'analysis');
    
    if exist(sprintf('%s/ROI_analysis',sub_dir), 'dir') ~= 7
        mkdir(sprintf('%s/ROI_analysis',sub_dir));
    end
    
    if exist(sprintf('%s/ROI_analysis/ROI_2voxels',sub_dir), 'dir') ~= 7
        mkdir(sprintf('%s/ROI_analysis/ROI_2voxels',sub_dir));
    end

    if exist(sprintf('%s/ROI_analysis/ROI_2voxels/F_parametric',sub_dir), 'dir') ~= 7
        mkdir(sprintf('%s/ROI_analysis/ROI_2voxels/F_parametric',sub_dir));
    end

    % open file
    fi = fopen(fullfile(sub_dir, 'ROI_analysis', 'ROI_2voxels', 'F_parametric', sprintf('F_parametric_beta_sub%02d.txt', sub_number(s))), 'w');

    % first load design for each participant
    D = mardo(fullfile(sub_dir, 'parametric_analysis', 'F_parametric', 'SPM.mat'));

    % for each roi
    for j = 1:length(roi_array)

        % then load ROI
        % Extract data
        Y = get_marsy(roi_array{j}, D, 'mean');

        % then Estimate
        % MarsBaR estimation
        E = estimate(D, Y);

        % then import contrast
        gc = get_contrasts(D);

        E2 = set_contrasts(E, gc);
        
        all_beta = betas(E);
        
        for r = 1:run_num %session
            for k = 1:1 %condition per session
                % then also export percent signal change             
                if (r == 1) && (j == 1) && (k == 1)
                    fprintf(fi, strcat('Session Condition ROI beta','\n') );
                end
                fprintf(fi, strcat([num2str(r) ' ' char(event_name(r,k)) ' ' rois{j} ' ' num2str(all_beta(8*(r-1)+k+1))], '\n'));
            end
        end
    end

    fclose(fi);
end