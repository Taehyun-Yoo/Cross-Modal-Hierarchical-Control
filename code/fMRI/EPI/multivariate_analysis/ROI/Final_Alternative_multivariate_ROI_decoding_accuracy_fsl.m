clear all
clc

base_dir = '/data/';    % Your project folder
project_dir = strcat(base_dir, 'Project1/MSHP/');

spm_dir = '/usr/local/MATLAB/R2023b/toolbox/spm12/';
addpath(spm_dir)

setenv('FSLDIR','~/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

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

hier = {'R', 'D'};
hier_num = 2;

roi_dir = strcat(project_dir, '/Code/EPI/multivariate_analysis/ROI/Prerequisite/');
roi_mask = dir(strcat(roi_dir, '*.nii'));

for s = 1:length(subjects) 
    subj_dir = fullfile(project_dir, 'analysis', sprintf('sub%02d', sub_number(s)));
    for h = 1: hier_num
        beta_dir = fullfile(subj_dir, 'analysis', 'multivariate_analysis_2', hier{h});
        output_dir = fullfile(beta_dir, 'between_alternative_final', 'ROI');

        fi = fopen(fullfile(output_dir, sprintf('fsl_%s_decoding_accuracy_sub%02d.txt', hier{h}, sub_number(s))), 'w');
        for r = 1:length(roi_mask)
            if r == 1
                fprintf(fi, strcat('ROI accuracy','\n'));
            end
            roi_tmp = strsplit(roi_mask(r).name, '.');
            
            roi_name = roi_tmp{1};

            acc_img = fullfile(beta_dir, 'between_alternative_final', 'mswres_accuracy_minus_chance.nii');
            [~,tmp] = system(sprintf('fslstats %s -k %s -M', acc_img, fullfile(roi_mask(r).folder, roi_mask(r).name)));

            fprintf(fi, strcat([roi_name ' ' tmp], '\n'));
        end
        fclose(fi);
    end
end