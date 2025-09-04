clear all
clc

%% Setup
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

hier_num = 2;
hier = {'R', 'D'};

for i = 1:length(sub_number)
    subjects{i} = sprintf('sub%02d', sub_number(i));
end

for s = 1:length(subjects)
    subj_dir = fullfile(project_dir, 'analysis', sprintf('sub%02d', sub_number(s)));
    % set path to native T1 nifti-file
    Subnative = fullfile(subj_dir, 'anatomical', 'ana.nii');
    % set path to y_*.nii file defining transformation into MNI-space
    Suby = fullfile(subj_dir, 'anatomical', 'y_ana.nii');
    
    % set output folder of textfile with coordinates
    for h = 1: hier_num
        beta_dir = fullfile(subj_dir, 'analysis', 'multivariate_analysis_2', hier{h});
        if exist(sprintf('%s/mask',beta_dir), 'dir') ~= 7
            mkdir(sprintf('%s/mask',beta_dir));
        end
    end

    %% Define MNI coordinates
    % coordinates from univariate analaysis
    mni_coords = [-45 -16 53
        -36 26 41]; %the mni coordinate of the brain region that you want

    mni_coords_labels = {'mvpa_R_ROI', 'mvpa_D_ROI'};
    
    % coordinates from multivariate analaysis
    % mni_coords = [-27 -1 56
    %     -30 29 44
    %     -30 53 8]; %the mni coordinate of the brain region that you want
    % 
    % mni_coords_labels = {'mvpa_R_ROI', 'mvpa_F_ROI', 'mvpa_D_ROI'};
    
    %% get native coordinates for each row of MNI-coordinates
    for iRow = 1:size(mni_coords, 1)
        native_coords(iRow,:) = calculate_subj_coord(Subnative, Suby, mni_coords(iRow,:));
    end
    
    %% create matrix of original MNI coordinates and coordinates in subject space with region labels
    % clear output_matrix
    % output_matrix(:,1) = mni_coords_labels;
    % output_matrix(:,2:4) = num2cell(mni_coords);
    % output_matrix(:,5:7) = num2cell(round(native_coords,2));
    % % convert output-matrix to table
    % output_table = cell2table(output_matrix,'VariableNames',...
    %     {'Region','MNI_x','MNI_y','MNI_z','Subj_x','Subj_y','Subj_z'});
    
    for h = 1:hier_num
        beta_dir = fullfile(subj_dir, 'analysis', 'multivariate_analysis_2', hier{h});
        list = dir(strcat(beta_dir, '/trial_RT_block_original/*_', hier{h}, '*_*/beta_0001.nii'));
        mask_img = fullfile(list(1).folder, 'mask.nii');
        output_folder = fullfile(beta_dir, 'mask');
        for iRow = 1:size(mni_coords, 1)
            [~,tmp] = system(sprintf('echo %s %s %s | std2imgcoord -img %s -std %s -vox', native_coords(iRow,1), native_coords(iRow,2), native_coords(iRow,3), mask_img, mask_img));
            tmp = strsplit(tmp, ' ');
            vox_loc_x = round(str2double(tmp{1}));
            vox_loc_y = round(str2double(tmp{2}));
            vox_loc_z = round(str2double(tmp{3}));
            
            [~,cmd] = system(sprintf('fslmaths %s -mul 0 -add 1 -roi %s 1 %s 1 %s 1 0 1 %s/%s_point.nii.gz -odt float', mask_img, string(vox_loc_x), string(vox_loc_y), string(vox_loc_z), output_folder, mni_coords_labels{iRow}));
            [~,cmd] = system(sprintf('fslmaths %s/%s_point.nii.gz -kernel sphere 9 -fmean -bin %s/final_%s.nii.gz', output_folder, mni_coords_labels{iRow}, output_folder, mni_coords_labels{iRow}));
            [~,cmd] = system(sprintf('gunzip %s/final_%s.nii.gz', output_folder, mni_coords_labels{iRow}));
        end
    end
    %% save to textfile
    %writetable(output_table,[output_folder '/native_coords.txt'],'Delimiter','\t');
end
