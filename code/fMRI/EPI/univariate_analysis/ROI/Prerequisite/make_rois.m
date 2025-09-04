clear all
clc

%% Use MarsBaR to make spherical ROIs

base_dir = '/data/';    % Your project folder
project_dir = strcat(base_dir, 'Project1/MSHP/');

spm_dir = '/usr/local/MATLAB/R2024b/toolbox/spm12/';
addpath(spm_dir)

addpath(strcat(spm_dir, 'toolbox/marsbar-0.45'));
marsbar('on')

%% Set general options

outDir = strcat(project_dir, 'Code/EPI/univariate_analysis/ROI/Prerequisite');
radius = 6; % mm

% coordinates are nvoxels rows by 3 columns for X,Y,Z
coords = [-27 -4 53
    -39 14 29
    -45 35 17];

% In terminal,
% fslmaths /data/Project1/MSHP/analysis/sub01/analysis/parametric_analysis/R_parametric/con_0001.nii -mul 0 -add 1 -roi 35 1 36 1 41 1 0 1 R_roi_peak.nii.gz -odt float
% fslmaths R_roi_peak.nii.gz -kernel sphere 9 -fmean R_roi_3voxels -odt float
% fslmaths R_roi_3voxels.nii.gz -bin R_roi_3voxels_mask.nii.gz

% In matlab,
% niftiFile = 'R_roi_3voxels_mask.nii'; % or .nii.gz
% niftiData = niftiread(niftiFile);
% niftiInfo = niftiinfo(niftiFile);
% matFile = 'R_roi_3voxels_mask.mat';
% save(matFile, 'niftiData', 'niftiInfo');

roi_label = {'R', 'F', 'D'};

% (alternatively, or better, you could put these in a text file and read
% them in using the dlmread function)

%% Error checking: directory exists, MarsBaR is in path
% if ~isdir(outDir)
%     mkdir(outDir);
% end

if ~exist('marsbar')
    error('MarsBaR is not installed or not in your matlab path.');
end

%% Make rois

% for i=1:size(coords,1)
%     thisCoord = coords(i,:);
% 
%     fprintf('Working on ROI %d/%d...', i, size(coords,1));
% 
%     roiLabel = char(roi_label(i));
% 
%    sphereROI = maroi_sphere(struct('centre', thisCoord, 'radius', sphereRadius));
% 
%    outName = fullfile(outDir, char(roi_label(i)));
% 
%    % save MarsBaR ROI (.mat) file
%    saveroi(sphereROI, [outName '.mat']);
% 
%    % save the Nifti (.nii) file
%    save_as_image(sphereROI, [outName '.nii']);
% 
%    fprintf('done.\n');
% 
% end

%% make .mat file and .nii file

% for pt_no = 1:size(coords, 1)
%     params = struct('centre', coords(pt_no, :), 'radius', radius);
%     roi = maroi_sphere(params);
% 
%     outName = fullfile(outDir, strcat(char(roi_label(pt_no)), '_roi_3voxels'));
%     saveroi(roi, [outName '.mat']);
% end
% 
roi_namearray = dir(fullfile(outDir, '*_ver2_2voxels_roi.mat'));

for pt_no = 1:size(coords, 1)
    roi_array{pt_no} = maroi(fullfile(outDir, roi_namearray(pt_no).name));
    roi = roi_array{pt_no};
    name = strtok(roi_namearray(pt_no).name, '.');

    save_as_image(roi, fullfile(outDir, [name '.nii']))
end

%% convert .nii file to .mat file
% marsbar -> Import -> Import ROIs from: cluster image
