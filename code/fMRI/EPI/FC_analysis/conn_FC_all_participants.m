clear all
clc
    
base_dir = '/data/';    % Your project folder
project_dir = strcat(base_dir, 'Project1/MSHP/');

roi_dir = strcat(project_dir, 'Code/EPI/FC_analysis/ROI/Prerequisite/anatomical_mask_3mm/');
sensory_roi_names = {'left_STG', 'left_primary_auditory_cortex', 'left_IPS', 'left_primary_visual_cortex', ...
    'right_STG', 'right_primary_auditory_cortex', 'right_IPS', 'right_primary_visual_cortex'};

result_dir = strcat(project_dir, 'analysis/conn_2nd_left_right/');

spm_dir = '~/spm12';
addpath(spm_dir)
tool_dir = strcat(spm_dir, '/toolbox/conn');
addpath(tool_dir)

tot_hier = {'R', 'F', 'D'};
tot_hier_num = 3;

tot_complexity = {'1', '2', '4'};
tot_complexity_num = 3;

tot_run_num = {};
tot_run_num{1} = 2;
tot_run_num{2} = 4;
tot_run_num{3} = 2;

tot_block_duration = {};
tot_block_duration{1} = 23;
tot_block_duration{2} = 39;
tot_block_duration{3} = 39;

tot_condition_name = {};
for h = 1:tot_hier_num
    tot_condition_name{h} = {strcat(tot_hier{h}, tot_complexity{1}), strcat(tot_hier{h}, tot_complexity{2}), strcat(tot_hier{h}, tot_complexity{3})};
end
tot_condition_num = 3;

TR=2;

init = 1;
total_sub_num = 53;    % Number of subjects
total_sub_number = [init:1:init+total_sub_num-1];

exc_sub_num = 11;
exc_sub_number = [8, 11, 13, 20, 29, 30, 34, 38, 39, 47, 50];

sub_num = total_sub_num - exc_sub_num;
sub_number = setdiff(total_sub_number, exc_sub_number);

NSUBJECTS = sub_num;

subjects = {};

for s = 1:length(sub_number)
    subjects{s} = sprintf('sub%02d', sub_number(s));
end

if exist(strcat(result_dir), 'dir') ~= 7
    mkdir(strcat(result_dir))
end

for h = 1:tot_hier_num
    
    if h == 1
        functional_roi_names = {'R_2voxels_roi'};
    elseif h == 2
        functional_roi_names = {'F_2voxels_roi'};
    else
        functional_roi_names = {'D_2voxels_roi'};
    end
    
    roi_names = [functional_roi_names, sensory_roi_names];

    run_num = tot_run_num{h};
    condition_name = tot_condition_name{h};

    if exist(strcat(result_dir, tot_hier{h}), 'dir') ~= 7
        mkdir(strcat(result_dir, tot_hier{h}))
    end
    
    fun_data = {};
    fun_files = {};
    functional = {};
    
    onsets = {};

    Realignment_file = {};

    for s = 1:length(sub_number)
        subj_dir = fullfile(project_dir, 'analysis', sprintf('sub%02d', sub_number(s)));
        
        anatomical{s} = fullfile(strcat(subj_dir, '/anatomical/wana.nii'));

        for r = 1:run_num
            if h == 1
                fun_data{s}{r} = fullfile(strcat(subj_dir,'/functional/run_', num2str(r), '/mwr*.nii'));
            elseif h == 2
                fun_data{s}{r} = fullfile(strcat(subj_dir,'/functional/run_', num2str(r+2), '/mwr*.nii'));
            else
                fun_data{s}{r} = fullfile(strcat(subj_dir,'/functional/run_', num2str(r+6), '/mwr*.nii'));
            end
            fun_files{s}{r} = dir(fun_data{s}{r});
            for j=1:length(fun_files{s}{r})
                fun_file = fun_files{s}{r}(j);
                input_dir = strcat(fun_file.folder, '/');      
                functional{s}{r}{j} = strcat(input_dir, fun_file.name);
            end
        end
        
        for r = 1:run_num
            if h == 1
                for c = 1:tot_condition_num
            % load txt file
                    fid = fopen(strcat(project_dir, 'RAW/main/behavioral_data/Onset/', subjects{s}, '/', num2str(r), '_', condition_name{c}, '.txt'), 'rt');
                    onsets{r}{s}{c} = cell2mat(textscan(fid, '%f', 'Delimiter', '\n') );
                    fclose(fid);
                end
                
            elseif h == 2
                for c = 1:tot_condition_num
            % load txt file
                    fid = fopen(strcat(project_dir, 'RAW/main/behavioral_data/Onset/', subjects{s}, '/', num2str(r+2), '_', condition_name{c}, '.txt'), 'rt');
                    onsets{r}{s}{c} = cell2mat(textscan(fid, '%f', 'Delimiter', '\n') );
                    fclose(fid);
                end
                
            else
                for c = 1:tot_condition_num
            % load txt file
                    fid = fopen(strcat(project_dir, 'RAW/main/behavioral_data/Onset/', subjects{s}, '/', num2str(r+6), '_', condition_name{c}, '.txt'), 'rt');
                    onsets{r}{s}{c} = cell2mat(textscan(fid, '%f', 'Delimiter', '\n') );
                    fclose(fid);
                end 
                
            end
        end

        for r = 1:run_num
            if h == 1
                run_dir = fullfile(subj_dir, 'functional', sprintf('run_%d', r));
            elseif h == 2
                run_dir = fullfile(subj_dir, 'functional', sprintf('run_%d', r+2));
            else
                run_dir = fullfile(subj_dir, 'functional', sprintf('run_%d', r+6));
            end
            Realignment_file(s,r) = {spm_select('FPList', run_dir, '^rp.*')};
        end
    end
    
    clear batch;
    
    batch = {};
    batch.filename = fullfile(strcat(result_dir, tot_hier{h}, '/conn_denoising.mat'));
    
    batch.Setup.isnew = 1;
    batch.Setup.nsubjects = NSUBJECTS;
    batch.Setup.RT = TR;
    batch.Setup.acquisitiontype = 1;

    batch.Setup.conditions.names = condition_name;
    
    for c = 1:tot_condition_num
        for s = 1:length(sub_number)
            for r = 1:run_num 
                batch.Setup.conditions.durations{c}{s}{r} = tot_block_duration{h};
                batch.Setup.conditions.onsets{c}{s}{r} = onsets{r}{s}{c};
            end
        end
    end
    
    batch.Setup.rois.names = roi_names;
    
    for v = 1:length(batch.Setup.rois.names)
        batch.Setup.rois.files{v}=fullfile(strcat(roi_dir, roi_names{v}, '.nii,1'));
    end

    batch.Setup.functionals=repmat({{}},[NSUBJECTS,1]);     % Point to functional volumes for each subject/session
    
    for s = 1:length(sub_number)
        for r = 1:run_num
            batch.Setup.functionals{s}{r} = functional{s}{r};
        end
    end
    
    for s = 1:length(sub_number)
        batch.Setup.structurals{s} = anatomical{s};
    end
    
    batch.Setup.voxelresolution= 3;
    batch.Setup.analysisunits = 1;
    batch.Setup.outputfiles = [0,1,0];
    
    batch.Setup.covariates.names = {'realignment'};
    batch.Setup.covariates.files{1} = repmat({{}},[NSUBJECTS,1]);

    for s = 1:length(sub_number)
        for r=1:run_num
            batch.Setup.covariates.files{1}{s}{r} = Realignment_file(s,r);
        end
    end
    
    batch.Setup.analyses = [1,2];                             % seed-to-voxel and ROI-to-ROI pipelines
    batch.Setup.done = 1;
    batch.Setup.overwrite = 'Yes';
    
    batch.Denoising.filter = [0.008, inf];                 % frequency filter (band-pass values, in Hz)
    batch.Denoising.detrending = 1;
    batch.Denoising.regbp = 1;
    batch.Denoising.confounds.names = {'White Matter', 'CSF', 'realignment', sprintf('Effect of %s', condition_name{1}), sprintf('Effect of %s', condition_name{2}), sprintf('Effect of %s', condition_name{3})};
    batch.Denoising.confounds.dimensions = {5, 5, 12, inf, inf, inf};
    batch.Denoising.counfound.deriv = {0, 0, 1, 1, 1, 1};
    batch.Denoising.done = 1;
    batch.Denoising.overwrite='Yes';
    
    batch.Analysis.type = 1;
    batch.Analysis.done = 1;
    batch.Analysis.overwrite = 1;

    conn_batch(batch);

end