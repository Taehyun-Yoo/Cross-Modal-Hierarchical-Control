clear all
clc

base_dir = '/data/';    % Your project folder
project_dir = strcat(base_dir, 'Project1/MSHP/');

spm_dir = '/usr/local/MATLAB/R2023b/toolbox/spm12/';
addpath(spm_dir)

MNI_brain_file = strcat(project_dir, 'MNI152_T1_1mm_brain.nii'); % MNI brain image file

init = 1;
total_sub_num = 53;    % Number of subjects
total_sub_number = [init:1:init+total_sub_num-1];

exc_sub_num = 7;
exc_sub_number = [13, 20, 29, 38, 39, 47, 50];

sub_num = total_sub_num - exc_sub_num;
sub_number = setdiff(total_sub_number, exc_sub_number);

run_num = 8;    % Number of Runs
run_num_1 = 2;
run_num_2 = 4;
run_num_3 = 2;

runs_num = [run_num_1 run_num_2 run_num_3];

hier_num = 3;
hier = {'R', 'F', 'D'};

voxel_size = [3 3 3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subjects = {};

for i = 1:length(sub_number)
    subjects{i} = sprintf('sub%02d', sub_number(i));
end

% Modify the below lines
% EPI_1, EPI_2, EPI_3, EPI_4, EPI_5, EPI_6, EPI_7, EPI_8, EPI_9, EPI_10, EPI_11, EPI_12, MPRAGE - each subject
% Label of MPRAGE.IMA in the LAST!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dicom_order = {'Response_1' 'Response_2' 'Feature_1' 'Feature_2' 'Feature_3' 'Feature_4' 'Dimension_1' 'Dimension_2' 'T1'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% DICOM Import
for s = 2:length(subjects)
    dicom_dir_day1 = strcat(project_dir, 'RAW/main/fMRI_data/', subjects{s}, '/Day1/Data');
    dicom_dir_day2 = strcat(project_dir, 'RAW/main/fMRI_data/', subjects{s}, '/Day2/Data');
    
    for f = 1:length(dicom_order)
       
        cd (dicom_dir_day2);
     
        if exist(dicom_order{f}, 'dir') == 7
            Day = '/Day2';
        else
            Day = '/Day1';
        end
        
        data{f} = fullfile(strcat(base_dir, subjects{s}, Day, '/Raw/Data/', dicom_order{f}, '/*.IMA'));
        files{f} = dir(data{f});
    end
    
    
    % Modify the below lines
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    run_1 = {};
    run_2 = {};
    run_3 = {};
    run_4 = {};
    run_5 = {};
    run_6 = {};
    run_7 = {};
    run_8 = {};
    anat = {};
    runs = {run_1, run_2, run_3, run_4, run_5, run_6, run_7, run_8, anat};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for f = 1:length(dicom_order)
        
        for j=1:length(files{f})
            file = files{f}(j);
            input_dir = strcat(file.folder, '/');
        
            runs{f}{length(runs{f})+1} = strcat(input_dir, file.name);
            
        end
        
    end
        
    % Modify the below lines
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for f = 1:length(dicom_order)
        runs{f} = runs{f}';
    end
    
    runfolders = {'run_1', 'run_2', 'run_3', 'run_4', 'run_5', 'run_6', 'run_7', 'run_8'}; 
    
    for h = 1:hier_num
        mkdir(sprintf('%sanalysis/%s/functional_concat/%s',project_dir,subjects{s},hier{h}))
        mkdir(sprintf('%sanalysis/%s/functional_concat/%s/all',project_dir,subjects{s},hier{h}))       
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    con_new_input = {};
    
    for h = 1:hier_num
        con_new_input{h} = strcat(project_dir, 'analysis/', subjects{s},'/functional_concat/',hier{h},'/all/');
    end
    
    for r = 1:run_num

        old_data{r} = fullfile(strcat(project_dir, 'analysis/', subjects{s},'/functional/run_', num2str(r), '/f*.nii'));
        old_files{r} = dir(old_data{r});

        if sum(ismember(1:sum(runs_num(1)),r))
            for j = 1:length(old_files{r})
                old_file = old_files{r}(j);
                cd(old_file.folder);
                copyfile(sprintf(old_file.name), con_new_input{1});
            end

        elseif sum(ismember(sum(runs_num(1))+1:sum(runs_num(1:2)),r))
            for j = 1:length(old_files{r})
                old_file = old_files{r}(j);
                cd(old_file.folder);
                copyfile(sprintf(old_file.name), con_new_input{2});
            end

        elseif sum(ismember(sum(runs_num(1:2))+1:sum(runs_num(1:3)),r))
            for j = 1:length(old_files{r})
                old_file = old_files{r}(j);
                cd(old_file.folder);
                copyfile(sprintf(old_file.name), con_new_input{3});
            end
        end
    end
    
    %% Variable specification
    
    input = {};
    rinput = {};
    mrinput = {};
    
    
    data = {};
    files = {};
    for h = 1:hier_num
        data{h} = fullfile(strcat(project_dir, 'analysis/', subjects{s},'/functional_concat/', hier{h}, '/all/f*.nii'));
        files{h} = dir(data{h});
    
        for j=1:length(files{h})
            file = files{h}(j);
            input_dir = strcat(file.folder, '/');
        
            input{h}{j} = strcat(input_dir, file.name);
            rinput{h}{j} = strcat(input_dir, 'r', file.name);
        end

        meaninput{h} = strcat(files{h}(1).folder, '/mean', files{h}(1).name, ',1');
    end
    
    
    %% Realignment

    matlabbatch = {};
    for h = 1:hier_num

        matlabbatch{h}.spm.spatial.realign.estwrite.data = {input{h}'}';
        %%
        matlabbatch{h}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{h}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{h}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{h}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch{h}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{h}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{h}.spm.spatial.realign.estwrite.eoptions.weight = '';

        matlabbatch{h}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{h}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{h}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{h}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{h}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    end
    spm_jobman('run',matlabbatch);

    %spm_jobman('run',matlabbatch);

    fprintf('Realignment \n')
    fprintf('Subject: %d/%d \n', s, sub_num)
    
    
    %% Coregistration

    matlabbatch = {};
    for h = 1:hier_num
        matlabbatch{h}.spm.spatial.coreg.estimate.ref = {strcat(project_dir, 'analysis/', subjects{s},'/anatomical/ana.nii,1')};
        matlabbatch{h}.spm.spatial.coreg.estimate.source = {meaninput{h}};
        matlabbatch{h}.spm.spatial.coreg.estimate.other = rinput{h}';
        matlabbatch{h}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{h}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{h}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{h}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    end
    spm_jobman('run',matlabbatch);

    fprintf('Coregistration \n')
    fprintf('Subject: %d/%d \n', s, sub_num)

end