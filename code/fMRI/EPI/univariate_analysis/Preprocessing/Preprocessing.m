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

exc_sub_num = 11;
exc_sub_number = [8, 11, 13, 20, 29, 30, 34, 38, 39, 47, 50];

sub_num = total_sub_num - exc_sub_num;
sub_number = setdiff(total_sub_number, exc_sub_number);

run_num = 8;    % Number of Runs

voxel_size = [3 3 3];
FWHM = [8 8 8];

subjects = {};

for i = 1:length(sub_number)
    subjects{i} = sprintf('sub%02d', sub_number(i));
end

% Modify the below lines
% EPI_1, EPI_2, EPI_3, EPI_4, EPI_5, EPI_6, EPI_7, EPI_8, MPRAGE - each subject
% Label of MPRAGE.IMA in the LAST!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dicom_order = {'Response_1' 'Response_2' 'Feature_1' 'Feature_2' 'Feature_3' 'Feature_4' 'Dimension_1' 'Dimension_2' 'T1'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DICOM Import
for s = 1:length(subjects) 
    dicom_dir_day1 = strcat(project_dir, 'RAW/main/fMRI_data/', subjects{s}, '/Day1/Data');
    dicom_dir_day2 = strcat(project_dir, 'RAW/main/fMRI_data/', subjects{s}, '/Day2/Data');
    
    for f = 1:length(dicom_order)
       
        cd (dicom_dir_day2);
     
        if exist(dicom_order{f}, 'dir') == 7
            Day = '/Day2';
        else
            Day = '/Day1';
        end
        
        matlabbatch = {};
        data{f} = fullfile(strcat(project_dir, 'RAW/main/fMRI_data/', subjects{s}, Day, '/Data/', dicom_order{f}, '/*.IMA'));
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
    
    mkdir(sprintf('%sanalysis/%s/analysis',project_dir,subjects{s}))
    mkdir(sprintf('%sanalysis/%s/anatomical',project_dir,subjects{s}))
    mkdir(sprintf('%sanalysis/%s/functional',project_dir,subjects{s}))
    
    for r = 1:run_num
        mkdir(sprintf('%sanalysis/%s/functional/%s',project_dir,subjects{s},runfolders{r}))
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for r = 1:length(runs)
        if r == length(dicom_order(1,:))
            %%
            matlabbatch{r}.spm.util.import.dicom.data = runs{r};
            %%
            matlabbatch{r}.spm.util.import.dicom.root = 'flat';
            matlabbatch{r}.spm.util.import.dicom.outdir = {strcat(project_dir, 'analysis/', subjects{s},'/anatomical')};
            matlabbatch{r}.spm.util.import.dicom.protfilter = '.*';
            matlabbatch{r}.spm.util.import.dicom.convopts.format = 'nii';
            matlabbatch{r}.spm.util.import.dicom.convopts.icedims = 0;
            %%
        else
            %%
            matlabbatch{r}.spm.util.import.dicom.data = runs{r};
            %%
            matlabbatch{r}.spm.util.import.dicom.root = 'flat';
            matlabbatch{r}.spm.util.import.dicom.outdir = {strcat(project_dir, 'analysis/', subjects{s},'/functional/run_', num2str(r))};
            matlabbatch{r}.spm.util.import.dicom.protfilter = '.*';
            matlabbatch{r}.spm.util.import.dicom.convopts.format = 'nii';
            matlabbatch{r}.spm.util.import.dicom.convopts.icedims = 0;
            %%
        end
    end
    
    spm_jobman('run',matlabbatch);
    
    
    %fprintf('DICOM IMPORT \n')
    %fprintf('Subject: %d/%d \n', s, sub_num)
    
    
    %% Rename
    
    data = fullfile(strcat(project_dir, 'analysis/', subjects{s}, '/anatomical/*.nii'));
    anafile = dir(data); 
    
    ana_oldname = strcat(anafile.folder, '/', anafile.name);
    
    movefile(ana_oldname, strcat(anafile.folder, '/ana.nii'));
    
    
    %% Variable specification
    
    input = {};
    rinput = {};
    wrinput = {};
    swrinput = {};
    mswrinput = {};
    
    
    data = {};
    files = {};
    for r = 1:run_num
        
        data{r} = fullfile(strcat(project_dir, 'analysis/', subjects{s},'/functional/run_', num2str(r), '/f*.nii'));
        files{r} = dir(data{r});
        
        for j=1:length(files{r})
            file = files{r}(j);
            input_dir = strcat(file.folder, '/');
            
            input{r}{j} = strcat(input_dir, file.name);
            rinput{r}{j} = strcat(input_dir, 'r', file.name);
            wrinput{r}{j} = strcat(input_dir, 'wr', file.name);
            swrinput{r}{j} = strcat(input_dir, 'swr', file.name);
            mswrinput{r}{j} = strcat(input_dir, 'mswr', file.name);    
        end
    
        meaninput{r} = strcat(files{r}(1).folder, '/mean', files{r}(1).name, ',1');
    end
    
    
    %% Realignment
    
    matlabbatch = {};
    for r = 1:run_num
        
        matlabbatch{r}.spm.spatial.realign.estwrite.data = {input{r}'}';
        %%
        matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{r}.spm.spatial.realign.estwrite.eoptions.weight = '';
    
        matlabbatch{r}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{r}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{r}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{r}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{r}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    end
    
    spm_jobman('run',matlabbatch);
    
    %fprintf('Realignment \n')
    %fprintf('Subject: %d/%d \n', s, sub_num)
    
    
    %% Coregistration
    
    matlabbatch = {};
    for r = 1:run_num
        matlabbatch{r}.spm.spatial.coreg.estimate.ref = {strcat(project_dir, 'analysis/', subjects{s},'/anatomical/ana.nii,1')};
        matlabbatch{r}.spm.spatial.coreg.estimate.source = {meaninput{r}};
        matlabbatch{r}.spm.spatial.coreg.estimate.other = rinput{r}';
        matlabbatch{r}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{r}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{r}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{r}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    end
    
    spm_jobman('run',matlabbatch);
    
    %fprintf('Coregistration \n')
    %fprintf('Subject: %d/%d \n', s, sub_num)
    
    
    %% Segmentation
    
    matlabbatch = {};
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {strcat(project_dir, 'analysis/', subjects{s},'/anatomical/ana.nii,1')};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {strcat(spm_dir, 'tpm/TPM.nii,1')};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {strcat(spm_dir, 'tpm/TPM.nii,2')};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {strcat(spm_dir,  'tpm/TPM.nii,3')};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {strcat(spm_dir, 'tpm/TPM.nii,4')};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {strcat(spm_dir, 'tpm/TPM.nii,5')};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {strcat(spm_dir, 'tpm/TPM.nii,6')};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];  % [inverse forward] --- deformation field
    
    
    spm_jobman('run',matlabbatch);
    
    %fprintf('Segmentation \n')
    %fprintf('Subject: %d/%d \n', s, sub_num)
    
    
    %% Normalisation
    
    matlabbatch = {};
    for r = 1:run_num
        %%
        matlabbatch{r}.spm.spatial.normalise.write.subj.def = {strcat(project_dir, 'analysis/', subjects{s},'/anatomical/y_ana.nii')};
        matlabbatch{r}.spm.spatial.normalise.write.subj.resample = rinput{r}';
        %%
        matlabbatch{r}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                                   78 76 85];
        matlabbatch{r}.spm.spatial.normalise.write.woptions.vox = voxel_size;
        matlabbatch{r}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{r}.spm.spatial.normalise.write.woptions.prefix = 'w';
    end
    
    spm_jobman('run',matlabbatch);
    
    %fprintf('Normalisation \n')
    %fprintf('Subject: %d/%d \n', s, sub_num)
    
    
    %% Smoothing
    
    matlabbatch = {};
    for r = 1:run_num
        matlabbatch{r}.spm.spatial.smooth.data = wrinput{r}';
        matlabbatch{r}.spm.spatial.smooth.fwhm = FWHM;
        matlabbatch{r}.spm.spatial.smooth.dtype = 0;
        matlabbatch{r}.spm.spatial.smooth.im = 0;
        matlabbatch{r}.spm.spatial.smooth.prefix = 's';
    end
    
    spm_jobman('run',matlabbatch);
    
    %fprintf('Smoothing \n')
    %fprintf('Subject: %d/%d \n', s, sub_num)

    %% Image Calculation
    
    for r = 1:run_num
        matlabbatch = {};
        for j=1:length(files{r})
            file = files{r}(j);
            matlabbatch{j}.spm.util.imcalc.input = {
                                                    swrinput{r}{j}
                                                    MNI_brain_file
                                                    };
            matlabbatch{j}.spm.util.imcalc.output = mswrinput{r}{j};
            matlabbatch{j}.spm.util.imcalc.outdir = {file.folder};
            matlabbatch{j}.spm.util.imcalc.expression = 'i1.*((i2)>0.1)';
            matlabbatch{j}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{j}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{j}.spm.util.imcalc.options.mask = 0;
            matlabbatch{j}.spm.util.imcalc.options.interp = 1;
            matlabbatch{j}.spm.util.imcalc.options.dtype = 4;
        end
        
        
        spm_jobman('run',matlabbatch);
    
        %fprintf('Image Calculation \n')
        %fprintf('Subject: %d/%d \n', s, sub_num)
        %fprintf('Run: %d/%d \n', r, run_num)
    
    end

end