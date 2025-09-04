function Alternative_decoding_accuracy_block_original_batch_RD(s)

    base_dir = '/data/';    % Your project folder
    project_dir = strcat(base_dir, 'Project1/MSHP/');
    
    spm_dir = '/usr/local/MATLAB/R2023b/toolbox/spm12/';
    addpath(spm_dir)
    tool_dir = '/usr/local/MATLAB/R2023b/toolbox/tdt_3.999H/';
    addpath(genpath(tool_dir))
    tool2_dir = strcat(project_dir, 'Code/EPI/multivariate_analysis/pilot/modified_toolbox');
    addpath(genpath(tool2_dir))
    
    MNI_brain_file = strcat(project_dir, 'MNI152_T1_1mm_brain.nii'); % MNI brain image file
    
    init = 1;
    total_sub_num = 53;    % Number of subjects
    total_sub_number = [init:1:init+total_sub_num-1];
    
    exc_sub_num = 11;
    exc_sub_number = [8, 11, 13, 20, 29, 30, 34, 38, 39, 47, 50];
    
    sub_num = total_sub_num - exc_sub_num;
    sub_number = setdiff(total_sub_number, exc_sub_number);

    FWHM = [8 8 8];

    %%%% mvpa type %%%%%%%%%%%%%%
    decoding_type ='searchlight'; %%% determines decoding method ('searchlight','ROI', or 'wholebrain')
    radius= 3; %% searchlight radius in 6 mm
    
    hier = {'R', 'D'};
    hier_num = 2;

    % subject specific directory
    subj_dir = fullfile(project_dir, 'analysis', sprintf('sub%02d', sub_number(s)));

    for h = 1: hier_num
        beta_dir = fullfile(subj_dir, 'analysis', 'multivariate_analysis_2', hier{h});
        if exist(sprintf('%s/between_alternative',beta_dir), 'dir') ~= 7
            mkdir(sprintf('%s/between_alternative',beta_dir));
            mkdir(sprintf('%s/between_alternative/trial_RT_block_original',beta_dir));
            
            list = dir(strcat(beta_dir, '/trial_RT_block_original/*_*_*_*_*/beta_0001.nii'));
            for l = 1:length(list)
                tmp = strsplit(list(l).folder, '/');
                filename = tmp{11};
                list(l).filename = filename;

                copyfile(fullfile(list(l).folder, list(l).name), sprintf('%s/between_alternative/trial_RT_block_original/%s.nii',beta_dir, list(l).filename))
            end
        end

        output_dir = fullfile(beta_dir, 'between_alternative');%% please change the directory if you want to put results in another folder
        trial_dir = fullfile(output_dir, 'trial_RT_block_original');
        
        cd (fullfile(trial_dir))

        clear cfg
        
        cfg = [];
        
        cfg = decoding_defaults(cfg);

        cfg.analysis = decoding_type;

        cfg.results.dir = output_dir;

        cfg.decoding.method = 'classification_kernel'; % classification using the kernel speedup as standard
        %cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; % linear classification
        %cfg.decoding.test.classification.model_parameters = '-q'; % linear classification
        %cfg.decoding.software = 'liblinear';
        %cfg.software = spm('ver');
 
        cfg.searchlight.radius = radius;
        cfg.searchlight.unit = 'voxels';
        %cfg.searchlight.spherical = 0;

        %cfg.files.mask = fullfile(list(1).folder, 'mask.nii');
        cfg.files.mask = fullfile(list(1).folder, 'mask.nii');
        %labelname1 = xclass_names{h}{1};
        %labelname2 = xclass_names{h}{2};
        %labelname3 = xclass_names{h}{3};
        %regressor_names = design_from_spm(beta_dir);
        
        alldesign = {};
        
        list1 = dir(fullfile(trial_dir, '*_*_1_*_*.nii'));
        list2 = dir(fullfile(trial_dir, '*_*_2_*_*.nii'));
        list3 = dir(fullfile(trial_dir, '*_*_3_*_*.nii'));
        list4 = dir(fullfile(trial_dir, '*_*_4_*_*.nii'));

        list_cell1 = struct2cell(list1);
        list_cell1 = list_cell1.';
        list_cell1 = list_cell1(:,[1]);
        for l = 1:length(list_cell1)
            tmp = strsplit(list_cell1{l,1}, '_');
            list_cell1{l,2} = str2double(tmp{1});
            %list_cell1{l,3} = tmp{3};
            if strcmp(tmp{2},strcat(hier{h}, '1'))
                list_cell1{l,4} = 1;
            else
                list_cell1{l,4} = 2;
            end
        end
        list_cell1 = horzcat(list_cell1(:,1), list_cell1(:,2), num2cell(ones([length(list_cell1), 1])), list_cell1(:,4));
        
        list_cell2 = struct2cell(list2);
        list_cell2 = list_cell2.';
        list_cell2 = list_cell2(:,[1]);
        for l = 1:length(list_cell2)
            tmp = strsplit(list_cell2{l,1}, '_');
            list_cell2{l,2} = str2double(tmp{1});
            %list_cell2{l,3} = tmp{3};
            if strcmp(tmp{2},strcat(hier{h}, '1'))
                list_cell2{l,4} = 1;
            else
                list_cell2{l,4} = 2;
            end
        end
        list_cell2 = horzcat(list_cell2(:,1), list_cell2(:,2), num2cell(2 * ones([length(list_cell2), 1])), list_cell2(:,4));
        
        list_cell3 = struct2cell(list3);
        list_cell3 = list_cell3.';
        list_cell3 = list_cell3(:,[1]);
        for l = 1:length(list_cell3)
            tmp = strsplit(list_cell3{l,1}, '_');
            list_cell3{l,2} = str2double(tmp{1});
            %list_cell3{l,3} = tmp{3};
            if strcmp(tmp{2},strcat(hier{h}, '1'))
                list_cell3{l,4} = 1;
            else
                list_cell3{l,4} = 2;
            end
        end
        list_cell3 = horzcat(list_cell3(:,1), list_cell3(:,2), num2cell(3 * ones([length(list_cell3), 1])), list_cell3(:,4));

        list_cell4 = struct2cell(list4);
        list_cell4 = list_cell4.';
        list_cell4 = list_cell4(:,[1]);
        for l = 1:length(list_cell4)
            tmp = strsplit(list_cell4{l,1}, '_');
            list_cell4{l,2} = str2double(tmp{1});
            %list_cell4{l,3} = tmp{3};
            if strcmp(tmp{2},strcat(hier{h}, '1'))
                list_cell4{l,4} = 1;
            else
                list_cell4{l,4} = 2;
            end
        end
        list_cell4 = horzcat(list_cell4(:,1), list_cell4(:,2), num2cell(4 * ones([length(list_cell4), 1])), list_cell4(:,4));
        
        alldesign = vertcat(list_cell1, list_cell2, list_cell3, list_cell4);
        save(sprintf('%s/train_test.mat', output_dir), 'alldesign')

        cfg.files.name = alldesign(:,1);
        cfg.files.chunk = cat(1, alldesign{:,2});
        cfg.files.label = cat(1, alldesign{:,3});
        cfg.files.xclass = cat(1, alldesign{:,4});

        cfg.scale.method = 'min0max1';
        cfg.scale.estimation = 'all';

        %cfg.scale.method = 'none'; % first disable scaling
        %cfg.scale.force_libsvm_no_scaling = 1; % then force that it stays disabled
        %cfg.scale.IKnowThatLibsvmCanBeSlowWithoutScaling = 1; % and acknowledge that we know that libsvm can be very slow without scaling

        % cfg.feature_selection.method = 'filter';
        % cfg.feature_selection.filter = 'F';
        % cfg.feature_selection.n_vox = 'automatic';

        %cfg.boot.n_boot = 10;
        %cfg.boot.balance_test = 1;
        cfg.design = make_design_xclass(cfg);
        cfg.design.unbalanced_data = 'ok';
        %cfg.basic_checks.DoubleFilenameEntriesOk = 1;

        cfg.results.output = {'accuracy','accuracy_minus_chance'}; % activate for alternative output
        cfg.results.overwrite = 1;
        results = decoding(cfg);
        
        matlabbatch = {};
   
        %%
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {strcat(subj_dir,'/anatomical/y_ana.nii')};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {strcat(output_dir,'/res_accuracy_minus_chance.nii')};
        %%
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                               78 76 85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

        spm_jobman('run',matlabbatch);
        

        matlabbatch = {};

        matlabbatch{1}.spm.util.imcalc.input = {
                                                strcat(output_dir,'/wres_accuracy_minus_chance.nii')
                                                MNI_brain_file
                                                };
        matlabbatch{1}.spm.util.imcalc.output = 'mwres_accuracy_minus_chance.nii';
        matlabbatch{1}.spm.util.imcalc.outdir = {output_dir};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*((i2)>0.1)';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;          
    
        spm_jobman('run',matlabbatch);

        matlabbatch = {};

        matlabbatch{1}.spm.spatial.smooth.data = {strcat(output_dir,'/wres_accuracy_minus_chance.nii')};
        matlabbatch{1}.spm.spatial.smooth.fwhm = FWHM;
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    
        spm_jobman('run',matlabbatch);
        
        matlabbatch = {};

        matlabbatch{1}.spm.util.imcalc.input = {
                                                strcat(output_dir,'/swres_accuracy_minus_chance.nii')
                                                MNI_brain_file
                                                };
        matlabbatch{1}.spm.util.imcalc.output = 'mswres_accuracy_minus_chance.nii';
        matlabbatch{1}.spm.util.imcalc.outdir = {output_dir};
        matlabbatch{1}.spm.util.imcalc.expression = 'i1.*((i2)>0.1)';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;          
    
        spm_jobman('run',matlabbatch);
    end
end