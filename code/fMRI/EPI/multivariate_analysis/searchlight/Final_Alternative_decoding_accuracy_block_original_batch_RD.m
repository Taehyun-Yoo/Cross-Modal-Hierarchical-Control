function Final_Alternative_decoding_accuracy_block_original_batch_RD(s)

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
        if exist(sprintf('%s/between_alternative_final',beta_dir), 'dir') ~= 7
            mkdir(sprintf('%s/between_alternative_final',beta_dir));
            %mkdir(sprintf('%s/between_alternative/trial_RT_block_original',beta_dir));
            
            list = dir(strcat(beta_dir, '/trial_RT_block_original/*_*_*_*_*/beta_0001.nii'));
            for l = 1:length(list)
                tmp = strsplit(list(l).folder, '/');
                filename = tmp{11};
                list(l).filename = filename;
                % if exist(sprintf('%s/between_alternative/trial_RT_block_original/%s.nii',beta_dir, list(l).filename), 'file') ~= 2
                %     copyfile(fullfile(list(l).folder, list(l).name), sprintf('%s/between_alternative/trial_RT_block_original/%s.nii',beta_dir, list(l).filename))
                % end
            end
        end

        if exist(sprintf('%s/between_alternative_mask/c1ana.nii',beta_dir), 'file') ~= 7
            copyfile(fullfile(subj_dir, 'anatomical/c1ana.nii'), sprintf('%s/between_alternative_final/c1ana.nii',beta_dir));
            
            matlabbatch = {};

            matlabbatch{1}.spm.util.imcalc.input = {
                                                    fullfile(list(1).folder, 'mask.nii')
                                                    strcat(beta_dir,'/between_alternative_final/c1ana.nii')
                                                    };
            matlabbatch{1}.spm.util.imcalc.output = 'gm_mask.nii';
            matlabbatch{1}.spm.util.imcalc.outdir = {strcat(beta_dir, '/between_alternative_final')};
            matlabbatch{1}.spm.util.imcalc.expression = 'i1.*((i2)>0.2)';
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;          
        
            spm_jobman('run',matlabbatch);
        end

        output_dir = fullfile(beta_dir, 'between_alternative_final');%% please change the directory if you want to put results in another folder
        trial_dir = fullfile(beta_dir, 'between_alternative_final', 'trial_RT_block_original');
        
        cd (fullfile(trial_dir))

        clear cfg
        
        cfg = [];
        
        cfg = decoding_defaults(cfg);

        cfg.analysis = decoding_type;

        cfg.results.dir = output_dir;

        cfg.decoding.method = 'classification_kernel'; % classification using the kernel speedup as standard
        cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; % linear classification
        %cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 100 -b 0 -q'; % linear classification
        %cfg.decoding.test.classification.model_parameters = '-q'; % linear classification
        %cfg.decoding.software = 'liblinear';
        %cfg.software = spm('ver');
 
        cfg.searchlight.radius = radius;
        cfg.searchlight.unit = 'voxels';
        cfg.searchlight.spherical = 1;

        %cfg.files.mask = fullfile(list(1).folder, 'mask.nii');
        cfg.files.mask = strcat(beta_dir, '/between_alternative_final/gm_mask.nii');
        %labelname1 = xclass_names{h}{1};
        %labelname2 = xclass_names{h}{2};
        %labelname3 = xclass_names{h}{3};
        %regressor_names = design_from_spm(beta_dir);
        
        lists = dir(fullfile(trial_dir, '*_*_*_*_*.nii'));
  
        list_cell = struct2cell(lists);
        list_cell = list_cell.';
        list_cell = list_cell(:,[1]);
        for l = 1:length(list_cell)
            tmp = strsplit(list_cell{l,1}, '_');
            list_cell{l,2} = str2double(tmp{1});
            list_cell{l,3} = str2double(tmp{3});
            if strcmp(tmp{2},strcat(hier{h}, '1'))
                list_cell{l,4} = 1;
            else
                list_cell{l,4} = 2;
            end
        end

        alldesign = {};
        alldesign = list_cell;
        save(sprintf('%s/train_test.mat', output_dir), 'alldesign')

        cfg.files.name = alldesign(:,1);
        cfg.files.chunk = cat(1, alldesign{:,2});
        cfg.files.label = cat(1, alldesign{:,3});
        cfg.files.xclass = cat(1, alldesign{:,4});
        
        cfg.scale.method = 'z';
        cfg.scale.estimation = 'across';
        cfg.scale.cutoff = [-3 3];

        %cfg.scale.method = 'min0max1';
        %cfg.scale.estimation = 'all';

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
        
        cfg.design.train = [];
        cfg.design.test = [];
        cv_num = 6;

        for cv = 1:cv_num
            if cv ~= 1
                cfg.design.label(:,cv) = cfg.design.label(:,1);
            end
            
            list_train = struct2cell(lists);
            list_train = list_train.';
            list_train = list_train(:,[1]);

            list_test = list_train;

            for l = 1:length(list_train)
                tmp = strsplit(list_train{l,1}, '_');
                % cv: 1 - train: 2, 4 in run 2 test: 1 in run 1
                if cv == 1
                    %train
                    if str2double(tmp{1}) == 2
                        if strcmp(tmp{2},strcat(hier{h}, '2')) || strcmp(tmp{2},strcat(hier{h}, '4'))
                            cfg.design.train(l,cv) = 1;
                        else
                            cfg.design.train(l,cv) = 0;
                        end
                    else
                        cfg.design.train(l,cv) = 0;
                    end
                    %test
                    if str2double(tmp{1}) == 1
                        if strcmp(tmp{2},strcat(hier{h}, '1'))
                            cfg.design.test(l,cv) = 1;
                        else
                            cfg.design.test(l,cv) = 0;
                        end
                    else
                        cfg.design.test(l,cv) = 0;
                    end
                % cv: 2 - train: 1, 4 in run 2 test: 2 in run 1
                elseif cv == 2
                    %train
                    if str2double(tmp{1}) == 2
                        if strcmp(tmp{2},strcat(hier{h}, '1')) || strcmp(tmp{2},strcat(hier{h}, '4'))
                            cfg.design.train(l,cv) = 1;
                        else
                            cfg.design.train(l,cv) = 0;
                        end
                    else
                        cfg.design.train(l,cv) = 0;
                    end
                    %test
                    if str2double(tmp{1}) == 1
                        if strcmp(tmp{2},strcat(hier{h}, '2'))
                            cfg.design.test(l,cv) = 1;
                        else
                            cfg.design.test(l,cv) = 0;
                        end
                    else
                        cfg.design.test(l,cv) = 0;
                    end
                % cv: 3 - train: 1, 2 in run 2 test: 4 in run 1
                elseif cv == 3
                    %train
                    if str2double(tmp{1}) == 2
                        if strcmp(tmp{2},strcat(hier{h}, '1')) || strcmp(tmp{2},strcat(hier{h}, '2'))
                            cfg.design.train(l,cv) = 1;
                        else
                            cfg.design.train(l,cv) = 0;
                        end
                    else
                        cfg.design.train(l,cv) = 0;
                    end
                    %test
                    if str2double(tmp{1}) == 1
                        if strcmp(tmp{2},strcat(hier{h}, '4'))
                            cfg.design.test(l,cv) = 1;
                        else
                            cfg.design.test(l,cv) = 0;
                        end
                    else
                        cfg.design.test(l,cv) = 0;
                    end
                % cv: 4 - train: 2, 4 in run 1 test: 1 in run 2
                elseif cv == 4
                    %train
                    if str2double(tmp{1}) == 1
                        if strcmp(tmp{2},strcat(hier{h}, '2')) || strcmp(tmp{2},strcat(hier{h}, '4'))
                            cfg.design.train(l,cv) = 1;
                        else
                            cfg.design.train(l,cv) = 0;
                        end
                    else
                        cfg.design.train(l,cv) = 0;
                    end
                    %test
                    if str2double(tmp{1}) == 2
                        if strcmp(tmp{2},strcat(hier{h}, '1'))
                            cfg.design.test(l,cv) = 1;
                        else
                            cfg.design.test(l,cv) = 0;
                        end
                    else
                        cfg.design.test(l,cv) = 0;
                    end
                % cv: 5 - train: 1, 4 in run 1 test: 2 in run 2
                elseif cv == 5
                    %train
                    if str2double(tmp{1}) == 1
                        if strcmp(tmp{2},strcat(hier{h}, '1')) || strcmp(tmp{2},strcat(hier{h}, '4'))
                            cfg.design.train(l,cv) = 1;
                        else
                            cfg.design.train(l,cv) = 0;
                        end
                    else
                        cfg.design.train(l,cv) = 0;
                    end
                    %test
                    if str2double(tmp{1}) == 2
                        if strcmp(tmp{2},strcat(hier{h}, '2'))
                            cfg.design.test(l,cv) = 1;
                        else
                            cfg.design.test(l,cv) = 0;
                        end
                    else
                        cfg.design.test(l,cv) = 0;
                    end
                % cv: 6 - train: 1, 2 in run 1 test: 4 in run 2 
                else
                    %train
                    if str2double(tmp{1}) == 1
                        if strcmp(tmp{2},strcat(hier{h}, '1')) || strcmp(tmp{2},strcat(hier{h}, '2'))
                            cfg.design.train(l,cv) = 1;
                        else
                            cfg.design.train(l,cv) = 0;
                        end
                    else
                        cfg.design.train(l,cv) = 0;
                    end
                    %test
                    if str2double(tmp{1}) == 2
                        if strcmp(tmp{2},strcat(hier{h}, '4'))
                            cfg.design.test(l,cv) = 1;
                        else
                            cfg.design.test(l,cv) = 0;
                        end
                    else
                        cfg.design.test(l,cv) = 0;
                    end
                end
            end
        end

        %cfg.design.label(:,2) = cfg.design.label(:,1);
        %cfg.design.train(:,2) = 1 - cfg.design.train(:,1);
        %cfg.design.test(:,2) = 1 - cfg.design.test(:,1);       
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