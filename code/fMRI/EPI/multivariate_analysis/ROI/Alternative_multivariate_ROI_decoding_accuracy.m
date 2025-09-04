clear all
clc

base_dir = '/data/';    % Your project folder
project_dir = strcat(base_dir, 'Project1/MSHP/');

spm_dir = '/usr/local/MATLAB/R2023b/toolbox/spm12/';
addpath(spm_dir)
tool_dir = '/usr/local/MATLAB/R2023b/toolbox/tdt_3.999H/';
addpath(genpath(tool_dir))
tool2_dir = strcat(project_dir, 'Code/EPI/multivariate_analysis/pilot/modified_toolbox');
addpath(genpath(tool2_dir))

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

%%%% mvpa type %%%%%%%%%%%%%%
decoding_type ='roi'; %%% determines decoding method ('searchlight','ROI', or 'wholebrain')

hier = {'R', 'D'};
hier_num = 2;

for s = 1:length(subjects) 
    subj_dir = fullfile(project_dir, 'analysis', sprintf('sub%02d', sub_number(s)));
    for h = 1: hier_num
        beta_dir = fullfile(subj_dir, 'analysis', 'multivariate_analysis_2', hier{h});
        
        if exist(sprintf('%s/between_alternative/mvpa_ROI',beta_dir), 'dir') ~= 7
            mkdir(sprintf('%s/between_alternative/mvpa_ROI',beta_dir));
        end

        roi_dir = strcat(beta_dir, '/mask/');
        roi_mask = dir(strcat(roi_dir, 'final_mvpa_*_ROI.nii'));
        %roi_mask(4:6,:) = [];

        trial_dir = fullfile(beta_dir, 'between_alternative', 'trial_RT_block_original');
        
        cd (fullfile(trial_dir))

        fi = fopen(fullfile(sprintf('%s/between_alternative/mvpa_ROI',beta_dir), sprintf('%s_decoding_accuracy_sub%02d.txt', hier{h}, sub_number(s))), 'w');
        for r = 1:length(roi_mask)
            if r == 1
                fprintf(fi, strcat('ROI accuracy','\n'));
            end
            tmp = strsplit(roi_mask(r).name, '_');
            tmp2 = strsplit(tmp{3}, '.');
            
            roi_name = strcat(tmp{2}, '_', tmp2{1});

            if exist(sprintf('%s/between_alternative/mvpa_ROI/%s',beta_dir, roi_name), 'dir') ~= 7
                mkdir(sprintf('%s/between_alternative/mvpa_ROI/%s',beta_dir, roi_name));
            end
        
            output_dir = fullfile(beta_dir, 'between_alternative', 'mvpa_ROI', roi_name);%% please change the directory if you want to put results in another folder

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
            
            %cfg.files.mask = fullfile(list(1).folder, 'mask.nii');
    
            cfg.files.mask = fullfile(roi_mask(r).folder, roi_mask(r).name);
            %labelname1 = xclass_names{h}{1};
            %labelname2 = xclass_names{h}{2};
            %labelname3 = xclass_names{h}{3};
            %regressor_names = design_from_spm(beta_dir);
            
            load(strcat(beta_dir, '/between_alternative/train_test.mat'));
    
            cfg.files.name = alldesign(:,1);
            cfg.files.chunk = cat(1, alldesign{:,2});
            cfg.files.label = cat(1, alldesign{:,3});
            cfg.files.xclass = cat(1, alldesign{:,4});
            
            cfg.scale.method = 'min0max1';
            cfg.scale.estimation = 'all';
    
            %cfg.scale.method = 'none'; % first disable scaling
            %cfg.scale.force_libsvm_no_scaling = 1; % then force that it stays disabled
            %cfg.scale.IKnowThatLibsvmCanBeSlowWithoutScaling = 1; % and acknowledge that we know that libsvm can be very slow without scaling
    
            %cfg.feature_selection.method = 'filter';
            %cfg.feature_selection.filter = 'F';
            %cfg.feature_selection.n_vox = 'automatic';
    
            %cfg.boot.n_boot = 10;
            %cfg.boot.balance_test = 1;
            cfg.design = make_design_xclass(cfg);
            cfg.design.unbalanced_data = 'ok';
            %cfg.basic_checks.DoubleFilenameEntriesOk = 1;
            
            cfg.results.output = {'accuracy','accuracy_minus_chance', 'predicted_labels', 'true_labels', 'decision_values'}; % activate for alternative output
            cfg.results.overwrite = 1;
            results = decoding(cfg);

            fprintf(fi, strcat([tmp{3} ' ' num2str(results.accuracy.output)], '\n'));
        end
        fclose(fi);
        close all
    end
end