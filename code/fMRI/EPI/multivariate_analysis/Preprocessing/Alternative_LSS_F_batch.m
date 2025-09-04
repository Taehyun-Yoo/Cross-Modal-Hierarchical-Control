function Alternative_LSS_F_block_original_batch(s)

    base_dir = '/data/';    % Your project folder
    project_dir = strcat(base_dir, 'Project1/MSHP/');
    
    spm_dir = '/usr/local/MATLAB/R2023b/toolbox/spm12/';
    addpath(spm_dir)
        
    init = 1;
    total_sub_num = 53;    % Number of subjects
    total_sub_number = [init:1:init+total_sub_num-1];
    
    exc_sub_num = 11;
    exc_sub_number = [8, 11, 13, 20, 29, 30, 34, 38, 39, 47, 50];
    
    sub_num = total_sub_num - exc_sub_num;
    sub_number = setdiff(total_sub_number, exc_sub_number);
    
    run_num = 4;        % Number of runs
    
    scan = 221;
    scans = repelem(scan, run_num);
    
    units = 'secs';     % scans or secs
    TR = 2;             % Repetition time
    
    %slice_num = 42;     % Number of slices per volume
    highpassfilter = 128;       % high pass filter
    
    
    % In case of slice timing corrected data
    % you should modified microtime_resulution and microtime_onset.
    microtime_resolution = 16;  % number of slices
    microtime_onset = 8;       % half of number of slices
    
    condition_name = {'1', '2', '3', '4'};    % name of conditions
    condition_num = 4;
    
    subjects = {};
    
    for i = 1:length(sub_number)
        subjects{i} = sprintf('sub%02d', sub_number(i));
    end

    % subject specific directory
    subj_dir = fullfile(project_dir, 'analysis', sprintf('sub%02d', sub_number(s)));

    % onset
    for r = 1:run_num
        for c = 1:length(condition_name)
        % load txt file
            fid = fopen(strcat(project_dir, 'RAW/main/behavioral_data/Onset/', subjects{s}, '/', num2str(r+2), '_F4_', condition_name{c}, '.txt'), 'rt');
            onsets{r}{c} = cell2mat(textscan(fid, '%f', 'Delimiter', '\n') ) + scan * TR * (r-1);
            fclose(fid);
        end
    end

    for r = 1:run_num
        for c = 1:length(condition_name)
        % load txt file
            fid = fopen(strcat(project_dir, 'RAW/main/behavioral_data/Onset/', subjects{s}, '/', num2str(r+2), '_F4_', condition_name{c}, '_RT.txt'), 'rt');
            RTs{r}{c} = cell2mat(textscan(fid, '%f', 'Delimiter', '\n') );
            fclose(fid);
        end
    end
    
    for r = 1:run_num
        for c = 1:length(condition_name)
        % load txt file
            fid = fopen(strcat(project_dir, 'RAW/main/behavioral_data/Onset/', subjects{s}, '/', num2str(r+2), '_F4_', condition_name{c}, '_Block.txt'), 'rt');
            blocks{r}{c} = cell2mat(textscan(fid, '%f', 'Delimiter', '\n') );
            fclose(fid);
        end
    end

    for r = 1:run_num
        fid = fopen(strcat(project_dir, 'RAW/main/behavioral_data/Onset/', subjects{s}, '/', num2str(r+2), '_F4_', '1234.txt'), 'rt');
        onsets_all{r} = cell2mat(textscan(fid, '%f', 'Delimiter', '\n') ) + scan * TR * (r-1);
        fclose(fid);
    end
    onsets_all_hierarchy = vertcat(onsets_all{1}, onsets_all{2}, onsets_all{3}, onsets_all{4});

    for r = 1:run_num        
        for c = 1:condition_num
            correct_num_sub_run_con{r}{c} = length(onsets{r}{c});
            for k = 1:correct_num_sub_run_con{r}{c}
                if c == 1
                    trial_in{r}{c}{k} = [num2str(r) '_F4_1_' num2str(blocks{r}{c}(k)) '_' num2str(k)];
                elseif c == 2
                    trial_in{r}{c}{k} = [num2str(r) '_F4_2_' num2str(blocks{r}{c}(k)) '_' num2str(k)];
                elseif c == 3
                    trial_in{r}{c}{k} = [num2str(r) '_F4_3_' num2str(blocks{r}{c}(k)) '_' num2str(k)];
                else
                    trial_in{r}{c}{k} = [num2str(r) '_F4_4_' num2str(blocks{r}{c}(k)) '_' num2str(k)];
                end
                trial{r}{c} = trial_in{r}{c}.';
            end
        end
        correct_num_sub_run{r} = sum([correct_num_sub_run_con{r}{:}]);
        onsets_total{r} = cat(1, onsets{r}{:});
        RTs_total{r} = cat(1, RTs{r}{:});
        trial_total{r} = cat(1, trial{r}{:});
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% First level model specification

    % Modify the below lines
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    for r = 1:run_num
        for k = 1:correct_num_sub_run{r}
            mkdir(sprintf('%s/analysis/multivariate_analysis_2/F/trial_RT_block_original/%s',subj_dir, trial_total{r}{k}))
        end
    end

    %if isfile(strcat(base_dir, subjects{s}, '/analysis/SPM.mat'))
    %    delete(strcat(base_dir, subjects{s}, '/analysis/*'))
    %end
    % For R2017a and previous releases, use the "exist" function.
    for r = 1:run_num
        for k = 1:correct_num_sub_run{r}
            if exist(strcat(subj_dir, '/analysis/multivariate_analysis_2/F/trial_RT_block_original/', trial_total{r}{k}, '/SPM.mat'), 'file') == 2
                delete(strcat(subj_dir, '/analysis/multivariate_analysis_2/F/trial_RT_block_original/', trial_total{r}{k}, '/*'))
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    for r = 1:run_num
        for c = 1:correct_num_sub_run{r}
            matlabbatch = {};
        
            matlabbatch{1}.spm.stats.fmri_spec.dir = {strcat(subj_dir, '/analysis/multivariate_analysis_2/F/trial_RT_block_original/', trial_total{r}{c})};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = units;
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = microtime_resolution;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = microtime_onset;
        
            run_dir = fullfile(subj_dir, 'functional_concat', 'F', 'all');
        
            % input image
            subj_scans =  spm_select('FPList', run_dir, '^rf.*');
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(subj_scans);

            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = trial_total{r}{c};
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = onsets_total{r}(c);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = RTs_total{r}(c);
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;

            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'others';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = setdiff(onsets_all_hierarchy, onsets_total{r}(c));
            if r == 1
                tmp = RTs_total{1};
                tmp(c) = [];
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = vertcat(tmp, RTs_total{2}, RTs_total{3}, RTs_total{4});
            elseif r == 2
                tmp = RTs_total{2};
                tmp(c) = [];
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = vertcat(RTs_total{1}, tmp, RTs_total{3}, RTs_total{4});
            elseif r == 3
                tmp = RTs_total{3};
                tmp(c) = [];
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = vertcat(RTs_total{1}, RTs_total{2}, tmp, RTs_total{4});
            else
                tmp = RTs_total{4};
                tmp(c) = [];
                matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = vertcat(RTs_total{1}, RTs_total{2}, RTs_total{3}, tmp);
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
    
            % other settings for each session
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {spm_select('FPList', run_dir, '^rp.*')};
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = highpassfilter;
    
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
            
            spm_jobman('run',matlabbatch);

            spm_fmri_concatenate(strcat(subj_dir, '/analysis/multivariate_analysis_2/F/trial_RT_block_original/', trial_total{r}{c}, '/SPM.mat'), scans);

        %% First level model estimation
            matlabbatch = {};

            matlabbatch{1}.spm.stats.fmri_est.spmmat = {strcat(subj_dir, '/analysis/multivariate_analysis_2/F/trial_RT_block_original/', trial_total{r}{c}, '/SPM.mat')};
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    
            spm_jobman('run',matlabbatch);
        
            fprintf('1st-level model specification \n')
            fprintf('Subject: %d/%d \n', s, sub_num)
        end
    end
end