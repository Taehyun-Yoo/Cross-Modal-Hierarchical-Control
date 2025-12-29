clear all
clc

spm('defaults','fmri');
spm_jobman('initcfg');


% Modify the below lines


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

subjects = {};

for i = 1:length(sub_number)
    subjects{i} = sprintf('sub%02d', sub_number(i));
end

% name of contrast

contrast_name = {'increasing complexity'};    

% second level contrast results folder
contrast_folder = {'increasing_complexity'}; 

mkdir(sprintf('%sanalysis/results_2nd/parametric_analysis_2_-inf/D_R_parametric_whole_brain',project_dir))
for c = 1:length(contrast_folder)
    mkdir(sprintf('%sanalysis/results_2nd/parametric_analysis_2_-inf/D_R_parametric_whole_brain/%s',project_dir,contrast_folder{c}))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %% Second level model specification
% for i = 1:length(contrast_name)
%     for s = 1:length(subjects)
%         second_contrasts{i}{s} = strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis/F_parametric/con_000', num2str(i), '.nii,1');
%     end
% end

% Estimation directory
contrast_dir = {};
Est_file = {};

for i=1:length(contrast_folder)
    contrast_dir{i} = strcat(project_dir, 'analysis/results_2nd/parametric_analysis_2_-inf/D_R_parametric_whole_brain/', contrast_folder{i});
    Est_file{i} = strcat(project_dir, 'analysis/results_2nd/parametric_analysis_2_-inf/D_R_parametric_whole_brain/', contrast_folder{i},'/SPM.mat');
end

for i=1:length(contrast_folder)
    matlabbatch = {};
    
    matlabbatch{1}.spm.stats.factorial_design.dir = {contrast_dir{i}};
    for s = 1:length(subjects)
        % F 대비 경로
        D_path = strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis_2_-inf/D_parametric/con_000', num2str(i), '.nii');
        % R 대비 경로
        R_path = strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis_2_-inf/R_parametric/con_000', num2str(i), '.nii');
    
        % 존재 확인 (없으면 명확히 에러)
        % if ~exist(F_path,'file')
        %     error('Missing F contrast for %s: %s', subjects{s}, F_path);
        % end
        % if ~exist(R_path,'file')
        %     error('Missing R contrast for %s: %s', subjects{s}, R_path);
        % end
    
        matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(s).scans = { ...
            [D_path, ',1'] ; ...
            [R_path, ',1'] };
    end
    matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca    = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova   = 0;
    matlabbatch{1}.spm.stats.factorial_design.cov             = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im      = 1;     % implicit mask = on
    matlabbatch{1}.spm.stats.factorial_design.masking.em      = {''};  % explicit mask (원하면 넣으세요)
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit  = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    
    
    %% Second level model estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {Est_file{i}};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    matlabbatch{2}.spm.stats.fmri_est.delete = 1;
    
    
    
    fprintf('2nd-level model specification \n')
    fprintf('Contrast: %d/%d \n', i, length(contrast_name))
    
    
    %% Second level contrast manager
    matlabbatch{3}.spm.stats.con.spmmat = {Est_file{i}};
    
    
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name    = 'D > R';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    % (옵션) R > F
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name    = 'R > D';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    % (옵션) 각 조건의 >0
    % matlabbatch{3}.spm.stats.con.consess{3}.tcon.name    = 'D > 0';
    % matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [1 0];
    % matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    % 
    % matlabbatch{3}.spm.stats.con.consess{4}.tcon.name    = 'F > 0';
    % matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 1];
    % matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.delete = 1;
    
    
    spm_jobman('run',matlabbatch);
    
    fprintf('2nd-level Contrast Manager \n')
    fprintf('Contrast: %d/%d \n', c, length(contrast_folder))

end
