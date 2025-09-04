clear all
clc

base_dir = '/data/';    % Your project folder
project_dir = strcat(base_dir, 'Project1/MSHP/');

PFC_dir = strcat(project_dir, 'mask/');

spm_dir = '/usr/local/MATLAB/R2023b/toolbox/spm12/';
addpath(spm_dir)

addpath(strcat(project_dir, '/useful_scripts/etc'));

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

hier = {'R', 'F', 'D'};

% first level statistics
multicom_cor_method = 'none';  % 'FWE' or 'none' or 'FDR'
thresh = 0.01;
extent = 0;

spm('defaults','fmri');
spm_jobman('initcfg');

for s = 1:length(sub_number)
    for h = 1:length(hier)
        pm_image_masked{s}{h} = strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis/', hier{h}, '_parametric/spmT_0001_masked.nii');   
        pm_image_binary{s}{h} = strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis/', hier{h}, '_parametric/spmT_0001_binary_mask.nii');   
        pm_image_cluster{s}{h} = strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis/', hier{h}, '_parametric/spmT_0001_global_cluster.nii');
    end
end

for s = 1:length(sub_number)
    for h = 1:length(hier)
        spm('defaults','FMRI');

        matlabbatch = {};
        matlabbatch{1}.spm.stats.results.spmmat = {strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis/', hier{h}, '_parametric/SPM.mat')};

        matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
        matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
        matlabbatch{1}.spm.stats.results.conspec.threshdesc = multicom_cor_method;
        matlabbatch{1}.spm.stats.results.conspec.thresh = thresh;
        matlabbatch{1}.spm.stats.results.conspec.extent = extent;
        matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec.mask.image.name = {strcat(PFC_dir, 'left_PFC.nii,1')};
        matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype = 0;
        matlabbatch{1}.spm.stats.results.units = 1;
        matlabbatch{1}.spm.stats.results.export{1}.ps = true;
        matlabbatch{1}.spm.stats.results.export{2}.pdf = true;
        matlabbatch{1}.spm.stats.results.export{3}.nary.basename = 'binary_mask';
        spm_jobman('run',matlabbatch);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        emptyflag = 0; % flag for noting if results are empty

        % generate filename
        filepath = xSPM.swd;
        filename = [xSPM.title '_' xSPM.STAT '_' xSPM.thresDesc '_k' num2str(xSPM.k)];

        % clean up filename
        filename = strrep(filename,' ','_');
        filename = strrep(filename,'0.','0-');
        removechars = {'.nii' '.img' '.' ',' '(' ')' '[' ']' '<' '>' '/' ':'};
        for i=1:length(removechars)
            filename = strrep(filename,removechars{i},'');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % generate SPM table
        TabDat = spm_list('List',xSPM,hReg);

        % get the table information and clean up
        d   = [TabDat.hdr;TabDat.dat];
        xyz = d(4:end,end);
        xyz2 = num2cell([xyz{:}]');

        % check whether there are clusters and if so, write out the results
        if isempty(xyz2)
            cell2csv([filepath '/' filename '_EMPTY.txt'],d,'\t');
            emptyflag = 1;

        else
            d(4:end,end:end+2) = xyz2;
            d(3,:)=[];

            % cell2csv from matlab file exchange
            cell2csv([filepath '/' filename '.txt'], d, '\t');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ~emptyflag % run the rest only for non-empty results

            % draw the overlay image
            image = slover([spm_dir 'canonical/single_subj_T1.nii']);
            image = add_spm(image);
            image.slices = [-48:4:84]; % these slices work for the above T1; otherwise adjust
            paint(image)

            % print image to file
            style = hgexport('factorystyle');
            style.Background = 'black';
            style.Format = 'png';
            hgexport(gcf, [filepath '/' filename '_slices.png'], style);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % save out the full image
            spm_write_filtered(xSPM.Z,xSPM.XYZ,SPM.xVol.DIM,SPM.xVol.M,'SPM-filtered',[filename]);
            
            % total_order = [];
            % size_table = size(TabDat.dat);
            % for k = 1:size_table(1,1)
            %     total_order = [total_order, TabDat.dat{k,5}];
            % end
            % [sortedvalue, sortorder] = sort(total_order);
            % global_cluster_position = find(sort(sortedvalue, 'descend') == total_order(1));
            % 
            % matlabbatch = {};
            % 
            % matlabbatch{1}.spm.util.imcalc.input = {
            %                                     pm_image_masked{s}{h}
            %                                     pm_image_binary{s}{h}
            %                                     };
            % matlabbatch{1}.spm.util.imcalc.output = pm_image_cluster{s}{h};
            % matlabbatch{1}.spm.util.imcalc.outdir = {strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis/', hier{h}, '_parametric')};
            % matlabbatch{1}.spm.util.imcalc.expression = sprintf('i1.*((i2)==%02d)', global_cluster_position);
            % matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            % matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            % matlabbatch{1}.spm.util.imcalc.options.mask = 1;
            % matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            % matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
            % 
            % spm_jobman('run',matlabbatch);
            
        end           
    end
end