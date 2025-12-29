clear all
clc

base_dir = '/data/';    % Your project folder
project_dir = strcat(base_dir, 'Project1/MSHP/');

PFC_dir = strcat(project_dir, 'mask/');

spm_dir = '/usr/local/MATLAB/R2023b/toolbox/spm12/';
addpath(spm_dir)

addpath(strcat(project_dir, '/useful_scripts/etc'));

result_dir = strcat(project_dir, 'analysis/results_2nd/SVD_analysis_yz_2');

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

text_filenames = {};
count_file = 0;

for s = 1:length(sub_number)
    for h = 1:length(hier)       
        file = dir(strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis_2_-inf/', hier{h}, '_parametric/*_by_left_PFC_T_p0-05_unc_k0.txt'));
        if exist(strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis_2_-inf/', hier{h}, '_parametric/',file.name), 'file') == 2
            count_file = count_file + 1;
            text_filenames{1,count_file} = strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis_2_-inf/', hier{h}, '_parametric/',file.name);                
            text_filenames{2,count_file} = h;
            text_filenames{3,count_file} = subjects{s};
        end        
    end
end

opts = {};
T = {};
cluster_threshold = {};
x_coordinates = {};
y_coordinates = {};
z_coordinates = {};

data = {};

for t = 1:length(text_filenames)
    
    opts{t} = detectImportOptions(text_filenames{1,t});
    T{t} = readtable(text_filenames{1,t}, opts{t});
    n{t} = ~isnan(T{t}.cluster);
    
    T{t} = T{t}(n{t}(:,1),:);

    x_coordinates{t} = T{t}.Var12;
    y_coordinates{t} = T{t}.Var13;
    z_coordinates{t} = T{t}.Var14;
    
    not_left = find(x_coordinates{t} >= 0);
    not_left = flip(not_left);
    
    for nl = 1:length(not_left)
        T{t}(not_left(nl),:) = [];
    end

    largest_cluster_size = max(T{t}.cluster_2);
    largest_cluster_size_position = find(T{t}.cluster_2 == largest_cluster_size);
    
    final_x_coordinate = x_coordinates{t}(largest_cluster_size_position);
    final_y_coordinate = y_coordinates{t}(largest_cluster_size_position);
    final_z_coordinate = z_coordinates{t}(largest_cluster_size_position);

    lateralized = 1;
    
    data_one{t} = {final_x_coordinate, final_y_coordinate, final_z_coordinate, text_filenames{2,t}, lateralized, text_filenames{3,t}};
    data = vertcat(data, data_one{t});    
end

cell2csv([result_dir '/total_data_left_0_05_k0.txt'], data, '\t');