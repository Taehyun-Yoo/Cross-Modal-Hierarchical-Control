clear all
clc

base_dir = '/data/';    % Your project folder
project_dir = strcat(base_dir, 'Project1/MSHP/');

PFC_dir = strcat(project_dir, 'mask/');

spm_dir = '/usr/local/MATLAB/R2023b/toolbox/spm12/';
addpath(spm_dir)

addpath(strcat(project_dir, '/useful_scripts/etc'));

result_dir = strcat(project_dir, 'analysis/results_2nd/SVD_analysis');

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
        file = dir(strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis/', hier{h}, '_parametric/*_by_left_PFC_T_p0-01_unc_k0.txt'));
        if exist(strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis/', hier{h}, '_parametric/',file.name), 'file') == 2
            count_file = count_file + 1;
            text_filenames{1,count_file} = strcat(project_dir, 'analysis/', subjects{s}, '/analysis/parametric_analysis/', hier{h}, '_parametric/',file.name);                
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
    
    left = find(x_coordinates{t} < 0);
    
    if length(left) > 0
        final_x_coordinate = x_coordinates{t}(left(1));
        final_y_coordinate = y_coordinates{t}(left(1));
        final_z_coordinate = z_coordinates{t}(left(1));
    else
        final_x_coordinate = 'not left';
        final_y_coordinate = 'not left';
        final_z_coordinate = 'not left';
    end

    lateralized = 1;
    
    data_one{t} = {final_x_coordinate, final_y_coordinate, final_z_coordinate, text_filenames{2,t}, lateralized, text_filenames{3,t}};
    data = vertcat(data, data_one{t});    
end

index = [];

for j = 1:length(data)
    if data{j} ~= 'not left'
        index(j) = 1;
    else
        index(j) = 0;
    end
end

index_array = find(index==1);

for j = 1:length(index_array)
    choosen_data(j,:) = data(index_array(j),:);
end

cell2csv([result_dir '/total_data_left_0_01_k0.txt'], choosen_data, '\t');