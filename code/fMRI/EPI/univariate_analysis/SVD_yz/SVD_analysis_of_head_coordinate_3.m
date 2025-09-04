clear all
clc

%% some warmup for the motivation :)

%use the rand function (uniform distribution), but it should work with
%other distributions as well
% t1 = rand(20,1);
% x  = 40 + t1*30;       % some right hemisphere positions
% 
% t1 = rand(20,1);
% x  = [x; -40 - t1*30];  % some left hemisphere positions
% 
% t1 = rand(40,1);
% y  = 20 + t1*40;       % rather frontal positions 
% 
% t1 = rand(40,1);
% z = 0.5*y+30 + t1*10;  % creating a tilted plane with noise
% 
% figure('Name','Simulated coordinates (clear left/right separation)');
% 
% plot3(x,y,z,' or');
% daspect([1 1 1]);
% xlabel('left  - right');
% ylabel('poste - anter');
% zlabel('infer - super');
% axis([-80 80 -20 80 0 100]);
% title('Simulated coordinates (clear left/right separation)');

base_dir = '/data/';    % Your project folder
project_dir = strcat(base_dir, 'Project1/MSHP/');
data_dir = strcat(project_dir, 'analysis/results_2nd/SVD_analysis_yz');

%% copying the simulated example into the matlab variable example

%example=zeros(size(x,1),6);
%example(:,1) = x;
%example(:,2) = y; 
%example(:,3) = z;

%% reading the file example.data (alternative to the previous point)
% 
importfilename_origin = strcat(data_dir, '/total_data_left_0_01_k0.txt');  % Hierarchy_contrast condition
%H_data = importdata(importfilename ,'');

fid = fopen(importfilename_origin, 'r');
H_data = textscan(fid, '%f %f %f %f %f sub%f', 'Delimiter', '\t'); % Adjust format and delimiter as needed
fclose(fid);

H_data = cell2mat(H_data);

for i = 1:length(H_data)
    k{i} = length(find(H_data(i,6) == H_data(:,6)));
end

index = [];

for j = 1:length(k)
    if k{j} == 3
        index(j) = 1;
    else
        index(j) = 0;
    end
end

% for j = 1:length(k)
%     index(j) = 1;     
% end

index_array = find(index==1);

for j = 1:length(index_array)
    choosen_H_data(j,:) = H_data(index_array(j),:);
end

final_H_data = sortrows(choosen_H_data, [4]);

% example = importdata(importfilename); % reading the coordinates
example = [final_H_data]; % reading the coordinates
n_example = size(final_H_data,1);
importfilename = 'Tcom_total.statdata';
% % just looking at the coordinates
%  example(:,1:3)
% 
 % selecting all points of the right hem and deleting them
%  rhidx=example(:,1)>0;
%  example(rhidx,:)=[];
%  example(:,1:3)
% 
% Mirror the points to the other hemisphere
%mexampl = example(:,1:3).*repmat([-1 1 1],size(example,1),1);
%example = [ example; [mexampl example(:,4:end)]];
%% First point of the SVD analysis (center the data)

% estimate the centerposition of all points and shiftthe origin
% to the centerposition
center_position = mean(example(:,2:3),1);
centered_pos    = example(:,2:3)-repmat(center_position,n_example,1);

%% SVD
% apply SVD
[U,S,V] = svd(centered_pos);
save('U.mat', 'U')
save('S.mat', 'S')
save('V.mat', 'V')
% delete the mirrored coordinates

%example(n_example+1:end,:)=[];

% U*S*V'+repmat(center_position,size(example,1),1)

singular_values = diag(S);             % 특이값 벡터
explained_variance = singular_values.^2;

total_variance = sum(explained_variance);
explained_variance_ratio = explained_variance / total_variance;

% take the first column for the statistics if the right hem is switched off
% take the second column if both hemisphere are analyzed
for_statistics = [centered_pos(1:n_example,:)*V example(:,4:end)];

[fipath,finame,~]=fileparts(importfilename_origin);
[fipath2,finame2,~]=fileparts(importfilename);
if isempty(fipath), fipath='.'; end;
fid=fopen([fipath '/' finame2 '.statdata'],'w');
fprintf(fid,'%2.3f %2.3f %8d %8d %8d\n',for_statistics');
fclose(fid);

%modification
set(gca,'FontName', 'Arial Narrow')
%set(gca,'defaultTextFontName', 'Helvetica')

%set(findobj(gcf, 'type', 'axes'), 'FontName', 'Helvetica', 'FontSize', 20);

subplot(1,1,1);  % right upper display
%NLidx=example(:,7)==0;  % non language data
%lhidx=example(NLidx,1)<0; % left hemisphereric data 
%lhexample=example(lhidx,:);
c1idx=example(:,4)==1;
c2idx=example(:,4)==2;
c3idx=example(:,4)==3;
%p.LineWidth = 3;
plot3(example(c1idx,1), example(c1idx,2),example(c1idx,3),'s', 'MarkerSize', 9, 'Color', '#676767', 'MarkerFaceColor','#2B6A6C');
hold on
%plot3(example(c1idx,1),example(c1idx,2),example(c1idx,3),' x', 'Color', '#8B488F');
%hold on
plot3(example(c2idx,1), example(c2idx,2),example(c2idx,3),'s', 'MarkerSize', 9, 'Color', '#676767', 'MarkerFaceColor','#F29724');
hold on
%plot3(example(c2idx,1),example(c2idx,2),example(c2idx,3),' x', 'Color', '#E6E007');
%hold on
plot3(example(c3idx,1), example(c3idx,2),example(c3idx,3),'s', 'MarkerSize', 9, 'Color', '#676767', 'MarkerFaceColor','#B80D48');
hold on
grid on;

cent = center_position;
cent_to_vector1 = cent + 20 * V(1:3);
cent_to_vector2 = cent + 20 * V(4:6);
cent_to_vector3 = cent + 20 * V(7:9);
% vectarrow(cent,cent_to_vector1)
% hold on
% vectarrow(cent,cent_to_vector2)
% hold on
% vectarrow(cent,cent_to_vector3)
% hold on

example = [ example; [mexampl example(:,4:end)]];

P=[mean(example(:,1)), mean(example(:,2)), mean(example(:,3))]; %mean
mX = [example(:,1)-P(1),example(:,1)-P(2),example(:,1)-P(3)]; %subtraction data, X - mean

t=[-100:100]; %t is range.
SS=V(:,2)*t;
%line equation
A=center_position(1)+SS(1,:);
B=center_position(2)+SS(2,:);
C=center_position(3)+SS(3,:);
%drawing
plot3(A,B,C, 'k-');

[n m] =size(example);
x0 = P(1);
y0 = P(2);
z0 = P(3);
a = V(1,2);
b = V(2,2);
c = V(3,2);

example(n_example+1:end,:)=[];

terror=0;
for i=1:37
    x = example(i,1);
    y = example(i,2);
    z = example(i,3);
    
    %get corss point on the line and that is of normal direction of 3d point
    t = -(a*x0 - a*x + b*y0 - b*y + c*z0 - c*z) / (a^2 + b^2 + c^2);    
    lx = x0 + t*a ;
    ly = y0 + t*b ;
    lz = z0 + t*c ;

    plot3(lx,ly,lz, '*', 'MarkerSize', 12, 'Color', '#2B6A6C'); %%point on the line - A
    %plot3(x,y,z,'k+'); %%point of the data - B
    
    plot3([x lx],[y ly],[z lz],':', 'Color', '#000000', 'linewidth', 0.3); %line A to B
    p.LineWidth = 0.5;
    d1 = ( (lx-x)^2+(ly-y)^2+(lz-z)^2 ); %uclidian distance between A and B
    terror= d1 + terror;
    
end

for i=38:74
    x = example(i,1);
    y = example(i,2);
    z = example(i,3);
    
    %get corss point on the line and that is of normal direction of 3d point
    t = -(a*x0 - a*x + b*y0 - b*y + c*z0 - c*z) / (a^2 + b^2 + c^2);    
    lx = x0 + t*a ;
    ly = y0 + t*b ;
    lz = z0 + t*c ;

    plot3(lx,ly,lz,'*', 'MarkerSize', 12, 'Color', '#F29724'); %%point on the line - A
    %plot3(x,y,z,'k+'); %%point of the data - B
    plot3([x lx],[y ly],[z lz],':', 'Color', '#000000', 'linewidth', 0.3); %line A to B
    
    d1 = ( (lx-x)^2+(ly-y)^2+(lz-z)^2 ); %uclidian distance between A and B
    terror= d1 + terror;
    
end

for i=75:111
    x = example(i,1);
    y = example(i,2);
    z = example(i,3);
    
    %get corss point on the line and that is of normal direction of 3d point
    t = -(a*x0 - a*x + b*y0 - b*y + c*z0 - c*z) / (a^2 + b^2 + c^2);    
    lx = x0 + t*a ;
    ly = y0 + t*b ;
    lz = z0 + t*c ;

    plot3(lx,ly,lz, '*', 'MarkerSize', 12, 'Color', '#B80D48'); %%point on the line - A
    %plot3(x,y,z,'k+'); %%point of the data - B
    plot3([x lx],[y ly],[z lz],':', 'Color', '#000000', 'linewidth', 0.3); %line A to B
    
    d1 = ( (lx-x)^2+(ly-y)^2+(lz-z)^2 ); %uclidian distance between A and B
    terror= d1 + terror;
    
end

axis([-60 20 -30 80 -30 80]);
view([-1 0 0]);

pbaspect([1 1 1])

axis equal

yticks(-30:55:80);
zticks(-30:55:80);

set(gca,'Color','#FFFFFF')
set(gca, 'linewidth', 1.1)
grid on