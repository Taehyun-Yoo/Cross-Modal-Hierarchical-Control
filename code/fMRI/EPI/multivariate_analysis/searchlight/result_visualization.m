clear all
clc

base_dir = '/data/';    % Your project folder
project_dir = strcat(base_dir, 'Project1/MSHP/');

spm_dir = '/usr/local/MATLAB/R2023b/toolbox/spm12/';
addpath(spm_dir)
tool_dir = '/usr/local/MATLAB/R2023b/toolbox/CanlabCore-master/';
addpath(genpath(tool_dir))

overlay_img = strcat(project_dir, 'MNI152_T1_1mm_brain.nii');

hier = {'R', 'D'};
hier_num = 2;

save_img_dir = {};
result_dir = {};
result_img = {};
result_cl = {};

for h = 1:hier_num
    save_img_dir{h} = strcat(project_dir, 'analysis/mvpa_2_results_2nd/searchlight/', hier{h}, '/between_alternative_final');
    result_dir{h} = strcat(save_img_dir{h}, '/whole_brain');
    result_img{h} = sprintf('%s/average_accuracy.nii', result_dir{h});
    result_cl{h} = fmri_data(result_img{h}, 'noverbose');
end

cl_list = {result_cl{1}, result_cl{2}};

R_accuracy_range = [25 27.500996];
D_accuracy_range = [5.890348 27.973179];
range_accuracy_lists = {R_accuracy_range, D_accuracy_range};

save_all_img_dir = strcat(project_dir, 'analysis/mvpa_2_results_2nd/searchlight/final_images');

obj = fmridisplay;
obj = montage(obj, 'saggital', 'slice_range', [-24 -54], 'spacing', -6, 'overlay', overlay_img);

fig = gcf;
fig.Position = [100 100 1400 600]; % figure 크기 조정
set(fig, 'Color', 'w');  % 배경 흰색

cmin = 5.890348;
cmax = 11.838826;
N = 256; % colormap resolution

cb_bottom = 0.02; % y 위치 통일
cb_height = 0.13; % colorbar 높이
cb_spacing = 0.01; % colorbar 간격
cb_left_start = 0.6; % montage 오른쪽 시작 x좌표
cb_width_total = 0.3;

for h = 1:hier_num
    %1) blob 추가
    % obj = addblobs(obj, cl_list{h}, ...
    %                'maxcolor', color_lists{h}{1}, ...
    %                'mincolor', color_lists{h}{5}, ...
    %                'cmaprange', [cmin cmax]);
    obj = addblobs(obj, cl_list{h}, ...
                   'maxcolor', color_lists{h}{5}, ...
                   'mincolor', color_lists{h}{1}, ...
                   'cmaprange', range_FWE_lists{h});
    if h == 1
        tmp = 3;
    elseif h == 2
        tmp = 1;
    else
        tmp = 2;
    end
    % 2) 각 hierarchy colorbar width 계산
    cb_width = (cb_width_total - (hier_num-1)*cb_spacing)/hier_num;
    cb_left = cb_left_start + (tmp-1)*(cb_width + cb_spacing);
    
    % 3) colorbar axes 생성
    cb_ax = axes('Position', [cb_left cb_bottom cb_width cb_height]);
    cb_ax.Visible = 'off';

    % 4) hierarchy 전용 colormap 생성
    cmap = [linspace(color_lists{h}{1}(1), color_lists{h}{5}(1), N)', ...
            linspace(color_lists{h}{1}(2), color_lists{h}{5}(2), N)', ...
            linspace(color_lists{h}{1}(3), color_lists{h}{5}(3), N)'];
    colormap(cb_ax, cmap);

    % 5) colorbar 생성 (가로)
    c = colorbar(cb_ax, 'Orientation', 'vertical');
    c.Ticks = linspace(0,1,5); % normalized 0~1
    %c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), linspace(cmin, cmax,5),'UniformOutput',false);
    c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), linspace(range_FWE_lists{tmp}(1), range_FWE_lists{tmp}(2),5),'UniformOutput',false);
end

saveas(gcf, sprintf('%s/all_parametric_saggital_colorbar_range.png', save_all_img_dir));

%1-2)coronal
obj = fmridisplay;
obj = montage(obj, 'coronal', 'slice_range', [-3 37], 'spacing', 8, 'overlay', overlay_img);
for h = 1:hier_num
    %obj = addblobs(obj, cl_list{h}, 'maxcolor', color_lists{h}{1}, 'mincolor', color_lists{h}{5}, 'cmaprange', [5.890348 11.838826]);
    obj = addblobs(obj, cl_list{h}, 'maxcolor', color_lists{h}{5}, 'mincolor', color_lists{h}{1}, 'cmaprange', range_FWE_lists{h});
end
saveas(gcf, sprintf('%s/all_parametric_coronal_ver2_colorbar_range.png', save_all_img_dir));

obj = fmridisplay;
obj = montage(obj, 'coronal', 'slice_range', [-3 37], 'spacing', 8, 'overlay', overlay_img);

fig = gcf;
fig.Position = [100 100 1400 600]; % figure 크기 조정
set(fig, 'Color', 'w');  % 배경 흰색

cmin = 5.890348;
cmax = 11.838826;
N = 256; % colormap resolution

cb_bottom = 0.02; % y 위치 통일
cb_height = 0.13; % colorbar 높이
cb_spacing = 0.01; % colorbar 간격
cb_left_start = 0.6; % montage 오른쪽 시작 x좌표
cb_width_total = 0.3;

for h = 1:hier_num
    % 1) blob 추가
    % obj = addblobs(obj, cl_list{h}, ...
    %                'maxcolor', color_lists{h}{1}, ...
    %                'mincolor', color_lists{h}{5}, ...
    %                'cmaprange', [cmin cmax]);
    obj = addblobs(obj, cl_list{h}, ...
                   'maxcolor', color_lists{h}{5}, ...
                   'mincolor', color_lists{h}{1}, ...
                   'cmaprange', range_FWE_lists{h});
    if h == 1
        tmp = 3;
    elseif h == 2
        tmp = 1;
    else
        tmp = 2;
    end
    % 2) 각 hierarchy colorbar width 계산
    cb_width = (cb_width_total - (hier_num-1)*cb_spacing)/hier_num;
    cb_left = cb_left_start + (tmp-1)*(cb_width + cb_spacing);
    
    % 3) colorbar axes 생성
    cb_ax = axes('Position', [cb_left cb_bottom cb_width cb_height]);
    cb_ax.Visible = 'off';

    % 4) hierarchy 전용 colormap 생성
    cmap = [linspace(color_lists{h}{1}(1), color_lists{h}{5}(1), N)', ...
            linspace(color_lists{h}{1}(2), color_lists{h}{5}(2), N)', ...
            linspace(color_lists{h}{1}(3), color_lists{h}{5}(3), N)'];
    colormap(cb_ax, cmap);

    % 5) colorbar 생성 (가로)
    c = colorbar(cb_ax, 'Orientation', 'vertical');
    c.Ticks = linspace(0,1,5); % normalized 0~1
    %c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), linspace(cmin, cmax,5),'UniformOutput',false);
    c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), linspace(range_FWE_lists{tmp}(1), range_FWE_lists{tmp}(2),5),'UniformOutput',false);
end

saveas(gcf, sprintf('%s/all_parametric_coronal_colorbar_range.png', save_all_img_dir));

% %1-3)axial
% obj = fmridisplay;
% obj = montage(obj, 'axial', 'slice_range', [42 17], 'spacing', -5, 'overlay', overlay_img);
% for h = 1:hier_num
%     %obj = addblobs(obj, cl_list{h}, 'maxcolor', color_lists{h}{1}, 'mincolor', color_lists{h}{5}, 'cmaprange', [5.890348 11.838826]);
%     obj = addblobs(obj, cl_list{h}, 'maxcolor', color_lists{h}{5}, 'mincolor', color_lists{h}{1}, 'cmaprange', range_FWE_lists{h});
% end
% saveas(gcf, sprintf('%s/all_parametric_axial_ver2_colorbar_range.png', save_all_img_dir));
% 
% obj = fmridisplay;
% obj = montage(obj, 'axial', 'slice_range', [42 17], 'spacing', -5, 'overlay', overlay_img);
% 
% fig = gcf;
% fig.Position = [100 100 1400 600]; % figure 크기 조정
% set(fig, 'Color', 'w');  % 배경 흰색
% 
% cmin = 5.890348;
% cmax = 11.838826;
% N = 256; % colormap resolution
% 
% cb_bottom = 0.02; % y 위치 통일
% cb_height = 0.13; % colorbar 높이
% cb_spacing = 0.01; % colorbar 간격
% cb_left_start = 0.6; % montage 오른쪽 시작 x좌표
% cb_width_total = 0.3;
% 
% for h = 1:hier_num
%     % 1) blob 추가
%     % obj = addblobs(obj, cl_list{h}, ...
%     %                'maxcolor', color_lists{h}{1}, ...
%     %                'mincolor', color_lists{h}{5}, ...
%     %                'cmaprange', [cmin cmax]);
%     obj = addblobs(obj, cl_list{h}, ...
%                    'maxcolor', color_lists{h}{5}, ...
%                    'mincolor', color_lists{h}{1}, ...
%                    'cmaprange', range_FWE_lists{h});
%     if h == 1
%         tmp = 3;
%     elseif h == 2
%         tmp = 1;
%     else
%         tmp = 2;
%     end
%     % 2) 각 hierarchy colorbar width 계산
%     cb_width = (cb_width_total - (hier_num-1)*cb_spacing)/hier_num;
%     cb_left = cb_left_start + (tmp-1)*(cb_width + cb_spacing);
% 
%     % 3) colorbar axes 생성
%     cb_ax = axes('Position', [cb_left cb_bottom cb_width cb_height]);
%     cb_ax.Visible = 'off';
% 
%     % 4) hierarchy 전용 colormap 생성
%     cmap = [linspace(color_lists{h}{1}(1), color_lists{h}{5}(1), N)', ...
%             linspace(color_lists{h}{1}(2), color_lists{h}{5}(2), N)', ...
%             linspace(color_lists{h}{1}(3), color_lists{h}{5}(3), N)'];
%     colormap(cb_ax, cmap);
% 
%     % 5) colorbar 생성 (가로)
%     c = colorbar(cb_ax, 'Orientation', 'vertical');
%     c.Ticks = linspace(0,1,5); % normalized 0~1
%     %c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), linspace(cmin, cmax,5),'UniformOutput',false);
%     c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), linspace(range_FWE_lists{tmp}(1), range_FWE_lists{tmp}(2),5),'UniformOutput',false);
% end
% 
% saveas(gcf, sprintf('%s/all_parametric_axial_colorbar_range.png', save_all_img_dir));

% %1-4)3D with mask
% for hm = 1:hier_in_mask_num
%     obj = fmridisplay;
% 
%     surface_handles = surface(mask_template_cl, 'hcp inflated left', ...
%         'sourcespace', 'MNI152NLin2009cAsym', ...
%         'targetsurface', 'fsLR_32k', ...
%          'color', [1 1 1]);
% 
%      set(surface_handles, 'FaceAlpha', 1);
% 
%     surface_handles = surface(mask_cl{hm}, 'hcp inflated left', ...
%         'sourcespace', 'MNI152NLin2009cAsym', ...
%         'targetsurface', 'fsLR_32k', ...
%          'color', [1 1 1]);
% 
%      set(surface_handles, 'FaceAlpha', 1);
% 
%     surface_activation_handles = surface(hier_in_mask_cl{hm}, 'hcp inflated left', ...
%         'sourcespace', 'MNI152NLin2009cAsym', ...
%         'targetsurface', 'fsLR_32k', ...
%          'color', [1 0 0]);
% 
%      set(surface_activation_handles, 'FaceAlpha', 1);
% 
%     hold on;
% 
%     render_on_surface(hier_in_mask_cl{hm}, surface_handles, ...
%     'color', [1 0 0], ...         % 빨강 activation
%     'scaledtransparency', true);
% 
%     set(surface_handles, 'FaceAlpha', 0.5);
% 
%     set(gcf, 'Color', 'w');  % 배경 흰색
% end

% for hm = 1:hier_in_mask_num
%     obj = fmridisplay;
% 
%     surface_mask = surface(mask_cl{hm}, 'hcp inflated left', ...
%         'sourcespace', 'MNI152NLin2009cAsym', ...
%         'targetsurface', 'fsLR_32k', ...
%          'color', [1 1 1]);
% 
%     set(surface_mask, 'FaceAlpha', 0.5);
% 
%     hold on;
%     surface_activation = surface(hier_in_mask_cl{hm}, 'hcp inflated left', ...
%         'sourcespace', 'MNI152NLin2009cAsym', ...
%         'targetsurface', 'fsLR_32k', ...
%          'color', [1 0 0]);
% 
%     set(surface_activation, 'FaceAlpha', 1);
% 
%     set(gcf, 'Color', 'w');  % 배경 흰색
% end

%1-4)3D with mask
for h = 1:hier_num
    scn_export_papersetup(400);
    surface_handles = surface(result_cl{h}, 'hcp inflated left', ...
        'sourcespace', 'MNI152NLin2009cAsym', ...
        'targetsurface', 'fsLR_32k');
    set(gcf, 'Color', 'w');  % 배경 흰색
    set(surface_handles, 'FaceAlpha', 1);
    saveas(gcf, sprintf('%s/mask_%s.png', save_all_img_dir, mask{hm}));
    close all
end

for hm = 1:hier_in_mask_num   
    if ~isempty(hier_in_mask_cl{hm})
        if hier_in_mask{hm} == 'R'
            tmp = 2;
            range = 1;
        elseif hier_in_mask{hm} == 'F'
            tmp = 3;
            range = 2;
        else
            tmp = 1;
            range = 3;
        end
        
        obj = fmridisplay;
        obj = montage(obj, 'axial', 'slice_range', [56 6], 'spacing', -10, 'overlay', overlay_img);
        obj = addblobs(obj, hier_in_mask_cl{hm}, 'maxcolor', color_lists{tmp}{5}, 'mincolor', color_lists{tmp}{1}, 'cmaprange', range_unc_lists{range});
        %poscm = colormap_tor(color_lists{tmp}{5}, color_lists{tmp}{1}); % orange to yellow
        %scn_export_papersetup(400);
        % surface_activation_handles = surface(hier_in_mask_cl{hm}, 'hcp inflated left', ...
        %     'sourcespace', 'MNI152NLin2009cAsym', ...
        %     'targetsurface', 'fsLR_32k', ...
        %      'pos_colormap', poscm);
        scn_export_papersetup(400);
        set(gcf, 'Color', 'w');  % 배경 흰색
        %set(surface_activation_handles, 'FaceAlpha', 1);
        %clim([3.301273 11.838826]);
        %clim(range_unc_lists{range});
        saveas(gcf, sprintf('%s/%s_in_mask_%s_axial.png', save_all_img_dir, hier_in_mask{hm}, mask{hm}));
        close all
    end
end


surface(result_cl{1}, 'inflated left', ...
            'sourcespace', 'MNI152NLin2009cAsym', ...
            'targetsurface', 'fsLR_32k');


canlab_results_fmridisplay(result_img{2}, 'compact', 'outline', 'outlinecolor', [.5 0 .5])
canlab_results_fmridisplay(result_img{2}, 'montagetype', 'full')
canlab_results_fmridisplay(result_img{2}, 'regioncenters', 'outline', 'outlinecolor', [.5 0 .5])



figure;
surface_handles = [addbrain('hires left')];
set(surface_handles(1), 'FaceAlpha', .2);
ax = gca; % axes handle 저장
render_on_surface(result_cl{1}, surface_handles, 'EdgeColor', 'none');


lateralization = {'left', 'right'};
later_num = 2;

%2)3D
obj = fmridisplay;
for h = 1:hier_num
    for l = 1:later_num
        result_mask = fmri_data(result_img{h}, 'noverbose');
        %tmp = canlab_results_fmridisplay(result_mask);
        %cl = region(result_mask);  % region 객체로 변환
        %color_tmp = color_lists{h};
        %obj = addblobs(obj, obj.activation_maps{1}, 'splitcolor', color_tmp(1:4));
        %obj = addblobs(obj, cl, 'splitcolor', color_tmp(1:4));
        %obj.activation_maps{h} = canlab_results_fmridisplay(result_mask).activation_maps;
        %saveas(gcf, sprintf('%s/%s_parametric_PFC_coronal_slabs.png', save_img_dir, hier{h}));
        %close all
        if h == 1
            tmp = 2;
        elseif h == 2
            tmp = 3;
        else
            tmp = 1;
        end

        poscm = colormap_tor(color_lists{tmp}{5}, color_lists{tmp}{1}); % orange to yellow

        surface(result_mask, ['inflated ', lateralization{l}], 'pos_colormap', poscm, 'clim', range_FWE_lists{h});
        figs = findall(0, 'Type', 'figure');  % 열려 있는 figure 전부 찾기
        for f = 1:length(figs)
            ax = findall(figs(f), 'Type', 'axes');
            for a = 1:length(ax)
                view(ax(a), [0 90]);
            end
        end
        set(gcf, 'Color', 'w');  % 배경 흰색

        saveas(gcf, sprintf('%s/%s_inflated_%s.png', save_all_img_dir, hier{h}, lateralization{l}));
        close all
    end
end
%surface(result_mask, 'hcp inflated left');

for h = 1:hier_num
    save_img_dir = strcat(project_dir, 'analysis/results_2nd/parametric_analysis_2_-inf/', hier{h}, '_parametric_whole_brain');
    result_dir = strcat(save_img_dir, '/increasing_complexity');
    result_img = sprintf('%s/T_clusters_%s.nii', result_dir, hier{h});
    result_mask = fmri_data(result_img, 'noverbose');
    result_cl = mask2clusters(result_img);

    %scn_export_papersetup(400);
    figure;
    surface_handles = [addbrain('hires left')];
    set(surface_handles(1), 'FaceAlpha', .2);
    ax = gca; % axes handle 저장
    render_on_surface(result_mask, surface_handles, 'EdgeColor', 'none');
    axes(ax); % view 적용할 axes로 다시 이동
    view(0, 90);
    zoom(0.5);
    %saveas(gcf, sprintf('/home/taehyun/Downloads/BrainNet-Viewer-master/results/rsa_2nd_results/%s/%s/%s/tfce_all_rsa(subcortical).png', folder_name{f}, half_name{h}, class_name{c}));
    close all
end

for h = 1:hier_num
    save_img_dir = strcat(project_dir, 'analysis/results_2nd/parametric_analysis_2_-inf/', hier{h}, '_parametric_whole_brain');
    result_dir = strcat(save_img_dir, '/increasing_complexity');
    result_img = sprintf('%s/T_clusters_%s.nii', result_dir, hier{h});
    result_mask = fmri_data(result_img, 'noverbose');
    result_cl = mask2clusters(result_img);

    %scn_export_papersetup(400);
    surface_handles = surface(result_mask, 'hcp inflated left', 'sourcespace', 'MNI152NLin2009cAsym', 'targetsurface', 'fsLR_32k');
    %saveas(gcf, sprintf('/home/taehyun/Downloads/BrainNet-Viewer-master/results/rsa_2nd_results/%s/%s/%s/tfce_all_rsa(subcortical).png', folder_name{f}, half_name{h}, class_name{c}));
    close all   
end

for h = 1:hier_num
    save_img_dir = strcat(project_dir, 'analysis/results_2nd/parametric_analysis_2_-inf/', hier{h}, '_parametric_whole_brain');
    result_dir = strcat(save_img_dir, '/increasing_complexity');
    result_img = sprintf('%s/T_clusters_%s.nii', result_dir, hier{h});
    my_display_obj = canlab_results_fmridisplay([], 'sagittal', 'noverbose');
    my_display_obj = addblobs(my_display_obj, region(result_img));
end

create_figure('lateral surfaces');

surface_handles = [addbrain('foursurfaces')];
surface_handles = surface(result_mask, 'foursurfaces', 'noverbose');






obj = fmridisplay;
for h = 1:hier_num
    for l = 1:1
        result_mask = fmri_data(result_img{h}, 'noverbose');
        %result_mask.dat(result_mask.dat < 5.8903) = NaN;

        surface(result_mask, ['MNI152NLin6Asym midthickness ' lateralization{l}]);
        figs = findall(0, 'Type', 'figure');  % 열려 있는 figure 전부 찾기
        for f = 1:length(figs)
            ax = findall(figs(f), 'Type', 'axes');
            for a = 1:length(ax)
                patches = findall(ax(a), 'Type', 'patch');
                for p = 1:length(patches)
                    % FaceAlpha 조절 (투명도)
                    set(patches(p), 'FaceAlpha', 0.8);
                    %set(patches(p), 'FaceColor', 'flat');
                    % FaceColor 변경 가능
                    % set(patches(p), 'FaceColor', [1 0 0]);  % 빨간색 예시
                end
                view(ax(a), [-80 30]);
                camlight(ax(a), "left");  % 빛 추가
                lighting(ax(a), 'gouraud');    % 부드러운 shading
                material(ax(a), 'dull');       % 무광 표면
            end
        end

        set(gcf, 'Color', 'w');  % 배경 흰색

        %saveas(gcf, sprintf('%s/%s_rendering.png', save_all_img_dir, hier{h}));
        close all
    end
end