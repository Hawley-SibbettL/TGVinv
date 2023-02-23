clearvars
close all
file = 'gb2_long2_1pc';
inv_folder = '\thesis\periodic\golden_test\reg mesh\';
tgv_lagrn = '15';
itr = 9; % iteration of inverion to plot    
cum_sens_limit = 0.02;

lagrn = 'ls0';
baseline_folder =  'D:\TGV_revision\thesis\mesh_define\';
diff_weight_flag = 0; % 0 is even fd weighting

%load baseline and sensitivity model
baseline = load( fullfile(baseline_folder, file) );
climits = [1.8, 3.6];
load('gb2_long2_sens') % load cumulative sensitivity



pathname = ['D:\TGV_revision\',inv_folder,'\TGV ',file,' lag',lagrn,'_tgv_lag',tgv_lagrn];
curr_dir = pwd;
cd(pathname)
if exist(['tgv_it_',num2str(itr),'_sp2.mat'],'file') == 2
    load(['tgv_it_',num2str(itr),'_sp2.mat']);
    sp1_flag = 0;
%     if sp2_itr_flag == 1; sp2_itr(end+1) = final.sp2_itr(end); end
elseif exist(['tgv_it_',num2str(itr),'_sp1.mat'],'file') == 2
    load(['tgv_it_',num2str(itr),'_sp1.mat']);
    sp1_flag = 1;
else
    cd(curr_dir)
    error(['Ending: filename ',num2str(itr) ,' not found'])
end
cd(curr_dir);



x = unique(final.param_x);
y = unique(final.param_y);
[xgrid, ygrid] = meshgrid(x,y);
res_image = reshape(log10(final.res_param1(:,end)),length(x),length(y))';
cum_sens_im = ones(length(x),length(y));    % used for sensitivity modelling later
sum_factor = ones(size(final.param_x)); % can be used to isolate only part of the model


cum_sens_thresh = cum_sens_image > cum_sens_limit;
cum_sens_thresh_param = reshape( cum_sens_thresh', [], 1);
cum_sens_NaN = ones( size( cum_sens_image ) );
cum_sens_NaN( ~cum_sens_thresh ) = NaN;

if sp1_flag == 0
    % calculate objective terms
    [cx, cy] = c_calc_periodic( unique(final.param_x), unique(final.param_y), diff_weight_flag);
    l1_obj = sqrt( sum( cum_sens_thresh_param .* sum_factor .* ( cx * final.res_param1(:, end) - final.px ).^2 + ( cy * final.res_param1(:, end) - final.py ).^2 ) / sum( cum_sens_thresh_param ) );
    tgv_obj = sqrt( sum( cum_sens_thresh_param.*sum_factor .* ( cx' * final.px ).^2 + ( cy' * final.py ).^2 + 0.5 * ( cy' * final.px + cx' * final.py).^2 ) );
    
    px_im = reshape(final.px, length(x), [])';
    py_im = reshape(final.py, length(x), [])';
    z_cen = 100*ones(size(px_im));
    pz_im = zeros(size(px_im));
    
    figure
    surf(x,y,sqrt(px_im.^2 + py_im.^2),'edgecolor','none')
    view([0,0,1])
    title(['p, itr = ',num2str(i)],'interpreter','none')
    colorbar
    set(gca,'ydir','reverse')
    caxis([min(min(sqrt(px_im.^2 + py_im.^2))),max(max(sqrt(px_im.^2 + py_im.^2)))]);
    hold on
    quiver3(unique(final.param_x), unique(final.param_y), z_cen, px_im, py_im, pz_im,'color','m');
    axis image


else
    l1_obj = 0;
    tgv_obj = 0;
end
% calculate model misfit
rms_difference =  ( res_image - log10(baseline.ip4di_direct) ).^2;
rms_difference = sqrt( sum( sum( rms_difference( cum_sens_thresh ) ) ) / sum ( cum_sens_thresh_param ) );
misfit_text = num2str(rms_difference);



figure
surf(x,y, res_image .* cum_sens_NaN ,'edgecolor','none')
view([0,0,1])
title(['log(m), itr = ',num2str(itr),', RMS data misfit = ',num2str(final.RMS(end)),', reconstruction misfit = ', misfit_text, ', l1 term = ', num2str( l1_obj ), ', TGV term = ', num2str( tgv_obj ) ],'interpreter','none')
set(gca,'ydir','reverse','fontsize',16,'clim', climits)
colorbar
axis image
colormap parula