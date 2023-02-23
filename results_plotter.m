% Used to plot and compare results for different inversions
% function results_plotter

clearvars
close all


%% Compare same frame of two inversions with the same colormap (must have
% same dimensions as designed for comparing l2/l1 of same thing etc)

% DDa4n10 sc_575
inv_file = 'gb1';
append_txt = '_1pc';
background_file = [inv_file, '_refgrid'];
inv_folder1 = 'thesis\periodic\golden_test\66pc_cool'; % l1/TV 
inv_folder2 = 'thesis\periodic\golden_test\66pc_cool'; % l2
inv_folder3 = 'thesis\periodic\golden_test\66pc_cool'; % TGV location
username = getenv('username');
tgv_lag = '20'; % TGV lamda_2 (*10)
l1_it = 5; % [] for final itr, double otherwise
l2_it = 5; % [] for final itr, double otherwise
tgv_it = '8'; % string
ref_flag = 0;
baseline_flag = 1; % Baseline image exists
true_baseline_flag = 0; % use original fine resolution baseline (1), or interpolated version (0)
metric_flag = 1; % enables/disables metric calculations
crop_plot_flag = 0; % To cut off low sensitivity regions.
cum_sens_limit = 0.005;

l2_type = 'l2'; % l12 or l22 for l1/l2 data term
lag = 'ls0';
% % Plots column of resistivty for borehole comparison
% plot_ind = ;
% figure; plot(res_image1(:,plot_ind),y,'-x','markersize',15);hold on;plot(res_image2(:,plot_ind),y,'-x','markersize',15);plot(res_image3(:,plot_ind),y,'-x','markersize',15);hold off;set(gca,'ydir','reverse');legend('l1','l2','TGV','location','best');set(gca,'fontsize',16);
% hold on; plot([0.5,3],[4.5,4.5],'k--');hold off % mark bh layers


% if ref_flag == 1; ref_text = '_ref'; else ref_text = ''; end;

% FIND LAST TGV FILENAME
if isempty(tgv_it)
    for j = 15:-1:1
        if exist(['D:\TGV_revision\',inv_folder3,'\TGV ',inv_file,append_txt,' lag',lag,'_tgv_lag',tgv_lag,'\tgv_it_',num2str(j),'_sp1.mat'],'file') == 2
            tgv_it = num2str(j);
            break
        end
    end
end

% % % Load in data
file1 = ['D:\TGV_revision\',inv_folder1,'\',inv_file,append_txt,'_l1_',lag,'lag'];
file2 = ['D:\TGV_revision\',inv_folder2,'\',inv_file,append_txt,'_',l2_type,'_',lag,'lag'];
file3 = ['D:\TGV_revision\',inv_folder3,'\TGV ',inv_file,append_txt,' lag',lag,'_tgv_lag',tgv_lag,'\tgv_it_',tgv_it,'_sp1'];


% file3 = ['C:\Users\',username,'\OneDrive - The University of Nottingham\UoN Box Migration\Coupled TGV results\',inv_folder,'\TGV ',inv_file,' lag',lag,'_tgv_lag',tgv_lag,ref_text,append_txt,'\tgv_it_',tgv_it,'_sp1'];
% file2 = file1;


contour_flag = 0;
trace_flag = 0; % = 1 for trace dataset

% Paul tracer test models
% trace_data = importdata(fullfile(trace_path,'trench51_res_10_av.csv'));

%%
% load_background
% baseline = log10(VFHS_slant_gf);




if baseline_flag == 1
    try
        load(['D:\TGV_revision\thesis\mesh_define\',background_file,'.mat']);            
%           load(['C:\Users\',username,'\OneDrive - The University of Nottingham\UoN Box Migration\TGV_revision\Models\',background_file,'.mat']);            
    catch
        disp('**** Baseline not found ***')
    end
    baseline = log10(ip4di_direct);
%     baseline = baseline()
end

load(file1,'final');final1 = final; clear final;load(file2,'final');final2 = final; clear final;load(file3,'final');final3 = final; clear final;
% load(['C:\Users\',username,'\OneDrive - The University of Nottingham\UoN Box Migration\Coupled TGV results\resInvColormap'])


% Calculate x and y dimensions of mesh (in pixels).
[xgrid, ygrid] = meshgrid(unique(final1.param_x), unique(final1.param_y));
xfind = find(final1.param_x == final1.param_x(1),2);
len_xdim = xfind(2) - xfind(1);
len_ydim = final1.num_param./len_xdim;
% Translate so that plot cells are center on parameter values, instead of
% holding those values on nodes
x = unique(final1.param_x); x = (x + [0; x(1:end-1)])/2;
y = unique(final1.param_y); y = (y + [0; y(1:end-1)])/2;

if isempty(l1_it)
    l1_it = size(final1.res_param1,2);
end
if isempty(l2_it)
    l2_it = size(final2.res_param1,2);
end


res_image1 = reshape(log10(final1.res_param1(:,l1_it)),len_xdim,len_ydim)';
res_image2 = reshape(log10(final2.res_param1(:,l2_it)),len_xdim,len_ydim)';
res_image3 = reshape(log10(final3.res_param1(:,end)),len_xdim,len_ydim)';

% Calculate sensitivities
if isfield(final1,'half_space_jac') || isfield(final2,'half_space_jac')
    if isfield(final1,'half_space_jac')
        cum_sens = mean(abs(final1.half_space_jac),1);
    else
        cum_sens = mean(abs(final2.half_space_jac),1);
    end
    cum_sens_image = reshape(cum_sens,len_xdim,len_ydim)'./max(cum_sens);
    figure(20)
    surf(x,y,cum_sens_image.*(cum_sens_image>cum_sens_limit),'edgecolor','none')
    view([0,0,1])
    title(['Cumulative sensitivity > ',num2str(cum_sens_limit)],'interpreter','none')
    set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)])
    colorbar
    colormap parula
    axis image
else
    cum_sens_im = ones(len_xdim,len_ydim);
end


% visualise image/calculate difference images
if baseline_flag == 1
    ratio_image1 = (res_image1 - baseline);
    ratio_image2 = (res_image2 - baseline);
    ratio_image3 = (res_image3 - baseline);
    %     mean_ratio1 = mean(mean(abs(ratio_image1)));
    %     mean_ratio2 = mean(mean(abs(ratio_image2)));
    %     disp(['mean_ratio1 = ',num2str(mean_ratio1),' |  mean_ratio2 = ',num2str(mean_ratio2)])
end

if baseline_flag == 1
    climits = [min([min(min(res_image1)),min(min(res_image2)),min(min(res_image3)),min(min(baseline))]), max([max(max(res_image1)),max(max(res_image2)),max(max(res_image3)),max(max(baseline))])];
else
    climits = [min([min(min(res_image1)),min(min(res_image2)),min(min(res_image3))]), max([max(max(res_image1)),max(max(res_image2)),max(max(res_image3))])];
end

% For high sensitivty plot only
if crop_plot_flag == 1
    res_image1(cum_sens_image<cum_sens_limit) = NaN;
    res_image2(cum_sens_image<cum_sens_limit) = NaN;
    res_image3(cum_sens_image<cum_sens_limit) = NaN;
    baseline(cum_sens_image<cum_sens_limit) = NaN;
end
%% Plotting split comparisons

z = ones(size(res_image1));
% climits(1) = 1;
x_interp = linspace(0,x(end)+0.5*(x(end)-x(end-1)),33);
y_interp = linspace(0,y(end) - 0.5*abs(y(end)-y(end-1)),13);
[xg_interp, yg_interp] = meshgrid(x_interp,y_interp);

F1 = scatteredInterpolant(final1.param_x,final1.param_y,log10(final1.res_param1(:,l1_it)));
F2 = scatteredInterpolant(final1.param_x,final1.param_y,log10(final2.res_param1(:,l2_it)));
F3 = scatteredInterpolant(final1.param_x,final1.param_y,log10(final3.res_param1(:,end)));

interp_im1 = F1(xg_interp,yg_interp);
interp_im2 = F2(xg_interp,yg_interp);
interp_im3 = F3(xg_interp,yg_interp);

if baseline_flag == 1
    Fbase = scatteredInterpolant(final1.param_x,final1.param_y,reshape(baseline',1,[])');
    interp_base = Fbase(xg_interp,yg_interp);
end

n_contours = 18;

% Contourplot
if contour_flag == 1
    
    figure(9)
    subplot(2,2,1)
    contour(x_interp,y_interp,interp_im1,'fill','on','levelstep',(climits(2)-climits(1))/n_contours)
    set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)],...
        'clim',climits)
    title([file1,', RMS = ', num2str(final1.RMS(l1_it))],'interpreter','none')
    colorbar
    axis image
    
    subplot(2,2,2)
    contour(x_interp,y_interp,interp_im2,'fill','on','levelstep',(climits(2)-climits(1))/n_contours)
    set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)],'clim',climits)
    colormap parula
    title([file2,', RMS = ', num2str(final2.RMS(l2_it))],'interpreter','none')
    colorbar
    axis image
    
    subplot(2,2,3)
    contour(x_interp,y_interp,interp_im3,'fill','on','levelstep',(climits(2)-climits(1))/n_contours)
    set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)],'clim',climits)
    colormap parula
    title([file3,', RMS = ', num2str(final3.RMS(end))],'interpreter','none')
    colorbar
    axis image
    
    if baseline_flag == 1
        subplot(2,2,4)
        contour(x_interp,y_interp,interp_base,'fill','on','levelstep',(climits(2)-climits(1))/n_contours)
        set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)],'clim',climits)
        colormap parula
        title(['baseline'],'interpreter','none')
        colorbar
        axis image
    end
end


% Pixel Plot
figure(10)
h_pix(1) = subplot(4,1,3);
surf(x,y,res_image1,'edgecolor','none')
view([0,0,1])
title([file1,' RMS = ', num2str(final1.RMS(l1_it))],'interpreter','none')
set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)],'clim',climits)
colorbar
axis image

h_pix(2) = subplot(4,1,2);
surf(x,y,res_image2,'edgecolor','none')
view([0,0,1])
title([file2,' RMS = ',num2str(final2.RMS(l2_it))],'interpreter','none')
set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)],'clim',climits)
colorbar
axis image

h_pix(3) = subplot(4,1,4);
surf(x,y,res_image3,'edgecolor','none')
view([0,0,1])
title([file3,' RMS = ',num2str(final3.RMS(end))],'interpreter','none')
set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)],'clim',climits)
colorbar
axis image

colormap parula

if baseline_flag == 1
    % Plot base image
    h_pix(4) =subplot(4,1,1);
    if ~exist('thresh_im','var') || true_baseline_flag ~= 1
        surf(x,y,baseline,'edgecolor','none')
    else
        x2 = unique(X_edge);
        y2 = unique(Y_edge);
        x2 = [(x2(1:end-1) + x2(2:end))/2; x2(end) + (x2(end) - x2(end-1))/2];
        y2 = [(y2(1:end-1) + y2(2:end))/2; y2(end) + (y2(end) - y2(end-1))/2];
        surf(x2,y2,log10(thresh_im),'edgecolor','none')
    end
    view([0,0,1])
    title(['original image: '],'interpreter','none')
    colorbar
    axis image
    set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)],'clim',climits)
    
    % Ratioimage
    climits_ratio = [min(min(abs(ratio_image1(:))),min(abs(ratio_image2(:)))) , max(max(max(ratio_image1)),max(max(ratio_image2)))];
    % Calculate cumulative % error
    percent_error1 = sort(100 - 100*reshape(abs(ratio_image1-baseline)./(baseline),1,[]));
    percent_error2 = sort(100 - 100*reshape(abs(ratio_image2-baseline)./(baseline),1,[]));
    percent_error3 = sort(100 - 100*reshape(abs(ratio_image3-baseline)./(baseline),1,[]));
    % cum_err1 = 0*percent_error1;
    % cum_err2 = cum_err1;
    % for q = 1:max(length(percent_error1))
    %    cum_err1(q) = sum(percent_error1(1:q));
    %    cum_err2(q) = sum(percent_error2(1:q));
    % end
    
    % Ratio Pixel Plot
    figure(18)
    h_rat(1) = subplot(2,2,1);
    surf(x,y,abs(ratio_image1),'edgecolor','none')
    view([0,0,1])
    title([file1, ' RATIO'],'interpreter','none')
    set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)],'clim',climits_ratio)
    colorbar
    axis image
    colormap parula
    
    h_rat(2) = subplot(2,2,2);
    surf(x,y,abs(ratio_image2),'edgecolor','none')
    view([0,0,1])
    title([file2,' RATIO'],'interpreter','none')
    set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)],'clim',climits_ratio)
    colorbar
    axis image
    
    h_rat(3) = subplot(2,2,3);
    surf(x,y,abs(ratio_image3),'edgecolor','none')
    view([0,0,1])
    title([file3,' RATIO'],'interpreter','none')
    set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)],'clim',climits_ratio)
    colorbar
    axis image
    
    
    %     % ratio histograms
    %     figure(3)
    %     subplot(3,1,1);
    %     h1 = histogram(percent_error1,20);
    %     h1_ax = gca;
    %     title('Percentage difference in parameters (file1)','interpreter','none')
    %
    %     subplot(3,1,2)
    %     h2 = histogram(percent_error2,20);
    %     title('Percentage difference in parameters (file2)','interpreter','none')
    %     h2_ax = gca;
    %
    %     subplot(3,1,3)
    %     h3 = histogram(percent_error3,20);
    %     title('Percentage difference in parameters (file3)','interpreter','none')
    %     h3_ax = gca;
    %
    %     set(h1_ax,'ylim',[0, max([max(get(h3, 'values')),max(get(h2, 'values')),max(get(h1, 'values'))])*1.1],'xlim',[0,max([get(h3,'binlimits'),get(h1,'binlimits'), get(h2,'binlimits')])])
    %     set(h2_ax,'ylim',[0, max([max(get(h3, 'values')),max(get(h2, 'values')),max(get(h1, 'values'))])*1.1],'xlim',[0,max([get(h3,'binlimits'),get(h1,'binlimits'), get(h2,'binlimits')])])
    %     set(h3_ax,'ylim',[0, max([max(get(h3, 'values')),max(get(h2, 'values')),max(get(h1, 'values'))])*1.1],'xlim',[0,max([get(h3,'binlimits'),get(h1,'binlimits'), get(h2,'binlimits')])])
    %
    
    %     cum_err1 = zeros(1,get(h1,'numbins'));
    %     cum_err2 = zeros(1,get(h2,'numbins'));
    %     cum_err3 = zeros(1,get(h3,'numbins'));
    %
    %     bins1 = get(h1, 'values');
    %     bins2 = get(h2, 'values');
    %     bins3 = get(h3, 'values');
    %     binlim1 = get(h1,'binedges');
    %     binlim2 = get(h2,'binedges');
    %     binlim3 = get(h3,'binedges');
    %     for q = 1:length(cum_err1)
    %         cum_err1(q) = sum(bins1(1:q));
    %         cum_err2(q) = sum(bins2(1:q));
    %         cum_err3(q) = sum(bins2(1:q));
    %     end
    % %
    %     figure(4)
    %     plot(binlim1(2:end),cum_err1,'b-')
    %     hold on
    %     plot(binlim2(2:end),cum_err2,'r-')
    %     plot(binlim3(2:end),cum_err3,'g-')
    %     hold off
    %     xlabel('% error')
    %     ylabel('cumulative numer of pixels')
    %     legend('file1','file2','file3','location','best')
    %
    
end

figure(10)

%% ########## Reconstruction ######
disp([inv_file,'    TGV ',tgv_lag])
if baseline_flag == 1
    
    % anticline
    % solid_mask = (baseline == baseline(end));
    % trans_mask =  ~solid_mask;
    
    pc_image1 = ratio_image1./res_image1;
    pc_image2 = ratio_image2./res_image2;
    pc_image3 = ratio_image3./res_image3;
    
    % elipse
    %     solid_mask = (baseline == 1); % constant centre of ellipse
    %     trans_mask = (baseline ~=1) & imdilate(baseline ~= baseline(1),strel('disk',0)); % gradiated boundaries of ellipse
    %
    %     % Bedrock
    %     solid_mask = round(100*baseline) == round(100*(baseline(1)));
    %     dilate_rad = 1;
    %     trans_mask = (round(100*baseline)~=round(100*baseline(1))) & imdilate(round(100*baseline)~=round(100*baseline(end,end)),strel('disk',dilate_rad));
    %
    % vf_g
    %     solid_mask = (baseline == baseline(1));
    %     dilate_rad = 0;
    %     trans_mask = imdilate((baseline~=baseline(1)) & (baseline~=baseline(end,end)),strel('disk',dilate_rad));
    %
    % Grad
    %     solid_mask = imdilate((xgrid<12) & (baseline~=baseline(1)) & (cum_sens_image > 0.05),strel('line',3,0));
    % %     solid_mask = (xgrid<12) & (baseline~=baseline(1)) & (cum_sens_image > 0.05);
    %     trans_mask = (baseline~=baseline(end,end)) & (xgrid > 12);
    
    
    %     trans_rms_pc1 = 100*sqrt(sum(sum(trans_mask.*pc_image1.^2)))./sum(sum(trans_mask));
    %     trans_rms_pc2 = 100*sqrt(sum(sum(trans_mask.*pc_image2.^2)))./sum(sum(trans_mask));
    %     trans_rms_pc3 = 100*sqrt(sum(sum(trans_mask.*pc_image3.^2)))./sum(sum(trans_mask));
    %
    %     solid_rms_pc1 = 100*sqrt(sum(sum(solid_mask.*pc_image1.^2)))./sum(sum(solid_mask));
    %     solid_rms_pc2 = 100*sqrt(sum(sum(solid_mask.*pc_image2.^2)))./sum(sum(solid_mask));
    %     solid_rms_pc3 = 100*sqrt(sum(sum(solid_mask.*pc_image3.^2)))./sum(sum(solid_mask));
    %
    %     disp(['Solid mask RMS pcs: ',num2str(solid_rms_pc1),'% ',num2str(solid_rms_pc2),'% ',num2str(solid_rms_pc3),'% '])
    %     disp(['Trans mask RMS pcs: ',num2str(trans_rms_pc1),'% ',num2str(trans_rms_pc2),'% ',num2str(trans_rms_pc3),'% '])
    
    
    if exist('cum_sens_image','var') && metric_flag == 1
        sens_mask = cum_sens_image > cum_sens_limit;
        sens_image1 = nonzeros(res_image1.*sens_mask);
        sens_image2 = nonzeros(res_image2.*sens_mask);
        sens_image3 = nonzeros(res_image3.*sens_mask);
        sens_bg = nonzeros(baseline.*sens_mask);
        
        % RMS percentage change
        sens_rms_pc1 = 100*sqrt(sum(sum(sens_mask.*pc_image1.^2)))./sum(sum(sens_mask));
        sens_rms_pc2 = 100*sqrt(sum(sum(sens_mask.*pc_image2.^2)))./sum(sum(sens_mask));
        sens_rms_pc3 = 100*sqrt(sum(sum(sens_mask.*pc_image3.^2)))./sum(sum(sens_mask));
        
        % Pearson correlation coefficiant
        cor_p1 = sum((sens_image1 - mean(sens_image1)).*(sens_bg - mean(sens_bg)))./sqrt(sum(sens_image1.^2 - mean(sens_image1).^2).*sum(sens_bg.^2 - mean(sens_bg).^2));
        cor_p2 = sum((sens_image2 - mean(sens_image2)).*(sens_bg - mean(sens_bg)))./sqrt(sum(sens_image2.^2 - mean(sens_image2).^2).*sum(sens_bg.^2 - mean(sens_bg).^2));
        cor_p3 = sum((sens_image3 - mean(sens_image3)).*(sens_bg - mean(sens_bg)))./sqrt(sum(sens_image3.^2 - mean(sens_image3).^2).*sum(sens_bg.^2 - mean(sens_bg).^2));
        
        %     % Pearson image
        cor_p1_im = ((res_image1.*sens_mask - mean(sens_image1)).*(baseline.*sens_mask - mean(sens_bg)));%./sqrt(sum(sens_image1.^2 - mean(sens_image1).^2).*sum(sens_bg.^2 - mean(sens_bg).^2));
        cor_p2_im = ((res_image2.*sens_mask - mean(sens_image2)).*(baseline.*sens_mask - mean(sens_bg)));%./sqrt(sum(sens_image2.^2 - mean(sens_image2).^2).*sum(sens_bg.^2 - mean(sens_bg).^2));
        cor_p3_im = ((res_image3.*sens_mask - mean(sens_image3)).*(baseline.*sens_mask - mean(sens_bg)));%./sqrt(sum(sens_image3.^2 - mean(sens_image3).^2).*sum(sens_bg.^2 - mean(sens_bg).^2));
        
        % pearson no mean image
        %     cor_p1_im = ((res_image1.*sens_mask).*(baseline.*sens_mask))./sqrt(sum(sens_image1.^2).*sum(sens_bg.^2));
        %     cor_p2_im = ((res_image2.*sens_mask).*(baseline.*sens_mask))./sqrt(sum(sens_image2.^2).*sum(sens_bg.^2));
        %     cor_p3_im = ((res_image3.*sens_mask).*(baseline.*sens_mask))./sqrt(sum(sens_image3.^2).*sum(sens_bg.^2));
        
        % Pearson with mean removed
        cor_m1 = sum(sens_image1.*sens_bg)./sqrt(sum(sens_image1.^2).*sum(sens_bg.^2));
        cor_m2 = sum(sens_image2.*sens_bg)./sqrt(sum(sens_image2.^2).*sum(sens_bg.^2));
        cor_m3 = sum(sens_image3.*sens_bg)./sqrt(sum(sens_image3.^2).*sum(sens_bg.^2));
        
        
        % SSIM - calculated from image as a spatial map, then summed over
        % sensitivy region for equivalence (note that localisation does mean it
        % will include info from pixles outside region)
        [~,ssim_im1] = ssim(res_image1,baseline,'exponents',[1,1,1]);%[0,0,1]
        [~,ssim_im2] = ssim(res_image2,baseline,'exponents',[1,1,1]);
        [~,ssim_im3] = ssim(res_image3,baseline,'exponents',[1,1,1]);
        ssim1 = sum(sum(ssim_im1.*sens_mask))./sum(sens_mask(:));
        ssim2 = sum(sum(ssim_im2.*sens_mask))./sum(sens_mask(:));
        ssim3 = sum(sum(ssim_im3.*sens_mask))./sum(sens_mask(:));
        
        
        %     % Pearson comparison with SSIM
        %     cor_p1s = corr2(res_image1(1:10,:),baseline(1:10,:));
        %     cor_p2s = corr2(res_image2(1:10,:),baseline(1:10,:));
        %     cor_p3s = corr2(res_image3(1:10,:),baseline(1:10,:));
        
        
        disp(['Metrics for regions with cumulative sensitivty > ',num2str(cum_sens_limit)])
        disp(['RMS percentage change: ',num2str(sens_rms_pc1),'%   ',num2str(sens_rms_pc2),'%   ',num2str(sens_rms_pc3),'%   '])
        disp(['Pearson correlation coefficiant r^2: ',num2str(cor_p1),'   ',num2str(cor_p2),'   ',num2str(cor_p3),'   '])
        disp(['Pearson (without mean): ',num2str(cor_m1),'   ',num2str(cor_m2),'   ',num2str(cor_m3),'   '])
        %     disp('Metrics for depth<6m (includes low sensitivity edges)')
        %     disp(['Pearson correlation coefficiant: ',num2str(cor_p1s),'   ',num2str(cor_p2s),'   ',num2str(cor_p3s),'   '])
        disp(['SSIM: ',num2str(ssim1),'   ',num2str(ssim2),'   ',num2str(ssim3),'   '])
        
        sens_mask2 = sens_mask;
        sens_mask2(sens_mask2 == 0) = 10e10;
        ssim_climits = real([min(min([sens_mask2(:).*ssim_im1(:),sens_mask2(:).*ssim_im2(:),sens_mask2(:).*ssim_im3(:)])),max(max([sens_mask(:).*ssim_im1(:),sens_mask(:).*ssim_im2(:),sens_mask(:).*ssim_im3(:)]))]);
        cor_p_climits = [min(min([sens_mask2(:).*cor_p1_im(:),sens_mask2(:).*cor_p2_im(:),sens_mask2(:).*cor_p3_im(:)])),max(max([sens_mask(:).*cor_p1_im(:),sens_mask(:).*cor_p2_im(:),sens_mask(:).*cor_p3_im(:)]))];
        
        figure(55)
        colormap parula
        subplot(2,2,1);
        surf(x,y,abs(ssim_im1).*sens_mask,'edgecolor','none')
        view([0,0,1])
        title('SSIM images')
        set(gca,'ydir','reverse','clim',ssim_climits)
        colorbar
        axis image
        
        subplot(2,2,2);
        surf(x,y,abs(ssim_im2).*sens_mask,'edgecolor','none')
        view([0,0,1])
        set(gca,'ydir','reverse','clim',ssim_climits)
        colorbar
        axis image
        
        subplot(2,2,3);
        surf(x,y,abs(ssim_im3).*sens_mask,'edgecolor','none')
        view([0,0,1])
        set(gca,'ydir','reverse','clim',ssim_climits)
        colorbar
        axis image
        
        %     figure(6)
        %     colormap jet
        %     subplot(2,2,1);
        %     surf(x,y,cor_p1_im.*sens_mask,'edgecolor','none')
        %     view([0,0,1])
        %     title('Pearson Images')
        %     set(gca,'ydir','reverse','clim',cor_p_climits)
        %     colorbar
        %     axis image
        %
        %     subplot(2,2,2);
        %     surf(x,y,cor_p2_im.*sens_mask,'edgecolor','none')
        %     view([0,0,1])
        %     set(gca,'ydir','reverse','clim',cor_p_climits)
        %     colorbar
        %     axis image
        %
        %     subplot(2,2,3);
        %     surf(x,y,cor_p3_im.*sens_mask,'edgecolor','none')
        %     view([0,0,1])
        %     set(gca,'ydir','reverse','clim',cor_p_climits)
        %     colorbar
        %     axis image
        
    end
    
    %     figure(1)
    %     h_rat(4) = subplot(2,2,4);
    %     surf(x,y,trans_mask*1 + solid_mask*2,'edgecolor','none')
    %     view([0,0,1])
    %     title('masked regions: 1=trans, 2=solid','interpreter','none')
    %     set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)])
    %     colorbar
    %     axis image
    
end

%plot central section
% layer3 levels
% depth = [0, 0.8, 0.8, 1.9, 1.9, 5.3, 5.3, 8];
% rho = ([50, 50, 300, 300, 100, 100, 400, 400]);
% layerG levels - load from layersG

x_slice = 23;%23 for half in old format, 21 centre of gauss150

% block_edge2 cross section
yBE = [0, 0.9, 0.9, 1.2, 1.2, 3, 3, 3.3, 3.3, max(y)];
yBE2 = [0, 7.5, 7.5, 12, 12, 24, 24, 28.5, 28.5, max(y)];
resBE = [10, 10, 33, 33, 100, 100, 33, 33, 10. 10];

figure(8)
%     plot(baseline(:,x_slice),unique(final1.param_y),'-k')

if baseline_flag== 1 && exist('rho','var')
    [~, thresh_x_slice] = min(abs(x_centre(1:end-1) - final1.param_x(x_slice)));
    plot(log10(mean(thresh_im(1:end-1, thresh_x_slice-1:thresh_x_slice+1),2)),y_centre(1:end-1),'-k')
    elseif baseline_flag == 1
    plot(baseline(:,x_slice),unique(final1.param_y),'-k')
end
%y
%     plot(log10(resBE), yBE2,'k-')
%     plot(log10(rho(1:end-1)),depth(1:end-1),'-k.','markersize',8)
hold  on
plot(res_image1(:,x_slice),unique(final1.param_y),'--^')
plot(res_image2(:,x_slice),unique(final1.param_y),'--o')
plot(res_image3(:,x_slice),unique(final1.param_y),'--s')
set(gca,'ydir','reverse')
xlabel('log_{10}\bf{m}')
ylabel('z (m)')
title(['section at x = ',num2str(final1.param_x(x_slice))])
legend('baseline','l_1','l_2','TGV','location','best')
hold off

y_slice = 6;

figure(7)
%     plot(baseline(:,x_slice),unique(final1.param_y),'-k')
if baseline_flag== 1 && exist('rho','var')
    [~, thresh_y_slice] = min(abs(y_centre(1:end-1) - final1.param_y(1 + y_slice*(length(x)-1))));
    plot(x_centre(1:end-1),log10(mean(thresh_im(thresh_y_slice-1:thresh_y_slice+1, 1:end-1),1)),'-k')%y
elseif baseline_flag == 1
    plot(unique(final1.param_x), baseline(y_slice,:),'-k')
end
hold  on
plot(unique(final1.param_x),res_image1(y_slice,:),'--^')
plot(unique(final1.param_x),res_image2(y_slice,:),'--o')
plot(unique(final1.param_x),res_image3(y_slice,:),'--s')
ylabel('log_{10}\bf{m}')
xlabel('x (m)')
title(['section at y = ',num2str(final1.param_y(1 + y_slice*(length(x)-1)))])
legend('baseline','l_1','l_2','TGV','location','best')
hold off


% Add lines for cross sections
%  hold on; plot3([x(x_slice),x(x_slice)],[min(y),max(y)],[5,5],'k--')


%% big block model: histograms of resistivty within the block
% block_bin = 30;
% block_x = [1.5,2.3];
%
% figure
% subplot(4,1,2)
% histogram(reshape(res_image2(2:5,15:30),1,[]),block_bin);
% xlim(block_x)
% subplot(4,1,3)
% histogram(reshape(res_image1(2:5,15:30),1,[]),block_bin);
% xlim(block_x)
% subplot(4,1,4)
% histogram(reshape(res_image3(2:5,15:30),1,[]),block_bin);
% xlim(block_x)
%
%% Clifton layer plotting

%     plot(baseline(:,x_slice),unique(final1.param_y),'-k')
figure(5)
clifton_im_list = {res_image2,res_image1,res_image3};
clifton_title_list = {'l_2','l_1','TGV'};

for z = 1:3
    
        clifton_im = clifton_im_list{z};
    
        subplot(3,1,z)
        plot(unique(final1.param_x),clifton_im(2,:),'--^')
        hold on
        plot(unique(final1.param_x),clifton_im(3,:),'--d')
        plot(unique(final1.param_x),clifton_im(4,:),'--v')

        ylabel('log_{10}\bf{m}')
        xlabel('x (m)')
        title(clifton_title_list{z})
        hold off
        set(gca,'ylim',[1.5,3])
        legend(['z = ',num2str(final1.param_y(1 + 2*(length(x)-1)))],['z = ',num2str(final1.param_y(1 + 3*(length(x)-1)))],['z = ',num2str(final1.param_y(1 + 4*(length(x)-1)))],'location','best')
end

figure(10)