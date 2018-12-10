% Calculates the best fit to the original image from a timeseries

clearvars
% close all

% Read in files
inv_file = 'new_diag_block2';
inv_folder = 'paper/simulations l2 data';
username = 'bgsvisluke2';
tgv_lag_list = {'15','30'}; % TGV lamda_2 (*10)

n = 8; % index of last tgv itr
ref_flag = 0;
append_txt = '_l1_data_gc3';
lag = '15';

% % % Load in data
if ref_flag == 1; ref_text = '_ref'; else ref_text = ''; end;
l1_file = ['C:\Users\',username,'\Box Sync\Coupled TGV results\',inv_folder,'\',inv_file,'_l1_',lag,'lag',ref_text,append_txt];
l2_file = ['C:\Users\',username,'\Box Sync\Coupled TGV results\',inv_folder,'\',inv_file,'_l2_',lag,'lag',ref_text,append_txt];
background_file = [inv_file, ref_text];% 'vf_2g';%

load(l1_file,'final');final_l1 = final; clear final;load(l2_file,'final');final_l2 = final; clear final;
load(['C:\Users\',username,'\Box Sync\Data\New TGV test models\',background_file,'.mat']);
baseline = log10(ip4di_im);


% PRE CALCULATE DIMENSIONS
% calculate x and y dimensions of mesh (in pixels).
[xgrid, ygrid] = meshgrid(unique(final_l1.param_x), unique(final_l1.param_y));
xfind = find(final_l1.param_x == final_l1.param_x(1),2);
len_xdim = xfind(2) - xfind(1);
len_ydim = final_l1.num_param./len_xdim;
% Translate so that plot cells are center on parameter values, instead of
% holding those values on nodes
x = unique(final_l1.param_x); x = (x + [0; x(1:end-1)])/2;
y = unique(final_l1.param_y); y = (y + [0; y(1:end-1)])/2;

% Find sensitivity
if isfield(final_l1,'half_space_jac')
    cum_sens = mean(abs(final_l1.half_space_jac),1);
elseif isfield(final_l2,'half_space_jac')
    cum_sens = mean(abs(final_l2.half_space_jac),1);
else
    cum_sens = 0*final_l1.res_param1(:,1);
end
cum_sens_image = reshape(cum_sens,len_xdim,len_ydim)'./max(cum_sens);
sens_mask = cum_sens_image > 0.05;
sens_bg = nonzeros(baseline.*sens_mask);

% CALCULATE METRICS

% L1, ignores uniform initial model
for i = 2:length(final_l1.RMS)
    
    % form image
    l1_im = reshape(log10(final_l1.res_param1(:,i)),len_xdim,len_ydim)';
    sens_image = nonzeros(l1_im.*sens_mask); % for pearson
    
    % Calculate metrics
    [~,ssim_im] = ssim(l1_im,baseline,'exponents',[0,0,1]);
    l1_ssim(i) = sum(sum(ssim_im.*sens_mask))./sum(sens_mask(:));
    l1_r2(i) = sum((sens_image - mean(sens_image)).*(sens_bg - mean(sens_bg)))./sqrt(sum(sens_image.^2 - mean(sens_image).^2).*sum(sens_bg.^2 - mean(sens_bg).^2)); % pearson r2
    l1_r2b(i) = sum(sens_image.*sens_bg)./sqrt(sum(sens_image.^2).*sum(sens_bg.^2)); % r2 without mean centering
    rms_l1(i) = final_l1.RMS(i);
    
    %         % L1 figure
    %         figure(1)
    %         subplot(ceil((length(final_l1.RMS)-1)/2),2,i-1)
    %         surf(x,y,l1_im,'edgecolor','none')
    %         view([0,0,1])
    %         title(['L1, itr = ',num2str(i),', RMS = ',num2str(rms_l1(i)),', SSIM = ',num2str(l1_ssim(i)),', r^2 = ',num2str(l1_r2(i))],'interpreter','none')
    %         set(gca,'ydir','reverse')
    %         colorbar
    %         axis image
    %         colormap jet
    
end

% L2, ignores uniform initial model
for i = 2:length(final_l2.RMS)
    
    % form image
    l2_im = reshape(log10(final_l2.res_param1(:,i)),len_xdim,len_ydim)';
    sens_image = nonzeros(l2_im.*sens_mask); % for pearson
    
    % Calculate metrics
    [~,ssim_im] = ssim(l2_im,baseline,'exponents',[0,0,1]);
    l2_ssim(i) = sum(sum(ssim_im.*sens_mask))./sum(sens_mask(:));
    l2_r2(i) = sum((sens_image - mean(sens_image)).*(sens_bg - mean(sens_bg)))./sqrt(sum(sens_image.^2 - mean(sens_image).^2).*sum(sens_bg.^2 - mean(sens_bg).^2)); % pearson r2
    l2_r2b(i) = sum(sens_image.*sens_bg)./sqrt(sum(sens_image.^2).*sum(sens_bg.^2)); % r2 without mean centering
    rms_l2(i) = final_l2.RMS(i);
    
    %         % L2 figure
    %         figure(2)
    %         subplot(ceil((length(final_l2.RMS)-1)/2),2,i-1)
    %         surf(x,y,l2_im,'edgecolor','none')
    %         view([0,0,1])
    %         title(['L2, itr = ',num2str(i),', RMS = ',num2str(rms_l2(i)),', r^2 = ',num2str(l2_ssim(i)),', SSIM = ',num2str(l2_r2(i))],'interpreter','none')
    %         set(gca,'ydir','reverse')
    %         colorbar
    %         axis image
    %         colormap jet
    
end

% TGV
j_max = length(tgv_lag_list);

for j = 1:j_max
    
    tgv_lag = tgv_lag_list{j};
    tgv_path = ['C:\Users\',username,'\Box Sync\Coupled TGV results\',inv_folder,'\TGV ',inv_file,' lag',lag,'_tgv_lag',tgv_lag,ref_text,append_txt];

    for tgv_it = 1:n
        
        
        if exist(fullfile(tgv_path,['tgv_it_',num2str(tgv_it),'_sp2.mat']),'file') == 2
            load(fullfile(tgv_path,['tgv_it_',num2str(tgv_it),'_sp2.mat']))
        elseif exist(fullfile(tgv_path,['tgv_it_',num2str(tgv_it),'_sp1.mat']),'file') == 2
            load(fullfile(tgv_path,['tgv_it_',num2str(tgv_it),'_sp1.mat']))
        else
            disp(['Ending: filename ',num2str(j) ,' not found'])
            break
        end
        
        
        load(fullfile(tgv_path, ['\tgv_it_',num2str(tgv_it),'_sp1']))
        final_tgv = final; clear final;
        
        % form image
        tgv_im = reshape(log10(final_tgv.res_param1(:,end)),len_xdim,len_ydim)';
        sens_image = nonzeros(tgv_im.*sens_mask); % for pearson
        
        % Calculate metrics, save RMS
        [~,ssim_im] = ssim(tgv_im,baseline,'exponents',[0,0,1]);
        tgv_ssim(tgv_it) = sum(sum(ssim_im.*sens_mask))./sum(sens_mask(:));
        tgv_r2(tgv_it) = sum((sens_image - mean(sens_image)).*(sens_bg - mean(sens_bg)))./sqrt(sum(sens_image.^2 - mean(sens_image).^2).*sum(sens_bg.^2 - mean(sens_bg).^2)); % pearson r2
        tgv_r2b(tgv_it) = sum(sens_image.*sens_bg)./sqrt(sum(sens_image.^2).*sum(sens_bg.^2)); % r2 without mean centering
        rms_tgv(tgv_it) = final_tgv.RMS(end);
        
        %         % TGV figure
        %         figure(3)
        %         subplot(ceil(n/2),2,tgv_it)
        %         surf(x,y,tgv_im,'edgecolor','none')
        %         view([0,0,1])
        %         title(['TGV, itr = ',num2str(tgv_it),', RMS = ',num2str(rms_TGV(tgv_it)),', SSIM = ',num2str(tgv_ssim(tgv_it)),', r^2 = ',num2str(l1_r2(tgv_it))],'interpreter','none')
        %         set(gca,'ydir','reverse')
        %         colorbar
        %         axis image
        %         colormap jet
    end
    
        tgv.ssim{j} = tgv_ssim;
        tgv.r2{j} = tgv_r2;
        tgv.r2b{j} = tgv_r2b;
        tgv.rms{j} = rms_tgv;
        %     tgv.ssim_max(j,1:3) = [max(tgv_ssim), rms_tgv(find(tgv_ssim == max(tgv_ssim),1)), find(tgv_ssim == max(tgv_ssim),1)];
        %     tgv.r2_max(j,1:3) = [max(tgv_r2), rms_tgv(find(tgv_r2 == max(tgv_r2),1)), find(tgv_r2 == max(tgv_r2),1)];
        %     tgv.r2b_max(j,1:3) = [max(tgv_r2b), rms_tgv(find(tgv_r2b == max(tgv_r2b),1)), find(tgv_r2b == max(tgv_r2b),1)];
        tgv.ssim_max(j,1:10) = [str2double(tgv_lag) ,max(tgv_ssim), rms_tgv(find(tgv_ssim == max(tgv_ssim),1)), find(tgv_ssim == max(tgv_ssim),1),max(tgv_r2), rms_tgv(find(tgv_r2 == max(tgv_r2),1)), find(tgv_r2 == max(tgv_r2),1),max(tgv_r2b), rms_tgv(find(tgv_r2b == max(tgv_r2b),1)), find(tgv_r2b == max(tgv_r2b),1)];
        
        max_list(j,1:10) = [str2double(tgv_lag) ,max(tgv_ssim), rms_tgv(find(tgv_ssim == max(tgv_ssim),1)), find(tgv_ssim == max(tgv_ssim),1),max(tgv_r2), rms_tgv(find(tgv_r2 == max(tgv_r2),1)), find(tgv_r2 == max(tgv_r2),1),max(tgv_r2b), rms_tgv(find(tgv_r2b == max(tgv_r2b),1)), find(tgv_r2b == max(tgv_r2b),1)];
    clear tgv_ssim tgv_r2 tgv_r2b rms_tgv
end
    
    % Keep results for each iteration
    
    
    % max_list headings:
    % lambda2/l1/l2, SSIM, SSIM itr, SSIM RMS, r^2, R^2 itr, R^2 RMS, R^2, etc
    j = j_max;

    
    
    l1.ssim = l1_ssim;
    l1.r2 = l1_r2;
    l1.r2b = l1_r2b;
    l1.rms = rms_l1;
%     l1.ssim_max = [max(l1_ssim), rms_l1(find(l1_ssim == max(l1_ssim),1)), find(l1_ssim == max(l1_ssim),1)];
%     l1.r2_max = [max(l1_r2), rms_l1(find(l1_r2 == max(l1_r2),1)), find(l1_r2 == max(l1_r2),1)];
%     l1.r2b_max = [max(l1_r2b), rms_l1(find(l1_r2b == max(l1_r2b),1)), find(l1_r2b == max(l1_r2b),1)];
    max_list(j+1,1:10) = [771 ,max(l1_ssim), rms_l1(find(l1_ssim == max(l1_ssim),1)), find(l1_ssim == max(l1_ssim),1),max(l1_r2), rms_l1(find(l1_r2 == max(l1_r2),1)), find(l1_r2 == max(l1_r2),1),max(l1_r2b), rms_l1(find(l1_r2b == max(l1_r2b),1)), find(l1_r2b == max(l1_r2b),1)];
    
    l2.ssim = l2_ssim;
    l2.r2 = l2_r2;
    l2.r2b = l2_r2b;
    l2.rms = rms_l2;
%     l2.ssim_max = [max(l2_ssim), rms_l2(find(l2_ssim == max(l2_ssim),1)), find(l2_ssim == max(l2_ssim),1)];
%     l2.r2_max = [max(l2_r2), rms_l2(find(l2_r2 == max(l2_r2),1)), find(l2_r2 == max(l2_r2),1)];
%     l2.r2b_max = [max(l2_r2b), rms_l2(find(l2_r2b == max(l2_r2b),1)), find(l2_r2b == max(l2_r2b),1)];
    max_list(j+2,1:10) = [772, max(l2_ssim), rms_l2(find(l2_ssim == max(l2_ssim),1)), find(l2_ssim == max(l2_ssim),1), max(l2_r2), rms_l2(find(l2_r2 == max(l2_r2),1)), find(l2_r2 == max(l2_r2),1),max(l2_r2b), rms_l2(find(l2_r2b == max(l2_r2b),1)), find(l2_r2b == max(l2_r2b),1)];
    
    % SSIM best images
    [~, Smax_ind] = max(max_list(1:length(tgv_lag_list),2));
    S_max = max_list([end-1:end, Smax_ind],1:4);
    % Pearson best images
    [~, Rmax_ind] = max(max_list(1:length(tgv_lag_list),5));
    R_max = max_list([end-1:end, Rmax_ind], [1, 5:7]);
    
    disp('Best iteration: SSIM metric')
    disp('    lag2      SSIM     RMS error  Iterataion')
    disp(S_max)
    disp('Best iteration: Pearson r.^2 metric')
    disp('    lag2      r.^2      RMS error  Iterataion')
    disp(R_max)
    
    % Retrospectively apply RMS  threshold for convergence
    rms_thresh = 0.9;
    l1_itr = find(rms_l1(2:end) < rms_thresh,1); % -1 for last itr before change, +1 for 2:end
    l2_itr = find(rms_l2(2:end) < rms_thresh,1); % -1 for last itr before change, +1 for 2:end
    for i = 1:length(tgv_lag_list)
        tgv_itrs(i) = find(tgv.rms{i} < rms_thresh,1) - 1; % -1 for last itr before change
    end
%     disp('RMS threshholded itr:')
%     disp(['L1: ',num2str(l1_itr),'      L2: ',num2str(l2_itr),'      tgv: ',num2str(tgv_itrs)])
%     
%     
%     % threshold in absolute RMS difference
%     rms_diff_thresh = 0.1;
%     l1_itr = find(abs(diff(rms_l1)) < rms_diff_thresh,1); %+1 for diff, -1 for last itr before change
%     l2_itr = find(abs(diff(rms_l2)) < rms_diff_thresh,1);
%     for i = 1:length(tgv_lag_list)
%         tgv_itrs(i) = find(abs(diff(tgv.rms{i})) < rms_diff_thresh,1);
%     end
% %     disp('RMS diff threshholded itr:')
% %     disp(['L1: ',num2str(l1_itr),'      L2: ',num2str(l2_itr),'      tgv: ',num2str(tgv_itrs)])
%     
%         % threshold in pc RMS %difference
%     rms_pc_thresh = 0.2;
%     l1_pc_diff = abs(diff(rms_l1(3:end)) - diff(rms_l1(2:end-1)))./abs(diff(rms_l1(2:end-1)));
%     l2_pc_diff = abs(diff(rms_l2(3:end)) - diff(rms_l2(2:end-1)))./abs(diff(rms_l2(2:end-1)));
% 
%     l1_itr = find(l1_pc_diff < rms_pc_thresh,1); %+1 for diff, -1 for last itr before change
%     l2_itr = find(l2_pc_diff < rms_pc_thresh,1);
%     
%     for i = 1:length(tgv_lag_list)
%         rms_tgv = tgv.rms{i};
%         tgv_pc_diff = abs(diff(rms_tgv(2:end)) - diff(rms_tgv(1:end-1)))./abs(diff(rms_tgv(1:end-1)));
%         tgv_itrs(i) = min([find(tgv_pc_diff < rms_pc_thresh,1),length(tgv_pc_diff)]);
%     end
% %     disp('RMS diff threshholded itr:')
% %     disp(['L1: ',num2str(l1_itr),'      L2: ',num2str(l2_itr),'      tgv: ',num2str(tgv_itrs)])
%     
% 
% % % Plot metrics
% % % SSIM plot
% % figure(4)
% % plot(2:1+length(l1_ssim),l1_ssim)
% % hold on
% % plot(1:n,tgv_ssim)
% % plot(2:1+length(l2_ssim),l2_ssim)
% % hold off
% % title('SSIM')
% % legend('l1','l2','TGV')
% % xlabel('iteration')
% % set(gca,'ylim',[0.8,1.05])
% %
% % % r2 plot
% % figure(5)
% % plot(2:1+length(l1_r2),l1_r2)
% % hold on
% % plot(1:n,tgv_r2)
% % plot(2:1+length(l2_r2),l2_r2)
% % hold off
% % title('r.^2')
% % legend('l1','l2','TGV')
% % xlabel('iteration')
% % set(gca,'ylim',[0.8,1.05])
% %
% % % r2 plot
% % figure(6)
% % plot(2:1+length(l1_r2b),l1_r2b)
% % hold on
% % plot(1:n,tgv_r2b)
% % plot(2:1+length(l2_r2b),l2_r2b)
% % hold off
% % title('r.^2 no mean')
% % legend('l1','l2','TGV')
% % xlabel('iteration')
% % set(gca,'ylim',[0.8,1.05])
% 




% Select the final
