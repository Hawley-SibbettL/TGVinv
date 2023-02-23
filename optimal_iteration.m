% Calculates the best fit to the original image from a timeseries

clearvars
% close all

% Read in files
inv_file = 'gb1';
background_file = [inv_file, '_refgrid'];
append_txt = '_1pc';
obj_flag = 0;

lag = 'ls0';
n = 7; % index of last tgv itr

inv_folder = 'thesis\periodic\golden_test';
tgv_lag_list = {'5','7.5','10','15','20'}; % TGV lamda_2 (*10)


% % % Load in data
l1_file = ['D:\TGV_revision\',inv_folder,'\',inv_file,append_txt,'_l1_',lag,'lag'];
l2_file = ['D:\TGV_revision\',inv_folder,'\',inv_file,append_txt,'_l2_',lag,'lag'];

load(l1_file,'final');final_l1 = final; clear final;load(l2_file,'final');final_l2 = final; clear final;
load(['D:\TGV_revision\thesis\mesh_define\',background_file,'.mat']);
baseline = log10(ip4di_direct);


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
sens_mask = cum_sens_image > 0.005;
sens_bg = nonzeros(baseline.*sens_mask);




% % Calculate the objective function terms
% % TV term: |grad(m)|
% grad_sumSq = (cx*log10(mesh.res_param1)).^2 + (cy*log10(mesh.res_param1)).^2;
% dx_im = reshape(cx*log10(mesh.res_param1), length(tmp_x), [])';
% dy_im = reshape(cy*log10(mesh.res_param1), length(tmp_x), [])';
% L2_reg = sqrt(sum(grad_sumSq(:)));
%
%
% TGV_t1_im = (sqrt((cx*log10(mesh.res_param1) - mesh.px).^2 + (cy*log10(mesh.res_param1) - mesh.py).^2));
% TGV_t2_im = (sqrt( (cx'*mesh.px).^2 + (cy'*mesh.py).^2 +  0.5*(cy'*mesh.py + cx'*mesh.px).^2));
% TGV_t1(i) = sum(TGV_t1_im(:));
% TGV_t2(i) = sum(TGV_t2_im(:));


[cx, cy] = c_calc_periodic(unique(final_l1.param_x), unique(final_l1.param_y), 1);




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
    l1_rmsError(i) = sqrt( mean( (l1_im(:) - baseline(:)).^2 ) );
    
    
    % calculate TV gradient
    grad_sumSq = (cx*log10(final_l1.res_param1(:,end))).^2 + (cy*log10(final_l1.res_param1(:,end))).^2;
    TV_reg = sum(sqrt(grad_sumSq(:)));
    
    
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
    [~,ssim_im] = ssim(l2_im,baseline,'exponents',[1,1,1]);
    l2_ssim(i) = sum(sum(ssim_im.*sens_mask))./sum(sens_mask(:));
    l2_r2(i) = sum((sens_image - mean(sens_image)).*(sens_bg - mean(sens_bg)))./sqrt(sum(sens_image.^2 - mean(sens_image).^2).*sum(sens_bg.^2 - mean(sens_bg).^2)); % pearson r2
    l2_r2b(i) = sum(sens_image.*sens_bg)./sqrt(sum(sens_image.^2).*sum(sens_bg.^2)); % r2 without mean centering
    rms_l2(i) = final_l2.RMS(i);
    l2_rmsError(i) = sqrt( mean( (l2_im(:) - baseline(:)).^2 ) );

    
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
    tgv_path = ['D:\TGV_revision\',inv_folder,'\TGV ',inv_file,append_txt,' lag',lag,'_tgv_lag',tgv_lag];
    
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
        [~,ssim_im] = ssim(tgv_im,baseline,'exponents',[1,1,1]);
        tgv_ssim(tgv_it) = sum(sum(ssim_im.*sens_mask))./sum(sens_mask(:));
        tgv_r2(tgv_it) = sum((sens_image - mean(sens_image)).*(sens_bg - mean(sens_bg)))./sqrt(sum(sens_image.^2 - mean(sens_image).^2).*sum(sens_bg.^2 - mean(sens_bg).^2)); % pearson r2
        tgv_r2b(tgv_it) = sum(sens_image.*sens_bg)./sqrt(sum(sens_image.^2).*sum(sens_bg.^2)); % r2 without mean centering
        rms_tgv(tgv_it) = final_tgv.RMS(end);
        tgv_rmsError(tgv_it) = sqrt( mean( (tgv_im(:) - baseline(:) ).^2 ) );
        
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
    tgv.rmsError{j} = tgv_rmsError;
    %     tgv.ssim_max(j,1:3) = [max(tgv_ssim), rms_tgv(find(tgv_ssim == max(tgv_ssim),1)), find(tgv_ssim == max(tgv_ssim),1)];
    %     tgv.r2_max(j,1:3) = [max(tgv_r2), rms_tgv(find(tgv_r2 == max(tgv_r2),1)), find(tgv_r2 == max(tgv_r2),1)];
    %     tgv.r2b_max(j,1:3) = [max(tgv_r2b), rms_tgv(find(tgv_r2b == max(tgv_r2b),1)), find(tgv_r2b == max(tgv_r2b),1)];
    tgv.ssim_max(j,1:10) = [str2double(tgv_lag) ,max(tgv_ssim), rms_tgv(find(tgv_ssim == max(tgv_ssim),1)), find(tgv_ssim == max(tgv_ssim),1),max(tgv_r2), rms_tgv(find(tgv_r2 == max(tgv_r2),1)), find(tgv_r2 == max(tgv_r2),1),max(tgv_r2b), rms_tgv(find(tgv_r2b == max(tgv_r2b),1)), find(tgv_r2b == max(tgv_r2b),1)];
    
    max_list(j,1:10) = [str2double(tgv_lag) ,max(tgv_ssim), rms_tgv(find(tgv_ssim == max(tgv_ssim),1)), find(tgv_ssim == max(tgv_ssim),1),max(tgv_r2), rms_tgv(find(tgv_r2 == max(tgv_r2),1)), find(tgv_r2 == max(tgv_r2),1),max(tgv_r2b), rms_tgv(find(tgv_r2b == max(tgv_r2b),1)), find(tgv_r2b == max(tgv_r2b),1)];
    
    
    
    % ------------ Calculate p -------------------------------------------
    if obj_flag == 1
        % asumes not using old neumann_bc ghost points
        num_param = final_tgv.num_param;
        domain_ind = 1:num_param;
        param_x = final_tgv.param_x;
        param_y = final_tgv.param_y;
        res_param = final_tgv.res_param1(:,end);
        tgv_lagrn(j) = str2num(tgv_lag)/10;
        
        grad_mx = cx*log10(res_param);
        grad_my = cy*log10(res_param);
        sp2_itr = 100;
        p_conv = 1e-5;
        
        
        %initialise p - includes bc ghost point terms
        px = zeros(num_param,1);
        py = px;
        
        
        % Used in main calculation
        gamma_p = 1e-12;               % cutoff for |x| -> 0
        gamma_c = 1e-15;
        % Used in performance metrics
        p1 = zeros(2*num_param,1);
        rms_p = zeros(1,sp2_itr+1);
        
        % main sp2 loop
        for i = 2:(sp2_itr+1)
            
            
            px_tmp = px;
            py_tmp = py;
            
            % calculate weights
            if  i == 2 % l2 1st iteration as initialisation
                Rm = speye(length(grad_mx));
                Rp = speye(length(px_tmp));
            else
                Rm = spdiags(1./sqrt((grad_mx - px_tmp).^2 + (grad_my - py_tmp).^2 + gamma_c^2), 0, length(grad_mx), length(grad_mx));
                Rp = spdiags(1./sqrt((cx'*px_tmp).^2 + (cy'*py_tmp).^2 + 0.5*(cx'*py_tmp + cy'*px_tmp).^2 + gamma_p.^2), 0, length(px_tmp), length(px_tmp) );
            end
            % Set up matrix equation for px, py and solve
            a11 = Rm + tgv_lagrn(j)*(cx*Rp*cx' + 0.5*cy*Rp*cy');
            a12 = tgv_lagrn(j)*0.5*cy*Rp*cx';
            a21 = tgv_lagrn(j)*0.5*cx*Rp*cy';
            a22 = Rm + tgv_lagrn(j)*(cy*Rp*cy' + 0.5*cx*Rp*cx');
            b1 = Rm*cx*log10(res_param);
            b2 = Rm*cy*log10(res_param);
            
            % Domain_ind first
            a11 = a11(domain_ind, domain_ind);
            a12 = a12(domain_ind, domain_ind);
            a21 = a21(domain_ind, domain_ind);
            a22 = a22(domain_ind, domain_ind);
            b1 = b1(domain_ind);
            b2 = b2(domain_ind);
            
            %                     A = [a11, a12; a21, a22];
            %                     b = [b1; b2];
            A = [a11, a12; a21, a22];
            b = [b1; b2];
            %                     condition = cond(A);
            %                     disp(['itr: ', num2str(i-1), ' | cond(A): ', num2str(condition,'%.2g')])
            p2 = A\b;
            %                     p2 = p2(domain_ind);
            clear A b b1 b2 a11 a12 a21 a22
            
            
            % spearate px, py
            px = p2(1:length(p2)/2);
            py = p2((length(p2)/2+1):length(p2));
            
            
            p2 = [px; py];
            
            
            
            
            rms_p(i-1) = sqrt(mean((p2 - p1).^2))/(0.5*sqrt(mean( (cx*log10(res_param)).^2 + (cy*log10(res_param)).^2)));
            %         disp(['rms dp = ',num2str(p_rms(i))])
            %                 rms_change(i) = abs((p_rms(i) - p_rms(i-1)))./(p_rms(i) + p_rms(i-1)); % percentage rms change per itr
            
            disp(['tgv iteration ',num2str(i-1),' : rms dp = ',num2str(rms_p(i-1))])
            
            
            if ((i > 4) && (rms_p(i-1) < p_conv)) %  stop if change smaller than p_conv after first few itr
                sp2_i = i;  % store number of iterations needed
                break
                %     elseif (i > 2 && (rms_p(i-1) > rms_p(i-2))) % stop if diverging
                %         sp2_i = i-1;  % store number of iterations needed
                %         mesh.px = p1(1:length(p1)/2);
                %         mesh.py = p1((length(p1)/2+1):length(p1));
                %         break
            end
            
            
            if i == (sp2_itr+1)
                sp2_i = i;
            end
            
            p1 = p2;
            
        end
        
        
        
        
        
        
        
        
        
        TGV_t1_im = (sqrt((cx*log10(res_param) - px).^2 + (cy*log10(res_param) - py).^2));
        TGV_t2_im = (sqrt( (cx'*px).^2 + (cy'*py).^2 +  0.5*(cy'*py + cx'*px).^2));
        TGV_t1(j) = sum(TGV_t1_im(:));
        TGV_t2(j) = sum(TGV_t2_im(:));
        TGV_reg(j) = TGV_t1(j) + tgv_lagrn(j).*TGV_t2(j);
        
    end
    
    clear tgv_ssim tgv_r2 tgv_r2b rms_tgv
end

if obj_flag == 1
    mkrsz = 10;
    lw = 1.2;
    fsz = 18;
    
    figure(100)
    plot(tgv_lagrn,TV_reg.*ones(1,length(tgv_lagrn)),'-','markersize',mkrsz, 'linewidth', lw)
    hold on
    plot(tgv_lagrn,TGV_reg,'o-','markersize',mkrsz, 'linewidth', lw)
    plot(tgv_lagrn,TGV_t1,'^-','markersize',mkrsz, 'linewidth', lw)
    plot(tgv_lagrn,TGV_t2,'v-','markersize',mkrsz, 'linewidth', lw)
    %     plot(tgv_lagrn,L2_reg.*ones(1,length(tgv_lagrn)),'v-','markersize',mkrsz, 'linewidth', lw)
    legend('TV = |\nabla m|_{L1}','TGV_\mu = TGV^{(1)}_{\mu} + \muTGV^{(2)}_{\mu}','TGV^{(1)}_{\mu} = |\nabla m - p|_{l1}','TGV^{(2)}_{\mu} = 0.5|\nabla p + \nabla p^T|_{l1}')% 'total TGV objective',
%     title('Impact of \mu on the terms of the TGV functional for static models','fontsize',fsz)
    xlabel('\mu','fontsize',fsz)
    ylabel(['Sum across all model cells'],'fontsize',fsz)
    set(gca,'fontsize',fsz)
    hold off
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
% RMS best images
[~, RMSmax_ind] = max(max_list(1:length(tgv_lag_list),5));
RMS_max = max_list([end-1:end, RMSmax_ind], [1, 5:7]);

disp('Best iteration: SSIM metric')
disp('    lag2      SSIM     RMS error  Iteration')
disp(S_max)
disp('Best iteration: Pearson r.^2 metric')
disp('    lag2      r.^2      RMS error  Iteration')
disp(R_max)

% Retrospectively apply RMS  threshold for convergence
%     rms_thresh = 0.9;
%     l1_itr = find(rms_l1(2:end) < rms_thresh,1); % -1 for last itr before change, +1 for 2:end
%     l2_itr = find(rms_l2(2:end) < rms_thresh,1); % -1 for last itr before change, +1 for 2:end
%     for i = 1:length(tgv_lag_list)
%         tgv_itrs(i) = find(tgv.rms{i} < rms_thresh,1) - 1; % -1 for last itr before change
%     end
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
