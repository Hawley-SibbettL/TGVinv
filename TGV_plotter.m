%% Plot results for different iterations
clearvars
close all

curr_dir = pwd;
file = 'gb2_long2_1pc';
inv_folder = '\thesis\periodic\golden_test\tgv pseudo solution init\reg mesh\wd repeat\strong_cool 2';
% folder = 'TGV_revision\ext_mesh\central_diff\dm1\strong_cool';
user = getenv('username');      
tgv_lagr = '15';    % As appears in file name (*10, string)
append_txt = '';
lagrn = 'ls0';
cum_sens_limit = 0.02;
diff_weight_flag = 0; % Assumes fd weighted evenly in x and y, for objective caluculation of final plots
n_max = 99; % maximum iteration number - used for subplots
init_model_flag = 0; % = 1 if first iteration is a loaded in model (as no data from that itr)

% plot properties
first_plot_itr = 3;
fsz = 14;
lw = 2;
mkrsz = 6;
linetype = 'k-o';

%
pathname = ['D:\TGV_revision\',inv_folder,'\TGV ',file,append_txt,' lag',lagrn,'_tgv_lag',tgv_lagr];
% curr = pwd;
% pathname = [curr,'\',inv_folder,'\TGV ',file,' lag',lagrn,'_tgv_lag',tgv_lagr,append_txt];

% if get errors with 'final.*' 
if init_model_flag == 1
    n = [1, 3:n_max];
else 
    n = 1:n_max;
end
n_im = length(n);


baseline_flag = 1; % = 1 for baseline
baseline_folder =  'D:\TGV_revision\thesis\mesh_define\';
p_plot_itr = [1:40]; % itr to plot p (if present)
lagrn_flag = 1;
sp2_itr_flag = 1;
rms_p_flag = 0; % Show convergence of p 
rms_flag = 1; % show rms misfit
fin_flag = 0; % Finningley cross sections
ratio_plot_flag = 0; % ratio plots (if baseline present)

lagrn = [];
rms = [];
sp2_itr = [];

% if lagrn_flag == 1; lagrn =0; end
% if sp2_itr_flag == 1; sp2_t = 0; end
sp1_flag = 0;

if baseline_flag == 1
    baseline = load( fullfile(baseline_folder, file) );
end




for i = n
    
    cd(pathname)
    if exist(['tgv_it_',num2str(i),'_sp2.mat'],'file') == 2
        load(['tgv_it_',num2str(i),'_sp2.mat']);
        if sp2_itr_flag == 1; sp2_itr(end+1) = final.sp2_itr(end); end
    elseif exist(['tgv_it_',num2str(i),'_sp1.mat'],'file') == 2 
        load(['tgv_it_',num2str(i),'_sp1.mat']);
        sp1_flag = 1;
    else
         % if get errors with 'final.*' then break
        cd(curr_dir)
        disp(['Ending: filename ',num2str(i) ,' not found'])

        % readjust number of models 
        if init_model_flag == 1
            n = [1, 3:i-1];
        else
            n = 1:i-1;
        end
        n_im = length(n);
        
        break
    end
    cd(curr_dir);
    

    
    
    % Extract info
    if lagrn_flag == 1; lagrn(end+1) = final.lagrn(end); end
    if rms_flag == 1; rms(end+1) = final.RMS(end); end
    
    x = unique(final.param_x);
    y = unique(final.param_y);
    [xgrid, ygrid] = meshgrid(x,y);
    res_image = reshape(log10(final.res_param1(:,end)),length(x),length(y))';
    cum_sens_im = ones(length(x),length(y));    % used for sensitivity modelling later
    sum_factor = ones(size(final.param_x)); % can be used to isolate only part of the model

    if isfield(final,'half_space_jac')
        cum_sens = mean(abs(final.half_space_jac),1);
        cum_sens_image = reshape(cum_sens,length(x),length(y))'./max(cum_sens);
        cum_sens_thresh = cum_sens_image > cum_sens_limit;
        cum_sens_thresh_param = reshape( cum_sens_thresh', [], 1);
        cum_sens_NaN = ones( size( cum_sens_image ) );
        cum_sens_NaN( ~cum_sens_thresh ) = NaN;
        cum_sens_NaN2 = reshape(cum_sens_NaN', [], 1);
    end
    
    % calculate objective terms
    if isfield(final, 'px')
        [cx, cy] = c_calc_periodic( unique(final.param_x), unique(final.param_y), diff_weight_flag);
        l1_obj = sum( cum_sens_thresh_param .* sum_factor .* sqrt(  ( cx * log10(final.res_param1(:, end))).^2 + ( cy * log10(final.res_param1(:, end))).^2 )  ) / sum( cum_sens_thresh_param );;
        tgv1_obj( i ) = sum( cum_sens_thresh_param .* sum_factor .* sqrt(  ( cx * log10(final.res_param1(:, end)) - final.px ).^2 + ( cy * log10(final.res_param1(:, end)) - final.py ).^2 )  ) / sum( cum_sens_thresh_param );
        tgv2_obj( i ) = sum( cum_sens_thresh_param.*sum_factor .* sqrt( ( cx' * final.px ).^2 + ( cy' * final.py ).^2 + 0.5 * ( cy' * final.px + cx' * final.py).^2 ) ) / sum( cum_sens_thresh_param );
    end
    
    % calculates model misfit
    if baseline_flag == 1
        % sensitivity handling

        rms_difference_tmp =  ( res_image - log10(baseline.ip4di_direct) ).^2;
        rms_difference( i ) = sqrt( sum( sum( rms_difference_tmp( cum_sens_thresh ) ) ) / sum ( cum_sens_thresh_param ) );
        misfit_text = num2str(rms_difference(i), '%.3f');
        
        % calculate correlation coefficient
        r_tmp = corrcoef(log10(baseline.ip4di_direct(:)), res_image(:) );
        r( i ) = r_tmp( 1, 2 ); % above function creates matrix solution
        ssimval( i ) = ssim( log10 ( baseline.ip4di_direct ), res_image );
        

    else
        misfit_text = '?';
    end
    
    
    
%     figure(20)
%     subplot(ceil(n_im/2),2,find(n==i))
    figure(i)
    surf(x,y, res_image .* cum_sens_NaN,'edgecolor','none')
    view([0,0,1])
    title(['RMS data misfit = ',num2str(final.RMS(end), '%.2f'),'; Model Misfit = ', misfit_text],'interpreter','none')
    set(gca,'ydir','reverse','fontsize',16)
    colorbar
    axis image
    colormap parula
    
%     if i == 1
%         figure(20)
%         subplot(ceil(n_im/2),2,1)
%         surf(x,y,res_image,'edgecolor','none')
%         view([0,0,1])
%         title(['log(m), itr = ','0',', RMS = ',num2str(final.RMS(1))],'interpreter','none')
%         set(gca,'ydir','reverse')
%         colorbar
%         axis image
%         colormap parula
%     end
    
    
    if baseline_flag == 1 && ratio_plot_flag == 1
        figure(136)
        subplot(ceil(n_im/2),2,find(n==i))
        surf(x,y,abs(res_image - log10(baseline.ip4di_direct) ) .* cum_sens_NaN,'edgecolor','none')
        view([0,0,1])
        title(['log(m), itr = ',num2str(i),', RMS = ',num2str(final.RMS(end)),' RATIO'],'interpreter','none')
        set(gca,'ydir','reverse')
        colorbar
        axis image
    end
    
    % p plots
    
    if sp1_flag == 0
        
        %     load(['tgv_it_',num2str(i),'_sp2'],'mesh')
        px_im = reshape(final.px,length(x),length(y))';
        py_im = reshape(final.py,length(x),length(y))';
        z_cen = 100*ones(size(px_im));
        pz_im = zeros(size(px_im));
        
        
        if isempty( find(p_plot_itr == i ) ) == 0
%             figure(137)
%             subplot(ceil(n_im/2),2,find(n==i))
            figure( i + 200 )
%             surf(x,y,sqrt(px_im.^2 + py_im.^2),'edgecolor','none')
            surf(x,y,sqrt(px_im.^2 + py_im.^2).*cum_sens_NaN,'edgecolor','none')
            view([0,0,1])
%             title(['p, itr = ',num2str(i)],'interpreter','none')
            colorbar
            set(gca,'ydir','reverse', 'fontsize', fsz)
            caxis([min(min(sqrt(px_im.^2 + py_im.^2))),max(max(sqrt(px_im.^2 + py_im.^2)))]);
            axis image
            hold on
%             quiver3(unique(final.param_x), unique(final.param_y), z_cen, px_im, py_im, pz_im,'color','m');
            quiver3( final.param_x, final.param_y, 100*ones(size(final.px)), final.px, final.py, zeros( size( final.px ) ) ,'color','m');
%             quiver3( cum_sens_NaN2.*final.param_x, cum_sens_NaN2.*final.param_y, 100*cum_sens_NaN2.*ones(size(final.px)), cum_sens_NaN2.*final.px, cum_sens_NaN2.*final.py, cum_sens_NaN2.*zeros( size( final.px ) ) ,'color','m');
        end
        
        p_im_final = sqrt(px_im.^2 + py_im.^2);
        
        pdiff = sqrt( ( final.px - cx * log10(final.res_param1(:, end)) ).^2 + (final.py - cy * log10(final.res_param1(:, end)) ).^2);
        pdiff_im = reshape( pdiff, length(x), length(y) )';
        
        figure(i + 300)
        surf(x,y,pdiff_im.*cum_sens_NaN,'edgecolor','none')
        view([0,0,1])
%                     title(['p, itr = ',num2str(i)],'interpreter','none')
        colorbar
        set(gca,'ydir','reverse', 'fontsize', fsz)
%         caxis([min(min(sqrt(px_im.^2 + py_im.^2))),max(max(sqrt(px_im.^2 + py_im.^2)))]);
        axis image
        
%         figure(18)
%         subplot(ceil(n_im/2),2,find(n==i))
%         surf(x,y,1*(sqrt(px_im.^2 + py_im.^2)>0.01),'edgecolor','none')
%         view([0,0,1])
%         title(['p, itr = ',num2str(i)],'interpreter','none')
%         set(gca,'ydir','reverse')
%         colorbar
%         colormap cool
%         axis image
        
        
        % p convergence plot
        if rms_p_flag == 1
            figure(135)
            subplot(ceil(n_im/2),2,find(n==i))
            plot(final.rms_p)
            xlabel('iteration')
            ylabel('rms dp/p')
            set( gca, 'fontsize', fsz) 
        end

        
    end
    
end

plot_ind = first_plot_itr:length(final.RMS);
itr_ind = 0:length(plot_ind) - 1;


% if baseline_flag == 1
%     figure(17)
%     surf(x,y,abs(baseline.ip4di_im),'edgecolor','none')
%     view([0,0,1])
%     title('baseline','interpreter','none')
%     set(gca,'ydir','reverse')
%     axis image
%     colorbar
% end

if lagrn_flag == 1
    
    figure(148)
    subplot( 2, 2, 1 )
    hold on
    semilogy( itr_ind( 2:end ), lagrn( plot_ind( 2:end )), linetype, 'linewidth', lw, 'markersize', mkrsz)
    set( gca, 'Yscale', 'log', 'fontsize', fsz )
    title('\lambda')

end

if sp2_itr_flag == 1
    figure(147)
    subplot(2,1,2)
    plot(n(1:length(sp2_itr)),sp2_itr,'-^')
    title('subproblem 2 iterations')
end

if rms_flag == 1
    
    figure(148)
    subplot(2,2,2)
    hold on
    semilogy( itr_ind , rms( plot_ind ) , linetype, 'linewidth', lw, 'markersize', mkrsz )
    title('rms data misfit')
    set( gca, 'fontsize', fsz , 'Yscale', 'log') 

    
    if baseline_flag == 1
    
        figure(148)        
        subplot( 2, 2, 3 )
        hold on
        semilogy( itr_ind, rms_difference( plot_ind ) , linetype, 'linewidth', lw, 'markersize', mkrsz )
        title( 'rms model misfit' )
        ylabel( 'rms model misfit' )
        xlabel( 'iteration' )
        set( gca, 'fontsize', fsz , 'Yscale', 'log') 

        figure(147)
        subplot( 2, 1, 1 )
        hold on
        semilogy( itr_ind, r( plot_ind ) , linetype, 'linewidth', lw, 'markersize', mkrsz )
        title( 'pearson r' )
        ylabel( 'pearson r' )
        xlabel( 'iteration' )
        set( gca, 'fontsize', fsz , 'Yscale', 'log') 

        figure(148)
        subplot( 2, 2, 4 )
        hold on
        semilogy( itr_ind, ssimval( plot_ind ), linetype, 'linewidth', lw, 'markersize', mkrsz )
        title( 'ssim' )
        ylabel( 'ssim' )
        xlabel( 'iteration' )
        set( gca, 'fontsize', fsz , 'Yscale', 'log') 

    end
    
end

if fin_flag == 1
    for x_slice = [23, 65]
        figure
        plot(res_image(:,x_slice), y,'g--^')
        set(gca,'ydir','reverse', 'fontsize', fsz )
        xlabel('log_{10}\bf{m}')
        ylabel('z (m)')
        title(['section at x = ',num2str(x(x_slice))])
    end
end

figure(49)
plot(  [0, length(final.RMS)], repmat( l1_obj, 2, 1) )
hold on
plot(tgv1_obj)
plot(tgv2_obj)
plot(tgv1_obj + final.tgv_lagrn*tgv2_obj)
legend('TV = |\nabla m|_{L1}','TGV_\mu = TGV^{(1)}_{\mu} + \muTGV^{(2)}_{\mu}','TGV^{(1)}_{\mu} = |\nabla m - p|_{l1}','TGV^{(2)}_{\mu} = 0.5|\nabla p + \nabla p^T|_{l1}')% 'total TGV objective',
% legend('TV = |\nabla m|_{L1}','TGV_\mu = TGV^{(1)}_{\mu} + \muTGV^{(2)}_{\mu}','TGV^{(1)}_{\mu} = |\nabla m - p|_{l1}','TGV^{(2)}_{\mu} = 0.5|\nabla p + \nabla p^T|_{l1}')% 'total TGV objective',
