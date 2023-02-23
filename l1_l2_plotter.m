
clearvars
close all


inv_file = 'gb2_long2_1pc';
inv_folder = '\thesis\periodic\golden_test\l2 itr1 init\reg mesh\';
% inv_folder = 'TGV_revision\ext_mesh\line_search_dm1\';
username = getenv('username');

type = 'l1'; % Inversion type: 'l1'; l12'; 'l22'
lag = 'ls0';
append_txt = '';%
ref_flag = 0;
i_final = 2; % if empty, then will plot last itr
final_plot_flag = 1; % == 1 if want further analysis of final result
baseline_flag = 1;
baseline_folder = 'D:\TGV_revision\thesis\mesh_define\';
diff_weight_flag = 0; % Assumes fd weighted evenly in x and y, for objective caluculation of final plots
cum_sens_limit = 0.02; %used to determine sensitivity region
first_plot_itr = 2; % first iteration of parameter charts to plot (depends on initialisation: 2 is first TV iteration)

% plot properties
fsz = 14;
lw = 2;
mkrsz = 6;
linetype = 'k-o';

if ref_flag == 1; ref_text = '_ref'; else ref_text = ''; end;
% curr = pwd;
%
% file = [curr,'\',inv_folder,'\'append_txt];
file = ['D:\TGV_revision\',inv_folder,'\',inv_file,append_txt,'_',type,'_',lag,'lag',ref_text];

load(file,'final')

n_im = length(final.RMS); % number of iterations in file

if isempty(i_final)
    i_final = n_im;
end

x = unique(final.param_x);
y = unique(final.param_y);
[xgrid, ygrid] = meshgrid(x,y);

if baseline_flag == 1
    baseline = load( fullfile(baseline_folder, inv_file) );
end

% sensitivity handling
if isfield(final,'half_space_jac')
    cum_sens = mean(abs(final.half_space_jac),1);
    cum_sens_image = reshape(cum_sens,length(x),length(y))'./max(cum_sens);
    cum_sens_thresh = cum_sens_image > cum_sens_limit;
    cum_sens_thresh_param = reshape( cum_sens_thresh', [], 1);
    cum_sens_NaN = ones( size( cum_sens_image ) );
    cum_sens_NaN( ~cum_sens_thresh ) = NaN;
else
    try
        thresh_data = load('D:\TGV_revision\thesis\periodic\golden_test\pseudo init\reg mesh\strong_cool 2\TGV gb2_long2_1pc lagls0_tgv_lag10\tgv_it_1_sp1');
        cum_sens = mean(abs(thresh_data.final.half_space_jac),1);
        cum_sens_image = reshape(cum_sens,length(x),length(y))'./max(cum_sens);
        cum_sens_thresh = cum_sens_image > cum_sens_limit;
        cum_sens_thresh_param = reshape( cum_sens_thresh', [], 1);
        cum_sens_NaN = ones( size( cum_sens_image ) );
        cum_sens_NaN( ~cum_sens_thresh ) = NaN;
    catch
        cum_sens_im = ones(length(x),length(y));
    end
end


for i = 1:n_im
    
    res_image = reshape(log10(final.res_param1(:,i)),length(x),length(y))';
    
    if baseline_flag == 1
        rms_difference_tmp =  ( res_image - log10(baseline.ip4di_direct) ).^2;
        rms_difference( i ) = sqrt( sum( sum( rms_difference_tmp( cum_sens_thresh ) ) ) / sum ( cum_sens_thresh_param ) );
        
        % calculate correlation coefficient
        r_tmp = corrcoef(log10(baseline.ip4di_direct(:)), res_image(:) );
        r( i ) = r_tmp( 1, 2 ); % above function creates matrix solution
        ssimval( i ) = ssim( log10 ( baseline.ip4di_direct ), res_image );
    end
    
    figure(21)
    subplot( ceil( n_im /2 ), 2, i )
    surf(x,y,res_image,'edgecolor','none')
    view([0,0,1])
    title(['log(m), itr = ',num2str(i),', RMS = ',num2str(final.RMS(i))],'interpreter','none')
    set(gca,'ydir','reverse')
    colorbar
    axis image
    colormap parula
end





% individual plot of final model, including objective function components
if final_plot_flag == 1
    
    sum_factor = ones( size( final.param_x ) ); % can be used to isolate only part of the model
    res_image = reshape( log10( final.res_param1( :, i_final ) ), length(x), length(y) )';
    

    if baseline_flag == 1

        misfit_text = num2str( rms_difference( i_final ) );
        
        figure(24)
        surf(x,y, log10(baseline.ip4di_direct).*cum_sens_NaN ,'edgecolor','none')
        view([0,0,1])
        title('baseline','interpreter','none')
        set(gca,'ydir','reverse')
        colorbar
        axis image
        colormap parula
        
        figure(25)
        surf(x,y, abs( log10(baseline.ip4di_direct) - res_image) .* cum_sens_NaN,'edgecolor','none')
        view([0,0,1])
        title('|baseline - solution|','interpreter','none')
        set(gca,'ydir','reverse')
        colorbar
        axis image
        colormap parula
        
    else
        misfit_text = '?';
    end
    
    [cx, cy] = c_calc_periodic( unique(final.param_x), unique(final.param_y), diff_weight_flag);
    
    l1_obj = sqrt( sum( cum_sens_thresh_param .* sum_factor .* ( cx * final.res_param1(:, i_final) ).^2 + ( cy * final.res_param1(:, i_final) ).^2 ) / sum( cum_sens_thresh_param ) );
    
    figure(23)
    surf(x,y,res_image .* cum_sens_NaN,'edgecolor','none')
    view([0,0,1])
    title(['log(m), itr = ', num2str( i_final ),', RMS Data Misfit = ', num2str( final.RMS(i_final) ), '; Model Misfit = ', misfit_text ],'interpreter','none') %  regularisation term = ', num2str( l1_obj ), 
    set( gca,'ydir','reverse', 'fontsize', fsz )
    colorbar
    axis image
    colormap parula
end

if baseline_flag == 1
    
%     figure(26)
%     scatter(log10( baseline.ip4di_direct( : ) ), res_image( : ) )
%     hold on
%     plot([1.5, 4], [1.5, 4])
%     axis image
%     xlabel('Synthetic model cell resistivity')
%     ylabel('Inverted model cell resistivity')
%     title( [ 'Pearson correlation coefficient r = ', num2str( r( i_final ) ), '; SSIM = ', num2str( ssimval( i_final ) ) ] )
%     hold off


    plot_ind = first_plot_itr:length(final.RMS);
    itr_ind = 0:length(plot_ind) - 1;
    
    figure(27)
    subplot( 2, 2, 4 )
    hold on
    semilogy( itr_ind, ssimval( plot_ind ), linetype, 'linewidth', lw, 'markersize', mkrsz )
    title( 'ssim' )
    ylabel( 'ssim' )
    xlabel( 'iteration' )
    set( gca, 'fontsize', fsz, 'Yscale', 'log' ) 
    
    figure( 27 )
    subplot( 2, 2, 3 )
    hold on
    semilogy( itr_ind, rms_difference( plot_ind ), linetype, 'linewidth', lw, 'markersize', mkrsz )
    title( 'rms model misfit' )
    ylabel( 'rms model misfit' )
    xlabel( 'iteration' )
    set( gca, 'fontsize', fsz, 'Yscale', 'log' ) 
    
    figure( 22 )
    semilogy( itr_ind, r( plot_ind ), linetype, 'linewidth', lw, 'markersize', mkrsz )
    title( 'pearson r' )
    ylabel( 'pearson r' )
    xlabel( 'iteration' )
    set( gca, 'fontsize', fsz, 'Yscale', 'log' ) 
    
    figure( 27 )
    subplot( 2, 2, 2 )
    hold on
    semilogy( itr_ind, final.RMS( plot_ind ), linetype, 'linewidth', lw, 'markersize', mkrsz )
    title( 'rms data misfit' )
    ylabel( 'rms data misfit' )
    xlabel( 'iteration' )
    set( gca, 'fontsize', fsz, 'Yscale', 'log' ) 

    
end


figure(27)
subplot( 2, 2, 1 )
hold on
semilogy( itr_ind( 2 : end ),  final.lagrn( plot_ind( 2 : end ) ), linetype, 'linewidth', lw, 'markersize', mkrsz )
ylabel('\lambda')
xlabel( 'iteration' )
title('\lambda')
set( gca, 'fontsize', fsz , 'Yscale', 'log') 



figure(21)




