%% Plot results for different iterations
clearvars
% close all

curr_dir = pwd;
file = 'gb2_long2_1pc';
inv_folder = '\thesis\periodic\golden_test\tgv l2 solution init\reg mesh\wd repeat\strong_cool 2';
% folder = 'TGV_revision\ext_mesh\central_diff\dm1\strong_cool';
user = getenv('username');
% mu_in = [{'7.5'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}, {'14'}, {'15'}, {'17.5'}, {'20'}, {'20'}, {'25'}, {'30'} ];    % As appears in file name (*10, string)
mu_in = [ {'8'}, {'10'}, {'11.5'},  {'13'},  {'15'}, {'16.5'}, {'18'}, {'20'}, {'22.5'}, {'25'}, {'27.5'}, {'30'} ];    % As appears in file name (*10, string)
% mu_in = [{'15'}];
append_txt = '';
lagrn_txt = 'ls0';
cum_sens_limit = 0.02;
diff_weight_flag = 0; % Assumes fd weighted evenly in x and y, for objective caluculation of final plots
n_max = 99; % maximum iteration number - used for subplots
init_model_flag = 0; % = 1 if first iteration is a loaded in model (as no data from that itr)
best_flag = 2; % Determines which model is used for parameters. [0] for last (converged model) [1] for lowest model misfit [2] for lowest data misfit

% plot properties
first_plot_itr = 2;
fsz = 14;
lw = 2;
mkrsz = 6;
linetype = 'k-o';


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

n_mu = length(mu_in);
i_count = 0;

for j = 1:n_mu
    
    
    tgv_lagr = mu_in{ j };
    pathname = ['D:\TGV_revision\',inv_folder,'\TGV ',file,append_txt,' lag',lagrn_txt,'_tgv_lag',tgv_lagr];

    for i = n
        
        rms_itr = 100;
        lagrn_itr = 100;
        
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
%             if init_model_flag == 1
%                 n = [1, 3:i-1];
%             else
%                 n = 1:i-1;
%             end
             n_im = i;
            
            break
        end
        cd(curr_dir);
        
        % Extract info
        if lagrn_flag == 1; lagrn(end+1) = final.lagrn(end); lagrn_itr = min([lagrn_itr, 0]); end
        if rms_flag == 1; rms(end+1) = final.RMS(end); rms_itr = min(rms_itr, 0); end
        
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
            l1_obj( i ) = sum( cum_sens_thresh_param .* sum_factor .* sqrt( ( cx * log10(final.res_param1(:, end))).^2 + ( cy * log10(final.res_param1(:, end))).^2 )  ) / sum( cum_sens_thresh_param );
            tgv1_obj( i ) = sum( cum_sens_thresh_param .* sum_factor .* sqrt(  ( cx * log10(final.res_param1(:, end)) - final.px ).^2 + ( cy * log10(final.res_param1(:, end)) - final.py ).^2 )  ) / sum( cum_sens_thresh_param );
            tgv2_obj( i ) = sum( cum_sens_thresh_param.*sum_factor .* sqrt( ( cx' * final.px ).^2 + ( cy' * final.py ).^2 + 0.5 * ( cy' * final.px + cx' * final.py).^2 ) ) / sum( cum_sens_thresh_param );
            p_sum( i ) = sum( cum_sens_thresh_param.*sum_factor .* sqrt( (final.px ).^2 + ( final.py ).^2 ) ) / sum( cum_sens_thresh_param );
        else % adds dummy values for easier interpretation
            [cx, cy] = c_calc_periodic( unique(final.param_x), unique(final.param_y), diff_weight_flag);
            l1_obj( i ) = 0;
            tgv1_obj( i ) = 0;
            tgv2_obj( i ) = 0;
            p_sum( i ) = 0;
        end
        
        % calculates model misfit
        % sensitivity handling
        
        rms_difference_tmp =  ( res_image - log10(baseline.ip4di_direct) ).^2;
        rms_difference( i ) = sqrt( sum( sum( rms_difference_tmp( cum_sens_thresh ) ) ) / sum ( cum_sens_thresh_param ) );
%         rms2( i ) = final.RMS( end );
        misfit_text = num2str(rms_difference(i));
        
        % calculate correlation coefficient
        r_tmp = corrcoef(log10(baseline.ip4di_direct(:)), res_image(:) );
        r( i ) = r_tmp( 1, 2 ); % above function creates matrix solution
        ssimval( i ) = ssim( log10 ( baseline.ip4di_direct ), res_image );
        i_count( end + 1 ) = i;
        
        if i == 1
            mu( j ) = final.tgv_lagrn; 
        end
        
        
    end
    
    if lagrn_itr == 100; lagrn_itr = 0; end
    if rms_itr == 100; rms_itr = 0; end
    
        % save tgv objective parameters for the best model
        if best_flag == 1
            [model_misfit( j ), solution_ind] = min( rms_difference(3:end) );          
            solution_ind = solution_ind + 2; % corrects for truncation
        elseif best_flag == 2
            [~, solution_ind] = min( rms(4:end) );
            solution_ind = solution_ind + 3; % corrects for truncation
            model_misfit(j) = rms_difference(solution_ind);
        else
            model_misfit(j) = rms_difference( end );
            solution_ind = length( rms_difference );
        end
        
        itr_index( j ) = solution_ind;
        data_misfit( j ) = rms( solution_ind - init_model_flag ); % -1 if itr 2 missing
        r_model( j ) = r( solution_ind );
        tgv1_obj_all( j ) = tgv1_obj( solution_ind );
        tgv2_obj_all( j ) = tgv2_obj( solution_ind );
        l1_obj_all( j ) = l1_obj( solution_ind );
        p_sum_all( j ) = p_sum( solution_ind );
        
        clear itr_index tgv1_obj tgv2_obj l1_obj p_sum r rms_difference
        rms = [];
end


plot_ind = first_plot_itr:length(final.RMS);
itr_ind = 0:length(plot_ind) - 1;



figure
subplot(1, 2 ,1)
plot( mu, l1_obj_all, 'k-o', 'markersize', mkrsz, 'linewidth', lw )
hold on
plot(mu, tgv1_obj_all, 'b-^', 'markersize', mkrsz, 'linewidth', lw)
plot(mu, tgv2_obj_all, 'r-v', 'markersize', mkrsz, 'linewidth', lw)
plot(mu, tgv1_obj_all + final.tgv_lagrn*tgv2_obj_all, 'm-d', 'markersize', mkrsz, 'linewidth', lw )
% plot(mu, p_sum_all, 'c-s', 'markersize', mkrsz, 'linewidth', lw )
legend('TV = |\nabla m|_{L1}','TGV_\mu = TGV^{(1)}_{\mu} + \muTGV^{(2)}_{\mu}','TGV^{(1)}_{\mu} = |\nabla m - p|_{l1}','TGV^{(2)}_{\mu} = 0.5|\nabla p + \nabla p^T|_{l1}', 'location', 'best')% 'total TGV objective',
title('Objective terms')
set(gca, 'fontsize', fsz)
ylim([0, 0.32])

% legend('TV = |\nabla m|_{L1}','TGV_\mu = TGV^{(1)}_{\mu} + \muTGV^{(2)}_{\mu}','TGV^{(1)}_{\mu} = |\nabla m - p|_{l1}','TGV^{(2)}_{\mu} = 0.5|\nabla p + \nabla p^T|_{l1}')% 'total TGV objective',

subplot( 2,2,2)
plot(mu, data_misfit, 'kd-','markersize', mkrsz, 'linewidth', lw )
title('RMS data misfit')
set(gca, 'fontsize', fsz, 'Yscale', 'log')
yticks([1.1:0.2:1.7, 2:0.5:3])
ylim([0.9, 2.2])
subplot(2,2,4)
plot(mu, model_misfit, 'kd-', 'markersize', mkrsz, 'linewidth', lw )
title('RMS model misfit')
set(gca, 'fontsize', fsz)
ylim([0.18,0.3])

