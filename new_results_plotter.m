% function new_results_plotter(filename, mu, lag, itr, baseline_flag)
% plots tgv results (final iterations)
% filename -     subfolder/datafilename of inversion results
% mu -   cell of mu value strings. 6 max (*10 as in filename convention)
% itr -     cell array giving the iteration to plot of each result. [] entries give
%           the final itr. If a single value, then that value will be used for all
%           entires. Large values default to [].
% lag -     cell array giving the damping parameters used. Can be scalar if the
%           same lambda was used for all inversions
% both vector orders go l1, l2, mu(1), mu(2), ... , mu(end)
% baseline_flag -   True if a baseillne image exists, false otherwise
function new_results_plotter(filename, mu, lag, itr, baseline_flag)
close all

% Control parameters - occasionally useful to change
pathname = ['C:\users\',getenv('username'),'\OneDrive - The University of Nottingham\UoN Box Migration\TGV_revision\'];
cum_sens_limit = 0.03;
append_txt = ''; % for old filenames, try not to use!




% preprocessing
n_im = length(mu) + 2;
[sub_path, filename, file_ext] = fileparts(filename);
% expand scalar iputs
if length(itr) == 1
    itr = repmat(itr, n_im, 1);
elseif isempty(itr)
    itr = repmat({[]}, n_im, 1);
end
if length(lag) == 1
    lag = repmat(lag, n_im, 1);
end
climits = [1e15, 0]; % preallocation; always overwritten
% plimits = [1e15, 0]; % preallocation; always overwritten
% open figures in order
for i = 1:ceil((n_im+1)/3)
    figure(i)
end


% load in baseline
if baseline_flag == 1
    try
        load(['C:\Users\',getenv('username'),'\OneDrive - The University of Nottingham\UoN Box Migration\TGV_revision\Models\',filename,'.mat']);
    catch
        disp('**** Baseline not found ***')
    end
    % open figures in order
    for i = 3 + (1:ceil((n_im+1)/3))
        figure(i)
    end
    % read in baseline
    baseline = log10(ip4di_im);
    % preallocate
    ssim_metric = zeros(n_im, 1);
    r_metric = zeros(n_im, 1);
    rms_error_metric = zeros(n_im, 1);
    climits_ratio = [1e15, 0];
    % need to do some work to show non-interpolated baseline
    %     baseline = log10(thresh_im);
end
% load resinv colormap
load(['C:\Users\',getenv('username'),'\OneDrive - The University of Nottingham\UoN Box Migration\Coupled TGV results\resInvColormap'])

%% load and plot files
for i = 1:n_im
    
    % load in file
    if i == 1 %l1
        filetype = [filename, '_l1_', lag{i}, 'lag', append_txt];
    elseif i == 2 %l2
        filetype = [filename, '_l2_', lag{i}, 'lag', append_txt];
    else
        % find tgv iteration
        if isempty(itr{i})
            for j = 15:-1:1
                sp1_file = fullfile([pathname, sub_path, '\TGV ', filename,' lag',lag{i},'_tgv_lag', mu{i-2},append_txt,'\tgv_it_', num2str(j),'_sp1.mat']);
                if exist(sp1_file,'file') == 2
                    itr{i} = j;
                    break
                end
            end
        end
        
        if j~= 1
            sp2_file = fullfile([pathname, sub_path, '\TGV ', filename,' lag',lag{i},'_tgv_lag', mu{i-2},append_txt,'\tgv_it_', num2str(itr{i} - 2),'_sp2.mat']);
        end
        filetype = ['/TGV ', filename, ' lag', lag{i}, '_tgv_lag', mu{i-2}, append_txt, '\tgv_it_', num2str(itr{i}),'_sp1'];
    end
    
    
    file = fullfile(pathname, sub_path, filetype);
    
    load(file, 'final')
    % set itr to final itr if not already specified
    if isempty(itr{i})
        itr{i} = final.itn;
    end
    
    if i == 1
        % Calculate x and y dimensions of mesh (in pixels).
        [xgrid, ygrid] = meshgrid(unique(final.param_x), unique(final.param_y));
        xfind = find(final.param_x == final.param_x(1),2);
        len_xdim = xfind(2) - xfind(1);
        len_ydim = final.num_param./len_xdim;
        % Translate so that plot cells are center on parameter values, instead of
        % holding those values on nodes
        x = unique(final.param_x); x = (x + [0; x(1:end-1)])/2;
        y = unique(final.param_y); y = (y + [0; y(1:end-1)])/2;
        
        % calculate sensitivities
        if isfield(final,'half_space_jac')
            if isfield(final,'half_space_jac')
                cum_sens = mean(abs(final.half_space_jac),1);
            end
            cum_sens_image = reshape(cum_sens,len_xdim,len_ydim)'./max(cum_sens);
            sens_mask = cum_sens_image > cum_sens_limit;
                        
            figure(20)
            surf(x,y,cum_sens_image.*(cum_sens_image>cum_sens_limit),'edgecolor','none')
            view([0,0,1])
            title(['Cumulative sensitivity > ',num2str(cum_sens_limit)],'interpreter','none')
            set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)])
            colorbar
            colormap jet
            axis image
        else
            cum_sens_im = ones(len_xdim,len_ydim);
        end
    end
    
    if i < 3
        res_image = reshape(log10(final.res_param1(:,itr{i})),len_xdim,len_ydim)';
        rms = final.RMS(itr{i});
    else
        res_image = reshape(log10(final.res_param1(:,end)),len_xdim,len_ydim)';
        rms = final.RMS(end);
    end
    
    % Plot resistivty.
    if i < 3
        fignum = 1;
        plotnum = i + 1;
    elseif i < 6
        fignum = 2;
        plotnum = i - 2;
    elseif i < 9
        fignum = 3;
        plotnum = i - 5;
    end
    
    climits = [min([min(res_image(:)), climits(1)]), max([max(res_image(:)), climits(2)])];
    
    figure(fignum)
    subplot(3,1,plotnum)
    surf(x,y,res_image,'edgecolor','none')
    view([0,0,1])
    title([filetype,';   RMS = ', num2str(rms)],'interpreter','none')
    set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)])
    colorbar
    axis image
    
    % plot p
    if i > 2  &&  itr{i} > 1
        fignum = fignum + 2;
        sp2 = load(sp2_file, 'final');
        
        px_im = reshape(sp2.final.px,length(x),length(y))';
        py_im = reshape(sp2.final.py,length(x),length(y))';
        pz_im = px_im.^2 + py_im.^2;
        z_cen = 100*ones(size(px_im));
        z_im = zeros(size(px_im));
%                     plimits = [min([min(pz_im(:)), plimits(1)]), max([max(pz_im(:)), plimits(2)])];
        if ~exist('plimits', 'var')
            plimits = [0, max(pz_im*1.5)];
        end
        
        figure(fignum)
        subplot(3, 1, plotnum)
        surf(x,y,sqrt(pz_im),'edgecolor','none')
        view([0,0,1])
        title(['p: ', sp2_file],'interpreter','none')
        colorbar
        set(gca,'ydir','reverse')
        caxis([min(min(sqrt(px_im.^2 + py_im.^2))),max(max(sqrt(px_im.^2 + py_im.^2)))]);
        hold on
        quiver3(unique(final.param_x), unique(final.param_y), z_cen, px_im, py_im, z_im,'color','m');
        axis image
        
    end
    
    
        
        
    % create and plot the baseline image and calculate fidelity metrics
    if baseline_flag ~= 0
        ratio_image = res_image - baseline;
        climits_ratio = [min(climits_ratio(1), min(abs(ratio_image(:)))), max(climits_ratio(2), max(abs(ratio_image(:))))];
        fignum = fignum + 5;
        
        sens_bg = nonzeros(baseline.*sens_mask);

        
        pc_image = nonzeros((ratio_image./res_image).*sens_mask);
        sens_image = nonzeros(res_image.*sens_mask); % extracts high sensitvity region
        [ssim_image, ~] = ssim(res_image, baseline, 'exponents', [1, 1, 1]);
        
        % rms percentage difference
        rms_error_metric(i) = 100*sqrt(sum((pc_image).^2))/length(pc_image);
        % pearson correlation coefficiant r^2
        r_metric(i) = sum(sens_image - mean(sens_image).*(sens_bg - mean(sens_bg)))./sqrt(sum(sens_image - mean(sens_image)).^2).*sqrt(sum(sens_bg - mean(sens_bg)).^2);
        ssim_metric(i) = sum(sum(ssim_image.*sens_mask))./sum(sens_mask(:)); 
        
        figure(fignum)
        subplot(3, 1, plotnum)
        surf(x,y,abs(ratio_image),'edgecolor','none')
        view([0,0,1])
        title([file, ' RATIO'],'interpreter','none')
        set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)])
        colorbar
        axis image
%         colormap summer
        
    end
    
end

if baseline_flag ~= 0
    % plot baseline image in subplot 1 fig 1
    climits = [min([min(baseline(:)), climits(1)]), max([max(baseline(:)), climits(2)])];
    figure(1)
    subplot(3,1,1)
    surf(x,y,baseline,'edgecolor','none')
    view([0,0,1])
    title('Original resistivity distribution','interpreter','none')
    set(gca,'ydir','reverse','xlim',[min(x), max(x)],'ylim',[min(y), max(y)])
    colorbar
    axis image
    set(gca, 'clim', climits)
    
    % create a table of metrics and then print it to cmd
    regularisation_method = [{'l1'}, {'l2'}, mu]';
    metrics = table(regularisation_method, rms_error_metric, r_metric, ssim_metric)

end




% update color axis to be uniform across all images
% plimits = [plimits(1), 1.5*plimits(2)]; % helps avoid saturation at top of color scale
plotnum = 2; % starts with l1 plot
for fignum = 1:ceil((n_im+1)/3)
    figure(fignum)
    while plotnum < 4
        subplot(3, 1, plotnum)
        set(gca, 'clim', climits)
        
        % p plot caxis
%         if fignum > 1
%             figure(fignum + 2)
%             subplot(3, 1, plotnum)
%             caxis(plimits)
%         end
%         
        plotnum = plotnum + 1;
        if (plotnum + 3*(fignum-1) - 1) > n_im% break after last image
            break
        end
    end
    plotnum = 1;
    
    if baseline_flag ~= 0
        figure(fignum + 5)
%         if fignum == 1; plotnum = 2; end
        while plotnum < 4
            subplot(3, 1, plotnum)
            set(gca, 'clim', climits_ratio)
            plotnum = plotnum + 1;
            if (plotnum + 3*(fignum-1) - 1) > n_im
                break
            end
        end
        plotnum = 1;
    end
end

if baseline_flag == 1
    disp(['Cumulative sensitivity inclusion threshold for metrics = ', num2str(cum_sens_limit)]);
end

figure(1)

end