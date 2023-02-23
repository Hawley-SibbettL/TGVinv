function [] =  p_init_analysis()


% Read in model
folder = 'D:\TGV_revision\thesis\mesh_define\gb1_refgrid'; % contains all p calculations for the named synthetic model
p_init = [2, 1, 3];
mu = [0.4, 1, 1.5, 2];

% load p
[p_dat, model, cx, cy, xdim] = load_p(folder, p_init, mu);

% calculate analysis measures
if length(p_init) > 1
    p_dist = p_analysis(p_dat, model, cx, cy, xdim, mu, p_init);
    
    p_header = {'rms difference', '0.4', '1.0', '1.5', '2.0'};
    p_rows = {'l2 - grad(m)'; 'l2 - gaussian'; 'grad(m) - gaussian'};
    p_dist_table = table(p_rows, p_dist(:, 1), p_dist(:,2), p_dist(:, 3), p_dist(:, 4));
    p_dist_table.Properties.VariableNames = p_header;
    
end

% find convergence threshold timings.
thresh_list = [1e-1, 1e-3, 1e-5];
p_dat = find_thresh(p_dat, thresh_list);

ydim = length(p_dat(1).px)/xdim;
dx = 0.375;
dy = 0.375;

% plot
leg_txt = '';
figure(1)
tiledlayout( ceil( length(mu) / 2 ), length(mu) - floor(length(mu) / 2 ), 'tilespacing', 'compact')
% marker = {'^-','v-','sq-'};
mkrsz = 10;
lw = 1;
fsz = 14;


for j = 1:length(mu)
    nexttile
    for i = 1:length(p_init)
        semilogy( p_dat( i, j ).rms_p , 'markersize', mkrsz, 'linewidth', lw)
            if length(p_init) == 1
                increase_itr = find( (( p_dat(i, j).rms_p(1:end-1) - p_dat(i, j).rms_p(2:end) ) < 0), 1, 'first');
            else
                xlim([0, 100])
                ylim([6e-6, 1])                
            end
        hold on
        leg_text{i} = p_label( p_dat( i, j ).init, p_dat( i, j ).mu );
        xlabel('iterations')
        ylabel('\Deltap / \nablam')
        set(gca, 'fontsize', fsz)

    end
%     if length(p_init == 1)
%         plot(increase_itr, p_dat(i, j).rms_p(increase_itr + 1), 'r*', 'markersize', mkrsz)
%         leg_text = [leg_text, {'end of monotonic regime'}];
%         yticks([10.^-[15:-3:0]])
%     end
    legend(leg_text , 'interpreter', 'latex')    
    hold off
    leg_text = '';
end
hold off


% calulate objective function terms
i = 1;
sum_factor = ones(size(p_dat(1,1).px)); % can be used to isolate only part of the model  
for j = 1:length( mu )
    l1_obj(j) = sum( sum_factor .* sqrt( ( cx * reshape( log10(model)', 1, [] )' - p_dat( i, j ).px ).^2 + ( cy * reshape( log10(model)', 1, [] )' - p_dat( i, j ).py ).^2 ) );
    tgv_obj(j) = sum( sum_factor .* sqrt( ( cx' * p_dat( i, j ).px ).^2 + ( cy' * p_dat( i, j ).py).^2 ) + 0.5 * ( cx' * p_dat( i, j ).py + cy' * p_dat( i, j ).px ).^2 );
end


    figure(10)
    plot( mu, l1_obj, '^-')
    hold on    
    plot( mu, tgv_obj, 'v-')
    plot( mu, l1_obj + mu.*tgv_obj, 'd-')
    legend('|\nabla\bf{m} - \bf{p}|', '0.5|\nabla\bf{p} + \nabla\bf{p}^T|', '|\nabla\bf{m} - \bf{p}| + 0.5\mu|\nabla\bf{p} + \nabla\bf{p}^T|')

% hold on
% plot( p_dat( 1, : ).mu, p_dat( 1, : ).tgv_obj )
% plot ( p_dat( 1, : ).mu, p_dat( 1, : ).l1_obj + p_dat( 1, : ).mu * p_dat( 1, : ).tgv_obj )
% hold off

end


function p_dist = p_analysis(p_dat, model, cx, cy, xdim, mu, p_init)
% Requires all mu and p_init values to be present
% p_dist(1,j) contains the rms distance between the first and second p initialisations for the jth value of mu 
% p_dist(2,j) contains the rms distance between the first and third p initialisations for the jth value of mu 
% p_dist(3,j) contains the rms distance between the second and third p initialisations for the jth value of mu 



model_vec = log10(reshape(model', [], 1));

grad_cx = cx*model_vec;
grad_cy = cy*model_vec;
ydim = length(p_dat(1).px)/xdim;
dx = 0.375;
dy = 0.375;
z_cen = 100*ones(size(model));
pz_im = zeros(size(model));

% calculate distance between initialisations.


for j = 1:length(mu)
    lp = length(mu);
    p_dist(1, j) = sqrt( mean( (p_dat( 2, j ).px - p_dat( 1, j ).px).^2 + (p_dat( 2, j ).py - p_dat( 1, j ).py).^2 ) );
    p_dist(2, j) = sqrt( mean( (p_dat( 3, j ).px - p_dat( 1, j ).px).^2 + (p_dat( 3, j ).py - p_dat( 1, j ).py).^2 ) );
    p_dist(3, j) = sqrt( mean( (p_dat( 3, j ).px - p_dat( 2, j ).px).^2 + (p_dat( 3, j ).py - p_dat( 2, j ).py).^2 ) );
end

% Plot all p plots
figure(2)
tiledlayout(2, 2,'tilespacing', 'compact')

for i = 1%:length(p_init)
    for j = 1:length(mu)
        
        nexttile
        
        px_im = reshape( p_dat( i, j ).px , xdim, [] )';
        py_im = reshape( p_dat( i, j ).py , xdim, [] )';
        
        surf( dx*(1:xdim), dy*(1:ydim), sqrt(px_im.^2 + py_im.^2),'edgecolor','none')
        view([0,0,1])
        title( p_label(p_dat( i, j ).init, p_dat( i, j ).mu), 'interpreter', 'latex')
        colorbar
        set(gca,'ydir','reverse')
        caxis([min(min(sqrt(px_im.^2 + py_im.^2))),max(max(sqrt(px_im.^2 + py_im.^2)))]);
        hold on
        quiver3( dx*(1:xdim) + dx/2, dy*(1:ydim) + dy/2, z_cen, px_im, py_im, pz_im,'color','m');
        axis image
        %             set(gca,'fontsize',fsz)
        hold off
        
        
    end
end

end



function p_dat = find_thresh(p_dat, thresh_list)



figure(3)
hold on

for i = 1:size(p_dat, 1)
    for j = 1:size(p_dat, 2)        
        for k = 1:length(thresh_list)
            p_dat( i, j).thresh(k) = find( [p_dat( i, j ).rms_p, 0] < thresh_list( k ) , 1);
        end
        semilogx( thresh_list, p_dat( i, j).thresh, '.' )
    end            
end

    
end








function label = p_label(p_init, mu)

if p_init == 1
    initialisation = '$\nabla m$';
elseif p_init == 2
    initialisation = 'l2 iteration';
elseif p_init == 3
    initialisation = 'gaussian noise';
elseif p_init ==0
    initialisation = 'zero';
end

label = ['p intialisation: ', initialisation, '; $\mu =$ ', num2str(mu) ];

end


function [pdat, model, cx, cy, xdim] = load_p(folder, p_init, mu)

% preallocations
n_init= length(p_init);
n_mu = length(mu);
k = 1;

for i = 1:length(p_init)
    for j = 1:length(mu)
        
        % load file
        try
            filename = ['p_init', num2str(p_init(i)), '_mu', num2str(mu(j)*10),'.mat'];
            data_in = load(fullfile(folder, filename));
        catch
            warning(['missing file ', filename])
            continue
        end
        if i == 1 && j == 1
            model = data_in.model;
            cx = data_in.cx;
            cy = data_in.cy;
            xdim = data_in.xdim; %used to reshape data once read in
        end
        
        pdat( i, j ).init = p_init(i);
        pdat( i, j ).mu = mu(j);
        pdat( i, j ).px = data_in.px;
        pdat( i, j ).py = data_in.py;
        pdat( i, j ).rms_p = data_in.rms_p;
        
    end
end

end