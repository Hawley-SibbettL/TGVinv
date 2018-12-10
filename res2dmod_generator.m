function res2dmod_generator()
close all
% Variable inputs
username ='bgsvisluke2';
in.filename = ['C:\Users\',username,'\OneDrive - The University of Nottingham\UoN Box Migration\Data\New TGV test models\res2D_input'];
in.n_electrode = 46;
in.electrode_spacing = 0.75;
% % high res
% in.n_electrode = 86;
% in.electrode_spacing = 10*(0.4);

in.depth_values = (0.15:0.15:11);

% Fixed inputs (in most cases)
in.cells_per_electrode = 4;
in.pseudo = 20; % pseudosection data levels
in.offset = 0; % If negative will manually determine the cells to the left and right of the electrodes
in.cells_per_row = (in.n_electrode - 1)*in.cells_per_electrode - (2*in.offset); %shrinks/grows symettrically with the offset
array_type = 3;
in.n_layers = length(in.depth_values);

% IP4DI grid
%load('C:\Users\bgsvisluke2\OneDrive - The University of Nottingham\UoN Box Migration\Coupled TGV results\Simulations\large gamma\Old\vf_gp_l1_05lag_ref','final')
% refined normalmesh
load(['C:\Users\',username,'\OneDrive - The University of Nottingham\UoN Box Migration\Coupled TGV results\Paper\simulations\gauss_fault_l1_15lag_fd'],'final') % ref/not ref to generate refined/normal images
% 2res
% load(['C:\Users\',username,'\OneDrive - The University of Nottingham\UoN Box Migration\Coupled TGV results\Paper\simulations 2res\new_diag_block3_2res_l1_15lag'],'final') % ref/not ref to generate refined/normal images
% Hres
% load(['C:\Users\',username,'\OneDrive - The University of Nottingham\UoN Box Migration\Coupled TGV results\Paper\simulations\block_edge2_Hres2_l22_15lag_fd_gd2'],'final') % ref/not ref to generate refined/normal images



% Create model grid - res_grid gives the resistivities on the lines
res_im = ones(in.n_layers+1,in.cells_per_row+1); % Last row+last column are not used: pcolor maps edge values to cells
bgr_res = 1;


% need grids for the cell edges (for plotting) and the centres (for
% resistivity values)
x_edge = 0:in.electrode_spacing/in.cells_per_electrode:(in.n_electrode-1)*in.electrode_spacing;
y_edge = [0, in.depth_values];
[X_edge, Y_edge] = meshgrid(x_edge,y_edge);
x_centre = [(x_edge(2:end) + x_edge(1:end-1))/2, 0];
y_centre = [(y_edge(2:end) + y_edge(1:end-1))/2, 0];
[X_centre, Y_centre] = meshgrid(x_centre, y_centre);

[Xfine, Yfine] = meshgrid(linspace(x_edge(1),x_edge(end),1000),linspace(y_edge(1),y_edge(end),1000));

% CALL MODEL GENERATING FILE (parameters found in individual functions)
[res_im, thresh_im, in.res_vec, threshold_bounds,ip4di_direct] = vf_g(res_im, X_centre, Y_centre, Xfine,Yfine,final);
% [res_im, thresh_im, in.res_vec, threshold_bounds] = ellipse_gen(res_im, X_centre, Y_centre, Xfine,Yfine);
% [res_im, rand_state] = bedrock_gen(res_im,  X_centre, Y_centre, Xfine, Yfine);
% [res_im, thresh_im, in.res_vec, threshold_bounds] = layers(res_im, X_centre, Y_centre, Xfine,Yfine);
% [res_im, thresh_im, in.res_vec, threshold_bounds] = slipped_fault(res_im, X_centre, Y_centre, Xfine,Yfine);


ip4di_direct = reshape(ip4di_direct,length(unique(final.param_x)),length(unique(final.param_y)))';


figure(1)
surf(X_edge,Y_edge,res_im,'edgecolor','none');view([0,0,1])
colorbar 
axis image
% shading(gca,'interp')
set(gca,'ydir','reverse')

if exist('thresh_im','var')
    figure(3)
    surf(X_edge,Y_edge,thresh_im,'edgecolor','none');view([0,0,1]);
    colorbar
    axis image
    % shading(gca,'interp')
    set(gca,'ydir','reverse')
end

% Interpolate onto Ip4DI grid for comparison (and to save)
[X_ip4di, Y_ip4di] = meshgrid(unique(final.param_x),unique(final.param_y));
ip4di_interpolant = scatteredInterpolant(repmat(x_centre,1,length(y_centre))', reshape(repmat(y_centre,length(x_centre),1), [], 1), reshape(thresh_im', [],1),'linear');
ip4di_im = reshape(ip4di_interpolant(final.param_x, final.param_y),length(unique(final.param_x)),[])';
% To get a thresholded ip4di_im
% if exist('threshold_bounds','var')
%     ip4di_im = imquantize(ip4di_im,threshold_bounds,in.res_vec);
% end

% save('C:\Users\bgsvisluke2\OneDrive - The University of Nottingham\UoN Box Migration\Data\New TGV test models\new_diag_block3_2res','ip4di_im','X-edge','Y_edge','depth','res_im','rho','thresh_im','x_centre','y_centre') %Too risky - use% command window

figure(4)
surf(X_ip4di, Y_ip4di, log10(ip4di_im)); view([0,0,1]);
colorbar 
axis image
% shading(gca,'interp')
set(gca,'ydir','reverse')

if exist('thresh_im','var')
    out_im = thresh_im;
    in.res_vec = unique(thresh_im);
else
    out_im = res_im;
    in.res_vec = unique(res_im);
end

if length(in.res_vec)>16
    error('Too many res levels for resinv')
end

in.res = repmat('z',in.n_layers,in.cells_per_row);

res_key = ['0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'];

for i = 1:in.cells_per_row
    for j = 1:in.n_layers
        for k = 1:length(in.res_vec)
            if out_im(j,i) == in.res_vec(k);
                in.res(j,i) = res_key(k);
            end
        end
    end
end

% Write to file
fid = fopen([in.filename,'.mod'],'w');

write_header(fid,in,1);

fprintf(fid,'%d\n%d\n%d\n%d',array_type,0,0,0); % 3 is for dipole dipole survey, zeros mark end of file
fclose(fid);




% ################# Used to create original template ######################
% Variable inputs
% in.filename = 'res2D_empty';
% in.n_electrode = 35;
% in.res_vec = [50, 100];  % Vector of resistivities in Ohm m
% in.electrode_spacing = 1;
% in.n_layers = 14;
% in.depth_values = [];
% in.pseudo = 10; % pseudosection data levels
% 
% % Fixed inputs (in most cases)
% in.cells_per_electrode = 4;
% in.cells_per_row = (in.n_electrode - 1)*in.cells_per_electrode;
% 
% % Generate empty grid for resistivities
% in.res = repmat('0',in.n_layers,in.cells_per_row);
% in.res(1) = '1';
% 
% 
% 
% % Write to file
% fid = fopen([in.filename,'.mod'],'w');
% 
% write_header(fid,in,0);
% 
% fprintf(fid,'%d\n%d\n%d\n%d',6,0,0,0); % six is for dipole dipole survey, zeros mark end of file
% fclose(fid);
depth = y_centre;
rho = thresh_im(:,1);



end



function write_header(fid,in,depth_flag)  % depth_flag=1 if manually specify depths

% title to electrode spacing
fprintf(fid, ['%s\n', '%d\n', '%d\n', '%d\n', '%1.2f\n'],in.filename,in.n_electrode,in.pseudo,0,in.electrode_spacing);

% Manual depth flag
if depth_flag == 1
    fprintf(fid,'%d\n',2);
else
    fprintf(fid,'%d\n',0);
end

% offset, n_columns, n_res, cells per electrode, resistance values, number
% of rows
fprintf(fid,['%d\n%d\n%d\n%d\n','%.4f',repmat(',%.4f',1,length(in.res_vec)-1),'\n%d\n'],in.offset,in.cells_per_row,length(in.res_vec),in.cells_per_electrode,in.res_vec,in.n_layers);

% depth values
if depth_flag == 1
    fprintf(fid,['%.4f',repmat(',%.4f',1,length(in.depth_values)-1),'\n'],in.depth_values);
end

% main resistivity block
fprintf(fid,[repmat('%c',1,in.cells_per_row),'\n'],in.res');

end

function [res_im, thresh_im, res_vec, threshold_bounds,ip4di] = vf_g(res_im,  X_centre, Y_centre, Xfine, Yfine, final)
% Vertical Gaussian fault (can have linear transition) plus gaussians

num_thresh = 14;

% Gaussian parameters
mu1 = [30,.6];%[30, 1];   % position of gaussian centre (m)
sigma1 = [5, 1];%[5, 1.5];      % std, in m. Note - will not be a true gaussian due to block lengthening
A1 = 100;         % Amplitude of gaussian inclusion above background (ohm m)
mu2 = [20, 1.8];
sigma2 = [6, 1];
A2 = 100;
mu3 = [12, 3];
sigma3 = [5, 1];     
A3 = 100;        
mu4 = [4, 1.8];
sigma4 = [3.5, 0.7];
A4 = 100;
mu5 = [16, 2.2];
sigma5 = [4, 0.8];
A5 = 0;

block_flag = 1;

% linear layer
% horizontal layer with linear transition from rho1-rho2
rho1 = 10;  % Upper res
rho2 = 10;  % Lower res
ub = 10;     % layer ub
lb = 30;   % layer lb
res_im(X_centre <= ub) = rho1;
res_im(X_centre >= lb) = rho2;
transition_zone = (X_centre <= lb)&(X_centre > ub);
% res_im(transition_zone) = ((rho1) + ((rho2)-(rho1)).*(X_centre(transition_zone)-ub)/(lb-ub));
res_im(transition_zone) = 10.^(log10(rho1) + (log10(rho2)-log10(rho1)).*(X_centre(transition_zone)-ub)/(lb-ub));

gaussian_inclusion = Gauss(mu1,sigma1,A1,Xfine,Yfine) + Gauss(mu2,sigma2,A2,Xfine,Yfine) + Gauss(mu3,sigma3,A3,Xfine,Yfine) + Gauss(mu4,sigma4,A4,Xfine,Yfine) + Gauss(mu5,sigma5,A5,Xfine,Yfine);
res_im = res_im + interp2(Xfine,Yfine,gaussian_inclusion,X_centre,Y_centre,'linear');
% res_im(Y_centre <= lb) = rho1; % for hard edge


bgr_res = [rho1, rho2];

% Discretise resistivity levels for res2Dmod (max 16)
threshold_bounds = multithresh((res_im),num_thresh);
% threshold_bounds = 10.^multithresh(log10(res_im),num_thresh);
res_vec = [threshold_bounds(1),(threshold_bounds(2:end) + threshold_bounds(1:end-1))/2,threshold_bounds(end)];

% % make sure background resitivities are reflected in the theshold
% for i = 1:length(bgr_res)
%     [~, bgr_ind] = min(abs(res_vec-bgr_res(i)));
%     res_vec(bgr_ind) = bgr_res(i);    
% end

thresh_im = imquantize(res_im,threshold_bounds,res_vec);

% Add solid inclusion afterwards
% res_im(interp2(Xfine,Yfine,Gauss(mu5,sigma5,A5,Xfine,Yfine),X_centre,Y_centre,'linear')>0.33*A5) = rho1 + A5;
% thresh_im(interp2(Xfine,Yfine,Gauss(mu5,sigma5,A5,Xfine,Yfine),X_centre,Y_centre,'linear')>0.33*A5) = rho1 + A5;

% Add a solid block
if block_flag == 1
    rho_b = 300;
    xbC = 21; % x block centre
    xbL = 3; % x block length
    ybC = 2;
    ybL = 2.2;
    xb = [xbC - xbL/2, xbC + xbL/2];
    yb = [ybC - ybL/2, ybC + ybL/2];
    res_im(X_centre > xb(1) & X_centre < xb(2) & Y_centre > yb(1) & Y_centre < yb(2)) = rho_b;
%     thresh_im(:,:) = rho1; % for block on solid background
    thresh_im(X_centre > xb(1) & X_centre < xb(2) & Y_centre > yb(1) & Y_centre < yb(2)) = rho_b;
end
% add a second, inner block (for gradual boundary)
% rho_b = 100;
% xbL2 = xbL - 1.1; % x block length
% ybL2 = ybL - 0.91;
% xb = 10*[xbC - xbL2/2, xbC + xbL2/2];
% yb = 10*[ybC - ybL2/2, ybC + ybL2/2];
% res_im(X_centre > xb(1) & X_centre < xb(2) & Y_centre > yb(1) & Y_centre < yb(2)) = rho_b;
% thresh_im(X_centre > xb(1) & X_centre < xb(2) & Y_centre > yb(1) & Y_centre < yb(2)) = rho_b;


% Add solid layers afterwards
% rho_b = 140;
% b = 0.8;
% res_im(Y_centre < ub) = rho_b;
% thresh_im(Y_centre < ub) = rho_b;

% res_vec = [res_vec, rho1+A5];
% threshold_bounds = [threshold_bounds, rho1 + A5 - 10];

% directly calculate onto ip4di mesh
ip4di = final.res_param1(:,1);
ip4di(final.param_x <= ub) = rho1;
ip4di(final.param_x >= lb) = rho2;
transition_zone = (final.param_x <= lb)&(final.param_x > ub);
ip4di(transition_zone) = 10.^(log10(rho1) + (log10(rho2)-log10(rho1)).*(X_centre(transition_zone)-ub)/(lb-ub));

ip4di = ip4di + Gauss(mu1,sigma1,A1,final.param_x,final.param_y) + Gauss(mu2,sigma2,A2,final.param_x,final.param_y) + Gauss(mu3,sigma3,A3,final.param_x,final.param_y) + Gauss(mu4,sigma4,A4,final.param_x,final.param_y) + Gauss(mu5,sigma5,A5,final.param_x,final.param_y);

if block_flag == 1
    ip4di(final.param_x > xb(1) & final.param_x < xb(2) & final.param_y > yb(1) & final.param_y < yb(2)) = rho_b;
end

threshold_bounds_ip4di = multithresh((ip4di),num_thresh);
% threshold_bounds = 10.^multithresh(log10(res_im),num_thresh);
res_vec_ip4di = [threshold_bounds(1),(threshold_bounds(2:end) + threshold_bounds(1:end-1))/2,threshold_bounds(end)];
ip4di = imquantize(ip4di,threshold_bounds_ip4di,res_vec_ip4di);

end

function [res_im, thresh_im, res_vec, threshold_bounds] = ellipse_gen(res_im,  X_centre, Y_centre, Xfine, Yfine)

num_thresh = 16;

rho1 = 30; % ellipse res
rho2 = 120; % background res
rho_step = 70; % res level at border of ellipse, linearly drops to rho2
mu = [16, 0]; % ellipse centre [x, y]
a = 8; % width
b = 1.7; % height
a_edge = 2; % width of boundary
b_edge = 2;

res_im = res_im*rho2;
ellipse =  ((X_centre - mu(1))/a).^2 + ((Y_centre - mu(2))/b).^2 < 1;
trans_zone = ( ((X_centre - mu(1))/a).^2 + ((Y_centre - mu(2))/b).^2 > 1 ) & ( ((X_centre - mu(1))/(a + a_edge)).^2 + ((Y_centre - mu(2))/(b + b_edge)).^2 < 1 );
trans_zone_ind = find(trans_zone);
trans_mask = zeros(size(res_im));

% find distance from transition zone edge to cells
for i = 1:length(trans_zone_ind)
    % find nearest ellipse point
    xd = (X_centre - X_centre(trans_zone_ind(i)));
    yd =  (Y_centre - Y_centre(trans_zone_ind(i)));
    dist = (xd.^2 + yd.^2);
    dist(~ellipse) = 1000;
    trans_mask(trans_zone_ind(i)) = min(min(dist));
end
trans_mask = trans_mask./max(max(trans_mask));
% res_im(trans_zone) = 10.^(log10(rho1) + trans_mask(trans_zone).*(log10(rho2) - log10(rho1)));
res_im(trans_zone) = 10.^(log10(rho_step) + trans_mask(trans_zone).*(log10(rho2) - log10(rho_step)));
res_im(ellipse) = rho1;

% Add gaussians if desired
mu1 = [19, 2];   
sigma1 = [1.5, 0.7];      
A1 = 0; 
mu2 = [12, 0.5];   
sigma2 = [2, 1];      
A2 = 0; 
gaussian_inclusion = Gauss(mu1,sigma1,A1,Xfine,Yfine) + Gauss(mu2,sigma2,A2,Xfine,Yfine);
res_im = res_im + interp2(Xfine,Yfine,gaussian_inclusion,X_centre,Y_centre,'linear');



bgr_res = [rho1, rho2];

% Discretise resistivity levels for res2Dmod (max 16)
threshold_bounds = 10.^multithresh(log10(res_im),num_thresh-1);
res_vec = [threshold_bounds(1),(threshold_bounds(2:end) + threshold_bounds(1:end-1))/2,threshold_bounds(end)];



% make sure background resitivities are reflected in the theshold
for i = 1:length(bgr_res)
    [~, bgr_ind] = min(abs(res_vec-bgr_res(i)));
    res_vec(bgr_ind) = bgr_res(i);    
end


% in.res_vec = bgr_res;
thresh_im = imquantize(res_im,threshold_bounds,res_vec);

% Add a solid ellipse, resistivity rho3
% rho3 = 0; % background res
% mu = [11, 1.3]; % ellipse centre [x, y]
% a = 2; % width
% b = 1; % height
% 
% % Add after imquantise to allow high contrast
% res_im(((X_centre - mu(1))/a).^2 + ((Y_centre - mu(2))/b).^2 < 1) = rho3;
% thresh_im(((X_centre - mu(1))/a).^2 + ((Y_centre - mu(2))/b).^2 < 1) = rho3;
% res_vec = [res_vec, rho3];
% threshold_bounds = [threshold_bounds, rho3 - 1];


end



function [res_im, rand_state] = bedrock_gen(res_im,  X_centre, Y_centre, Xfine, Yfine)

% jagged bedrock
rho1 = 60;
rho2 = 50;
init_depth = 2;
step = 0.4;
seclength = 2;
% [~, depth_ind] = min(abs(in.depth_values - init_depth));
bedrock_depth = zeros(size(res_im,2),1);
bedrock_depth(1:seclength) = init_depth;
rand_state = rng;
load('C:\Users\bgsvisluke2\OneDrive - The University of Nottingham\UoN Box Migration\Data\New TGV test models\bedrock_aw2_ref','rand_state'); rng(rand_state);
rand_vec = rand(size(res_im,2),1);


for i = seclength+1:seclength:length(bedrock_depth)
    if rand_vec(i) < 0.33
        bedrock_depth(i:i+seclength-1) = bedrock_depth(i-seclength);
    elseif rand_vec(i) < 0.66
        bedrock_depth(i:i+seclength-1) = bedrock_depth(i-seclength) - step;
    else
        bedrock_depth(i:i+seclength-1) = bedrock_depth(i-seclength) + step;
    end
        
    if bedrock_depth(i) < 0.5
        bedrock_depth(i:i+seclength) = 0.5;
    end
end
if length(bedrock_depth)>size(res_im,2)
    bedrock_depth(end - (size(res_im,2) )) = [];
end

res_im(Y_centre < repmat(bedrock_depth,1,size(res_im,1))') = rho1;
res_im(Y_centre >= repmat(bedrock_depth,1,size(res_im,1))') = rho2;

% Add weathering
weather_depth = 0.8;
uniform_weathering  = (Y_centre >= repmat(bedrock_depth,1,size(res_im,1))') & (Y_centre < repmat(bedrock_depth + weather_depth,1,size(res_im,1))');
bw_thresh = mean(bedrock_depth);
broken_weathering = (Y_centre >= repmat(bedrock_depth,1,size(res_im,1))') & (Y_centre < bw_thresh);

% Choose type of weathering here
res_im(broken_weathering) = mean([rho1,rho2]);

end

function [res_im, thresh_im, res_vec, threshold_bounds] = layers(res_im,  X_centre, Y_centre, Xfine, Yfine)
% 
% % Solid layers
% % set up layer depths
% depth = [0.8, 1.9];
% rho = [50, 300, 100];
% num_thresh = length(rho);
% 
% res_im = res_im*rho(1);
% 
% for i = 1:length(depth)
%     res_im(Y_centre > depth(i)) = rho(i+1);
% end
% for i = 1:length(depth)
%     res_im(Y_centre > depth(i)) = rho(i+1);
% end
% 
% 
% % Discretise resistivity levels for res2Dmod (max 16)
% threshold_bounds = [rho(1)-0.1, rho+0.1];
% res_vec = rho;
% thresh_im = res_im;


% layersG

num_thresh = 15;

% depth = [0.9, 2.2];
depth = [0.2, 3.7];
% rho = log10([50, 300, 150]);
rho = log10([250, 250, 80]); % two layer
stdG = 1.6;

depth_y = Y_centre(1:end-1,1);
rho_y = ones(size(depth_y));
rho_y(depth_y<=depth(1)) = rho(1);
rho_y(depth_y>depth(1) & depth_y<=depth(2)) = rho(2);
rho_y(depth_y>depth(2)) = rho(3);
% GL = exp(-(depth_y - depth(2)).^2./(2*(stdG).^2));
% GL = (rho(2)-rho(3))*GL/max(GL);
% rho_y(depth_y>depth(2)) =  rho(3) + GL(depth_y>depth(2));
rho_y = 10.^rho_y;

% Fix to 16 levels
rho_y = round(rho_y,0);
resLevels = unique(rho_y);
lvlDiff = length(resLevels)-16;

if lvlDiff>0
    [lvl, lvlInd] = sort(abs(resLevels-10.^rho(3)));
    for q = 2:lvlDiff+1
        rho_y(rho_y == resLevels(lvlInd(q))) = round(10.^rho(3),0);
    end
end


figure(95)
plot(depth_y,rho_y)
figure(96)
plot(depth_y,log10(rho_y),'r')


res_im = res_im.*[rho_y;rho_y(end)];

thresh_im = res_im;

res_vec = unique(rho_y);
threshold_bounds = [];


% % Discretise resistivity levels for res2Dmod (max 16)
% threshold_bounds = 10.^multithresh(log10(res_im),num_thresh);
% res_vec = [threshold_bounds(1),(threshold_bounds(2:end) + threshold_bounds(1:end-1))/2,threshold_bounds(end)];
% 

% make sure background resitivities are reflected in the theshold
% for i = 1:length(rho)
%     [~, bgr_ind] = min(abs(res_vec-rho(i)));
%     in.res_vec(bgr_ind) = rho(i);    
% end

% thresh_im = imquantize(res_im,threshold_bounds,res_vec);


end

function [res_im, thresh_im, res_vec, threshold_bounds] = slipped_fault(res_im,  X_centre, Y_centre, Xfine, Yfine)
    % Generated 2 layers for a slipped fault, each with a linear transition
    % to the background. the layers are on the rhs of the image


    num_thresh = 15;
    
    % Background resistivity
    rho_bg = 200;
    res_im(:,:) = rho_bg;
    layer_bound = 2;
    
    
    % Top layer
    rho_up = 200;
    ub_left = 21;
    ub_right = 21;
    res_im(X_centre > ub_right & Y_centre < layer_bound) = rho_up;
    upper_trans = (X_centre < ub_right & X_centre > ub_left & Y_centre < layer_bound);
    res_im(upper_trans) = 10.^(log10(rho_up) + (log10(rho_bg)-log10(rho_up)).*(X_centre(upper_trans)-ub_right)/(ub_left-ub_right));;
    
    % Lower layer
    rho_dn = 700;
    rho_dn_int = 10.^log10((rho_dn + rho_bg)/2);
    lb_left = 9;
    lb_int = 9;
    lb_right = 20;
    res_im(X_centre > lb_right & Y_centre > layer_bound) = rho_dn;
    lower_transL = (X_centre < lb_int & X_centre > lb_left & Y_centre > layer_bound);
    res_im(lower_transL) = 10.^(log10(rho_dn_int) + (log10(rho_bg)-log10(rho_dn_int)).*(X_centre(lower_transL)-lb_int)/(lb_left-lb_int));    
    lower_transR = (X_centre < lb_right & X_centre > lb_int & Y_centre > layer_bound);
    res_im(lower_transR) = 10.^(log10(rho_dn) + (log10(rho_dn_int)-log10(rho_dn)).*(X_centre(lower_transR)-lb_right)/(lb_int-lb_right));    

    
    
    % Discretise resistivity levels for res2Dmod (max 16)
%     threshold_bounds = 10.^multithresh(log10(res_im(1:2:end,1:2:end)),num_thresh);
    threshold_bounds = sort([rho_bg, linspace((rho_dn_int),(rho_dn),13)]);
    res_vec = [threshold_bounds(1),(threshold_bounds(2:end) + threshold_bounds(1:end-1))/2,threshold_bounds(end)];
    thresh_im = imquantize(res_im,threshold_bounds,res_vec);


end




function gaussian = Gauss(mu,sigma,A,Xfine,Yfine)

    gaussian = A*exp(-(((Xfine) - mu(1)).^2./(2*sigma(1).^2) + ((Yfine) - mu(2)).^2/(2*sigma(2).^2)));


end
