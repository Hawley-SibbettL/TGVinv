% Used to define a forward model on a grid (for no topography)
% Note: image is displayed on a regularised grid, but wih correct
% co-ordinates (which represent he centre of the blocks)
% Read mesh data from grid
% clear all
load('fwd_mesh.mat')

% Plot (imagsc, square grid- ideally with custom ticks) 

% Calculate x and y dimensions of mesh (in pixels)
xfind = find(mesh.param_x == mesh.param_x(1),2);
len_xdim = xfind(2) - xfind(1);
len_ydim = mesh.num_param./len_xdim;

% mesh.init_res_param = 1:mesh.num_param; % Use to check plotting order ok
res_image = reshape(mesh.init_res_param,len_xdim,len_ydim)';
[xgrid, ygrid] = meshgrid(unique(mesh.param_x), unique(mesh.param_y));

figure(1)
imagesc(res_image); colorbar; colormap cool;
ax = gca;
set(ax,'Xtick',(1:len_xdim) ,'Ytick',(1:len_ydim) ,...
    'XtickLabel',num2str(mesh.param_x(1:len_xdim)),'Yticklabel',num2str(mesh.param_y(1:len_xdim:end)));
grid on

%%
% Insert model code here! (maybe use flags to pick from standard options in
% the future)

% Layered model - Finningley
% res_image(1:2,:) = 62;
% res_image(3:4,:) = 15;
% res_image(5:7,:) = 125;
% res_image(8:end,:) = 50;

% Finningley varients
% res_image(2,:) = 30.5;  % Wuu1
% res_image(3,:) = 30.5;  % Wud1
% res_image(8,:) = 80;  % Wdd1
% res_image(7,:) = 80;    % Wdu1
% res_image(4,:) = 10.^1.5; res_image(5,:) = 10.^1.8; % Wmu2



% Load from dd records 
% load_background
% res_image = VFHS_sqtrans_gf2;
%%

% Plot model
figure(2)
imagesc(log10(res_image)); colorbar; colormap cool;
ax = gca;
set(ax,'Xtick', (1:len_xdim), 'Ytick', (1:len_ydim),...
    'XtickLabel',num2str(mesh.param_x(1:len_xdim)),'Yticklabel',num2str(mesh.param_y(1:len_xdim:end)));
grid on

% Pixel Plot
figure(3)
surf(unique(mesh.param_x),unique(mesh.param_y),log10(res_image),'edgecolor','none')
set(gca,'ydir','reverse')
view([0,0,1])
colormap cool
colorbar 
axis image


% Save model
mesh.init_res_param = reshape(res_image',[],1);
mesh.res_param1 = mesh.init_res_param;
mesh.res_param2 = mesh.res_param1;

save('fwd_model')


% figure(3)
% imagesc(reshape(mesh.res_param1,len_xdim,len_ydim)'); grid on; colorbar; colormap cool;
% ax = gca;
% set(ax,'Xtick',(1:len_xdim) ,'Ytick',(1:len_ydim) ,...
%     'XtickLabel',num2str(mesh.param_x(1:len_xdim)),'Yticklabel',num2str(mesh.param_y(1:len_ydim:end)));
