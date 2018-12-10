% Checks the cells used in the roughness matrix
close all

%create smoothness matrix
mesh=smooth_mtx_surface4(input,mesh);

% Plot resistance image

% Calculate x and y dimensions of mesh (in pixels)
xfind = find(mesh.param_x == mesh.param_x(1),2);
len_xdim = xfind(2) - xfind(1);
len_ydim = mesh.num_param./len_xdim;

res_image = reshape(1:mesh.num_param,len_xdim,len_ydim);

figure(1)
imagesc(mesh.param_x(1:len_xdim),mesh.param_y(1:len_xdim:end),res_image);
colorbar; colormap cool;
% ax = gca;
% set(ax,'Xtick',mesh.param_x(1:len_xdim),'Ytick',mesh.param_y(1:len_xdim:end),...
%  'XtickLabel',num2str(mesh.param_x(1:len_xdim)),'Yticklabel',num2str(mesh.param_y(1:len_xdim:end)));
% %'Xtick',(1:len_xdim) ,'Ytick',(1:len_ydim) ,...
%     
% grid on
hold on

% Cycle through rows of cx, cy showing the points at which the matrix has
% value
current_minus = find(mesh.cx(1,:) < 0);
current_plus = find(mesh.cx(1,:) > 0);

cdata = plot([mesh.param_x(current_minus), mesh.param_x(current_plus)],[mesh.param_y(current_minus),  mesh.param_y(current_plus)],'kv');

% for i = 2:mesh.num_param
%     pause(0.01)
%     current_minus = find(mesh.cx(i,:) < 0);
%     current_plus = find(mesh.cx(i,:) > 0);
%     
%     if isempty(current_minus) && isempty(current_plus)
%         disp(['Empty row! i = ', num2str(i)]) 
%     end
%     
%     set(cdata,'Xdata',[mesh.param_x(current_minus), mesh.param_x(current_plus)],...
%         'Ydata',[mesh.param_y(current_minus),  mesh.param_y(current_plus)])
% end

current_minus = find(mesh.cy(1,:) < 0);
current_plus = find(mesh.cy(1,:) > 0);

cdata = plot([mesh.param_x(current_minus), mesh.param_x(current_plus)],[mesh.param_y(current_minus),  mesh.param_y(current_plus)],'kv');

for i = 1:mesh.num_param
    pause(0.01)
    current_minus = find(mesh.cy(i,:) < 0);
    current_plus = find(mesh.cy(i,:) > 0);
    
    if isempty(current_minus) && isempty(current_plus)
        disp(['Empty row! i = ', num2str(i)]) 
    end
    
    set(cdata,'Xdata',[mesh.param_x(current_minus), mesh.param_x(current_plus)],...
        'Ydata',[mesh.param_y(current_minus),  mesh.param_y(current_plus)])
end

