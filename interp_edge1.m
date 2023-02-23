
% interpolates function u onto the edges of the grid defined by
% x,y = param_x,param_y. Homogeneous neumann bc are used to find the edges
% on the domain boundary
% each u is [u, ux, uy]
% uh; uv - cell centres, with ghost points added to the [upper/lower; left/right] boundaries respectively
% uch; ucv - cell corners, including points on the [upper/lower; left/right] boundaries respectively
% ueh; uev - horizontal and vertical cell edges respectively. Includes boundary edges
function [uh, uv, uch, ucv, ueh, uev] = interp_edge1(u, x, y)

% The unequal cell spacings make this a little bit trickier than for a
% uniform mesh, as the interpolation requires the cell spacings at each
% point.

% key values
tmp_x = unique(x);
tmp_y = unique(y);
xL = length(tmp_x);
yL = length(tmp_y);

% uh, uv - just add ghost points in each dimension
[x_ghostx, x_ghosty, ux_ghost] = param_generator([2*tmp_x(1) - tmp_x(2); tmp_x; 2*tmp_x(end) - tmp_x(end-1)], tmp_y, u, [1, 0]);
[y_ghostx, y_ghosty, uy_ghost] = param_generator(tmp_x, [2*tmp_y(1) - tmp_y(2); tmp_y; 2*tmp_y(end) - tmp_y(end-1)], u, [0, 1]);

uh = [uy_ghost, y_ghostx, y_ghosty];
uv = [ux_ghost, x_ghostx, x_ghosty];

% calculate uv, uh, corner grids

% add ghost points at either end of the original mesh 
[x_ghost1, ~, ux_ghost1] = param_generator([2*tmp_x(1) - tmp_x(2); tmp_x], tmp_y, u, [1, 0]);
[x_ghost2, ~, ux_ghost2] = param_generator([tmp_x; 2*tmp_x(end) - tmp_x(end-1)], tmp_y, u, [2, 0]);
[~, y_ghost1, uy_ghost1] = param_generator(tmp_x, [2*tmp_y(1) - tmp_y(2); tmp_y], u, [0, 1]);
[~, y_ghost2, uy_ghost2] = param_generator(tmp_x, [tmp_y; 2*tmp_y(end) - tmp_y(end-1)], u, [0, 2]);

% interpolate x,y onto edges
x_edge = (x_ghost1 + x_ghost2)/2;
y_edge = (y_ghost1 + y_ghost2)/2;

% parameter meshes for new grids
% [x_edge, ~] = param_generator(tmp_x_edge, tmp_y);
% [~, y_edge] = param_generator(tmp_x, tmp_y_edge);
% [xc, yc] = param_generator(tmp_x_edge, tmp_y_edge);

ux_edge = ux_ghost1.*((x_edge - x_ghost1)./(x_ghost2 - x_ghost1)) + ux_ghost2.*((x_ghost2 - x_edge)./(x_ghost2 - x_ghost1));
uy_edge = uy_ghost1.*((y_edge - y_ghost1)./(y_ghost2 - y_ghost1)) + uy_ghost2.*((y_ghost2 - y_edge)./(y_ghost2 - y_ghost1));

[uev_x, uev_y] = param_generator(unique(x_edge), tmp_y);
[ueh_x, ueh_y] = param_generator(tmp_x, unique(y_edge));

ueh = [uy_edge, ueh_x, ueh_y];
uev = [ux_edge, uev_x, uev_y];

% now do the same thing for the corner points

% expand x_edge, y_edge to include the boundary points
[xc_ghost1, ~, uxc_ghost1] = param_generator([2*tmp_x(1) - tmp_x(2); tmp_x], unique(y_edge), uy_edge, [1, 0]);
[xc_ghost2, ~, uxc_ghost2] = param_generator([tmp_x; 2*tmp_x(end) - tmp_x(end-1)], unique(y_edge), uy_edge, [2, 0]);
[~, yc_ghost1, uyc_ghost1] = param_generator(unique(x_edge), [2*tmp_y(1) - tmp_y(2); tmp_y], ux_edge, [0, 1]);
[~, yc_ghost2, uyc_ghost2] = param_generator(unique(x_edge), [tmp_y; 2*tmp_y(end) - tmp_y(end-1)], ux_edge, [0, 2]);

% interpolate x,y onto edges
xc = (xc_ghost1 + xc_ghost2)/2;
yc = (yc_ghost1 + yc_ghost2)/2;

ucx = uxc_ghost1.*((xc - xc_ghost1)./(xc_ghost2 - xc_ghost1)) + uxc_ghost2.*((xc_ghost2 - xc)./(xc_ghost2 - xc_ghost1));
ucy = uyc_ghost1.*((yc - yc_ghost1)./(yc_ghost2 - yc_ghost1)) + uyc_ghost2.*((yc_ghost2 - yc)./(yc_ghost2 - yc_ghost1));

uc = (ucx + ucy)./2;

% remove top/bottom boundary points for ucv
ucv = [uc(xL+2:end-xL-1), xc(xL+2:end-xL-1), yc(xL+2:end-xL-1)];
% remove left/ right boundary points for uch
uc_x = uc;
uc_x(1:xL+1:end) = []; % remove x(-1/2) edges
uc_x(xL:xL:end) = []; % remove x(end+1/2) edges, takinf reduced size into account
% repeat for x, y coordinates
xc_x = xc; xc_x(1:xL+1:end) = []; xc_x(xL:xL:end) = [];
yc_x = yc; yc_x(1:xL+1:end) = []; yc_x(xL:xL:end) = [];

uch = [uc_x, xc_x, yc_x];



end

% tmp_x = unique(x);
% tmp_y = unique(y);
% 
% xL = length(tmp_x);
% yL = length(tmp_y);

% % 
% % y_edge(1:xL) = 1.5*tmp_y(1) - 0.5*tmp_y(2); % 
% % y_edge = [y_edge'; (1.5*tmp_y(end) - 0.5*tmp_y)*ones(length(xL), 1)];
% 
% 
% % define ghost points extending mesh, then calculate interpolated points
% x_ghost1 = ones((xL+1)*yL, 1); % preallocate zeros array
% x_start_bc = 1:xL+1:(xL+1)*yL;
% x_end_bc = (xL+1:xL+1:(xL+1)*yL);
% x_ghost1(~x_ghost1(x_start_bc)) = x; % x0:xn
% x_ghost1(x_start_bc) = 2*tmp_x(1) - tmp_x(2);
% 
% x_ghost2 = zeros((xL+1)*yL, 1); % preallocate zeros array
% x_ghost2(~x_end_bc) = x;% x1:xn+1;
% x_ghost2(x_end_bc) = 2*tmp_x(1) - tmp_x(2);
% 
% y_ghost1 = [repmat(2*tmp_y(1) - tmp_y(2), xL, 1); y]; % y0:yn
% y_ghost2 = [y; repmat(2*tmp_y(end) - tmp_y(end-1), xL, 1)]; % y1:yn+1
%     
% x_edge = (x_ghost1 + x_ghost2)/2; % interpolated edge points
% y_edge = (y_ghost1 + y_ghost2)/2;
% x_diff = x_ghost2 - x_ghost1; % needed for interpolation calc
% y_diff = y_ghost2 - y_ghost1;
% 
% % now interpolate
% % u{i+1\2) = u{i}*(x{i+1/2}-x{i}/(x{i+1}-x{i})) + u{i+1}*{(x{i+1} - x{i+1/2})/(x{i+1} - x{i})}
% 
% % add u ghost points for calculation
% ux_ghost1 = zeros(size(x_edge));
% ux_ghost1(~x_start_bc) = u;
% ux_ghost1(x_start_bc) = u(1:xL:end); % from neuman bc
% ux_ghost2 = zeros(size(x_edge));
% ux_ghost2(~x_end_bc) = u;
% ux_ghost2(x_end_bc) = u(xL:xL:end); % from neuman bc
% 
% uy_ghost1 = [u(1:xL); u];
% uy_ghost2 = [u; u(end-xL+1:end)];
% 
% % interpolate u onto edges
% ux_edge = ux_ghost1.*((x_edge - x_ghost1)./x_diff) + ux_ghost2.*((x_ghost2 - x_edge)./x_diff);
% uy_edge = uy_ghost1.*((y_edge - y_ghost1)./x_diff) + uy_ghost2.*((y_ghost2 - y_edge)./y_diff);
% 
% 
% 
% % now calculate corner points in the same way. First, ghost points need
% % adding to a grid of the edge coordinates
% 
% 
% % add an extra y row to x arrays to match the dimensions of the corner grid
% xc_ghost1 = [x_ghost1; xc_ghost1(1:xL+1)];
% xc_ghost2 = [x_ghost2; xc_ghost2(1:xL+1)];
% 
% xc_diff = xc_ghost2 - xc_ghost1;
% x_corner = [x_edge, x_edge(1:xL+1)];
% 
% % add an extra x point to each row of the y arrays to the match corner grid
% % dimensions
% yc_start_bc = (1:xL+1:(xL+1)*(yL+1));
% yc_end_bc = (xL+1:xL+1:(xL+1)*(yL+1));
% 
% yc_ghost1 = zeros((xL+1)*(yL+1), 1); 
% yc_ghost1(~yc_start_bc) = y_ghost1; % x0:xn
% yc_ghost1(xy_start_bc) = y_ghost1(1:xL:end);
% 
% yc_ghost2 = zeros((xL+1)*(yL+1), 1); 
% yc_ghost2(~yc_end_bc) = y_ghost2; % x0:xn
% yc_ghost2(xy_end_bc) = y_ghost2((xL:xL:end));
% 
% yc_diff = yc_ghost2 - yc_ghost1; 
% 
% y_corner = zeros((xL+1)*(yL+1), 1); 
% y_corner(~yc_start_bc) = y_edge; % x0:xn
% y_corner(xy_start_bc) = y_edge(1:xL:end);
% 
% % Use the bc to get the values of the ghost points on the edges
% uy_corner1 = zeros((xL+1)*(yL+1) ,1);
% uy_corner1(~uy_start_bc) = uy_edge;
% uy_corner1(uy_start_bc) = uy_edge(1:xL:end);
% 
% uy_corner2 = zeros((xL+1)*(yL+1) ,1);
% uy_corner2(~uy_end_bc) = uy_edge;
% uy_corner2(uy_end_bc) = uy_edge(xL:xL:end);
% 
% ux_corner1 = [ux_edge(1:xL+1); ux_edge];
% ux_corner2 = [ux_edge; ux_edge(end-xL+1:end)];
% 
% 
% ucx = uy_corner1*((x_corner - xc_ghost1)./xc_diff) + uy_corner2*((xc_ghost2 - x_corner)./xc_diff);
% ucy = ux_corner1*((y_corner - yc_ghost1)./yc_diff) + ux_corner2*((yc_ghost2 - y_corner)./yc_diff);
% 
% uc = (ucx + ucy)./2;
% 


