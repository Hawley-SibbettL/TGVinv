% Performs half point-difference calculation
% p2 = [px; py]
% m - log10(resistivty), mu - tgv_lagrn
function [p2] = p_calc(mesh, lagrn ,tgv_lagrn, gamma_c, gamma_p)

px = mesh.px;
py = mesh.py;

% #####################################################################
% Here calculate variables as in ctc_calc
% key values
m = log10(mesh.res_param1);
tmp_x = unique(mesh.param_x);
tmp_y = unique(mesh.param_y);
xL = length(tmp_x);
yL = length(tmp_y);
dx = (tmp_x(2) - tmp_x(1)); % cell width 
tmp_xe = mesh.xe;
tmp_ye = mesh.ye;


% Create shifted versions of m, x, y, px, py
[~, ~, mx_ghost1] = param_generator([tmp_x(1) - dx; tmp_x(1:end-1)], tmp_y, m, [3, 0]);
[~, ~, mx_ghost2] = param_generator([tmp_x(2:end); tmp_x(end) + dx], tmp_y, m, [4, 0]);
[~, ~, my_ghost1] = param_generator(tmp_x, [2*tmp_y(1) - tmp_y(2); tmp_y(1:end-1)], m, [0, 3]);
[~, ~, my_ghost2] = param_generator(tmp_x, [tmp_y(2:end); 2*tmp_y(end) - tmp_y(end-1)], m, [0, 4]);

[x_ghost1, ~, px_ghost1] = param_generator([tmp_x(1) - dx; tmp_x(1:end-1)], tmp_y, px, [3, 0]);
[x_ghost2, ~, px_ghost2] = param_generator([tmp_x(2:end); tmp_x(end) + dx], tmp_y, px, [4, 0]);
[~, y_ghost1, px_ghost_t] = param_generator(tmp_x, [2*tmp_y(1) - tmp_y(2); tmp_y(1:end-1)], px, [0, 3]);
[~, y_ghost2, px_ghost_b] = param_generator(tmp_x, [tmp_y(2:end); 2*tmp_y(end) - tmp_y(end-1)], px, [0, 4]);

[~, ~, py_ghost1] = param_generator([tmp_x(1) - dx; tmp_x(1:end-1)], tmp_y, m, [3, 0]);
[~, ~, py_ghost2] = param_generator([tmp_x(2:end); tmp_x(end) + dx], tmp_y, m, [4, 0]);
[~, ~, py_ghost_t] = param_generator(tmp_x, [2*tmp_y(1) - tmp_y(2); tmp_y(1:end-1)], m, [0, 3]);
[~, ~, py_ghost_b] = param_generator(tmp_x, [tmp_y(2:end); 2*tmp_y(end) - tmp_y(end-1)], m, [0, 4]);


% calculate finite differences in x
diff_x1 = mesh.param_x - x_ghost1; % x - x(0:end-1)
diff_x2 = x_ghost2 - mesh.param_x; % x(2:end+1) - x
diff_y1 = mesh.param_y - y_ghost1;
diff_y2 = y_ghost2 - mesh.param_y;

% calculate parameters on cell edges

% generate cell vectors on the cell edges
[x_vedge, y_vedge] = param_generator(tmp_xe, tmp_y);
[x_hedge, y_hedge] = param_generator(tmp_x, tmp_ye);

% Interpolate parameters onto the edges
% first create expanded ghost point vectors for the vertical edges
[xv1, ~, mv1] = param_generator([tmp_x(1) - dx; tmp_x], tmp_y, m, [1, 0]);
[xv2, ~, mv2] = param_generator([tmp_x; tmp_x(end) + dx], tmp_y, m, [2, 0]);
[~, ~, pxv1] = param_generator([tmp_x(1) - dx; tmp_x], tmp_y, px, [1, 0]);
[~, ~, pxv2] = param_generator([tmp_x; tmp_x(end) + dx], tmp_y, px, [2, 0]);
[~, ~, pyv1] = param_generator([tmp_x(1) - dx; tmp_x], tmp_y, py, [1, 0]);
[~, ~, pyv2] = param_generator([tmp_x; tmp_x(end) + dx], tmp_y, py, [2, 0]);
% perform interpolation
m_vedge = mv2.*(x_vedge - xv1)./(xv2 - xv1) + mv1.*(xv2 - x_vedge)./(xv2 - xv1);
px_vedge = pxv2.*(x_vedge - xv1)./(xv2 - xv1) + pxv1.*(xv2 - x_vedge)./(xv2 - xv1);
py_vedge = pyv2.*(x_vedge - xv1)./(xv2 - xv1) + pyv1.*(xv2 - x_vedge)./(xv2 - xv1);


% for the horizontal edges
[~, yh1, mh1] = param_generator(tmp_x, [2*tmp_y(1) - tmp_y(2); tmp_y], m, [0, 1]);
[~, yh2, mh2] = param_generator(tmp_x, [tmp_y; 2*tmp_y(end) - tmp_y(end-1)], m, [0, 2]);
[~, ~, pxh1] = param_generator(tmp_x, [2*tmp_y(1) - tmp_y(2); tmp_y], px, [0, 1]);
[~, ~, pxh2] = param_generator(tmp_x, [tmp_y; 2*tmp_y(end) - tmp_y(end-1)], px, [0, 2]);
[~, ~, pyh1] = param_generator(tmp_x, [2*tmp_y(1) - tmp_y(2); tmp_y], py, [0, 1]);
[~, ~, pyh2] = param_generator(tmp_x, [tmp_y; 2*tmp_y(end) - tmp_y(end-1)], py, [0, 2]);

m_hedge = mh2.*(y_hedge - yh1)./(yh2 - yh1) + mh1.*(yh2 - y_hedge)./(yh2 - yh1);
px_hedge = pxh2.*(y_hedge - yh1)./(yh2 - yh1) + pxh1.*(yh2 - y_hedge)./(yh2 - yh1);
py_hedge = pyh2.*(y_hedge - yh1)./(yh2 - yh1) + pyh1.*(yh2 - y_hedge)./(yh2 - yh1);

% shifted edge vectors containing either end for use in calculations
% first shift the mesh
[xe_ghost1, ~] = param_generator(tmp_xe(1:end-1), tmp_y);
[xe_ghost2, ~] = param_generator(tmp_xe(2:end), tmp_y);
[~, ye_ghost1] = param_generator(tmp_x, tmp_ye(1:end-1));
[~, ye_ghost2] = param_generator(tmp_x, tmp_ye(2:end));

% then calculate new parameter values
px_vedge1 = px_vedge; px_vedge2 = px_vedge; py_vedge1 = py_vedge; py_vedge2 = py_vedge;
m_vedge1 = m_vedge; m_vedge2 = m_vedge;
px_vedge1(xL+1:xL+1:end) = []; py_vedge1(xL+1:xL+1:end) = []; m_vedge1(xL+1:xL+1:end) = [];
px_vedge2(1:xL+1:end) = []; py_vedge2(1:xL+1:end) = []; m_vedge2(1:xL+1:end) = [];

px_hedge1 = px_hedge; px_hedge2 = px_hedge; py_hedge1 = py_hedge; py_hedge2 = py_hedge;
m_hedge1 = m_hedge; m_hedge2 = m_hedge;
px_hedge1(end-xL+1:end) = []; py_hedge1(end-xL+1:end) = []; m_hedge1(end-xL+1:end) = [];
px_hedge2(1:xL) = []; py_hedge2(1:xL) = []; m_hedge2(1:xL) = [];

% half point differences at cell centres for later use
diff_x_cell = (xe_ghost2 - xe_ghost1); % x(i = 1/2:end+1/2) - x(-1/2:end-1/2)
diff_y_cell = (ye_ghost2 - ye_ghost1); % y(j = 1/2:end+1/2) - y(-1/2:end-1/2)


% find corner points from the edges

[hc_ghost1, ~, mhc_ghost1] = param_generator([tmp_x(1) - dx; tmp_x], tmp_ye, m_hedge, [1, 0]); % x, m(0:end, -0.5:end+0.5)
[hc_ghost2, ~, mhc_ghost2] = param_generator([tmp_x; tmp_x(end) + dx], tmp_ye, m_hedge, [2, 0]); % x, m(1:end+1, -0.5:end+0.5)
[~, ~, pxhc_ghost1] = param_generator([tmp_x(1) - dx; tmp_x], tmp_ye, px_hedge, [1, 0]); % x, m(0:end, -0.5:end+0.5)
[~, ~, pxhc_ghost2] = param_generator([tmp_x; tmp_x(end) + dx], tmp_ye, px_hedge, [2, 0]); % x, m(1:end+1, -0.5:end+0.5)
[~, ~, pyhc_ghost1] = param_generator([tmp_x(1) - dx; tmp_x], tmp_ye, py_hedge, [1, 0]); % x, m(0:end, -0.5:end+0.5)
[~, ~, pyhc_ghost2] = param_generator([tmp_x; tmp_x(end) + dx], tmp_ye, py_hedge, [2, 0]); % x, m(1:end+1, -0.5:end+0.5)

[~, vc_ghost1, mvc_ghost1] = param_generator(tmp_xe, [2*tmp_y(1) - tmp_y(2); tmp_y], m_vedge, [0, 1]); % y, m(-0.5:end+0.5, 0:end)
[~, vc_ghost2, mvc_ghost2] = param_generator(tmp_xe, [tmp_y; 2*tmp_y(end) - tmp_y(end-1)], m_vedge, [0, 2]); % y, m(-0.5:end+0.5, 1:end+1)
[~, ~, pxvc_ghost1] = param_generator(tmp_xe, [2*tmp_y(1) - tmp_y(2); tmp_y], px_vedge, [0, 1]); % y, m(-0.5:end+0.5, 0:end)
[~, ~, pxvc_ghost2] = param_generator(tmp_xe, [tmp_y; 2*tmp_y(end) - tmp_y(end-1)], px_vedge, [0, 2]); % y, m(-0.5:end+0.5, 1:end+1)
[~, ~, pyvc_ghost1] = param_generator(tmp_xe, [2*tmp_y(1) - tmp_y(2); tmp_y], py_vedge, [0, 1]); % y, m(-0.5:end+0.5, 0:end)
[~, ~, pyvc_ghost2] = param_generator(tmp_xe, [tmp_y; 2*tmp_y(end) - tmp_y(end-1)], py_vedge, [0, 2]); % y, m(-0.5:end+0.5, 1:end+1)


% interpolate x,y onto corners
[xc, yc] = param_generator(tmp_xe, tmp_ye);% % x, y(-0.5:end+0.5, -0.5:end+0.5)
% xc = (hc_ghost1 + hc_ghost2)/2; 
% yc = (vc_ghost1 + vc_ghost2)/2; 

% Separately, interpolate the adjacent x and y points onto the corners
mhc = mhc_ghost2.*((xc - hc_ghost1)./(hc_ghost2 - hc_ghost1)) + mhc_ghost1.*((hc_ghost2 - xc)./(hc_ghost2 - hc_ghost1)); 
pxhc = pxhc_ghost2.*((xc - hc_ghost1)./(hc_ghost2 - hc_ghost1)) + pxhc_ghost1.*((hc_ghost2 - xc)./(hc_ghost2 - hc_ghost1));
pyhc = pyhc_ghost2.*((xc - hc_ghost1)./(hc_ghost2 - hc_ghost1)) + pyhc_ghost1.*((hc_ghost2 - xc)./(hc_ghost2 - hc_ghost1)); 

mvc = mvc_ghost2.*((yc - vc_ghost1)./(vc_ghost2 - vc_ghost1)) + mvc_ghost1.*((vc_ghost2 - yc)./(vc_ghost2 - vc_ghost1)); 
pxvc = pxvc_ghost2.*((yc - vc_ghost1)./(vc_ghost2 - vc_ghost1)) + pxvc_ghost1.*((vc_ghost2 - yc)./(vc_ghost2 - vc_ghost1)); 
pyvc = pyvc_ghost2.*((yc - vc_ghost1)./(vc_ghost2 - vc_ghost1)) + pyvc_ghost1.*((vc_ghost2 - yc)./(vc_ghost2 - vc_ghost1)); 


% take the mean of the x and y interpolations to get the corner point value
mc = (mvc + mhc)./2;
pxc = (pxvc + pxhc)./2;
pyc = (pyvc + pyhc)./2;

% Now form the corner point vectors needed for ctc. 1/2 is at (i-0.5:end-0.5/i+0.5:end+0.5). (t/b) (top/bottom) are (j=-1/2:end-1/2) and (j=1/2:end+1/2).
% preallocate
mc_t1 = mc; pxc_t1 = pxc; pyc_t1 = pyc;
mc_b1 = mc; pxc_b1 = pxc; pyc_b1 = pyc;
mc_t2 = mc; pxc_t2 = pxc; pyc_t2 = pyc;
mc_b2 = mc; pxc_b2 = pxc; pyc_b2 = pyc;

% remove bottom/top rows respectively
mc_t1(end - xL:end) = []; pxc_t1(end - xL:end) = []; pyc_t1(end - xL:end) = [];
mc_t2(end - xL:end) = []; pxc_t2(end - xL:end) = []; pyc_t2(end - xL:end) = [];
mc_b1(1:xL+1) = []; pxc_b1(1:xL+1) = []; pyc_b1(1:xL+1) = [];
mc_b2(1:xL+1) = []; pxc_b2(1:xL+1) = []; pyc_b2(1:xL+1) = [];

% remove l/r edge points
mc_t1(xL+1:xL+1:end) = []; pxc_t1(xL+1:xL+1:end) = []; pyc_t1(xL+1:xL+1:end) = [];
mc_b1(xL+1:xL+1:end) = []; pxc_b1(xL+1:xL+1:end) = []; pyc_b1(xL+1:xL+1:end) = [];
mc_t2(1:xL+1:end) = []; pxc_t2(1:xL+1:end) = []; pyc_t2(1:xL+1:end) = [];
mc_b2(1:xL+1:end) = []; pxc_b2(1:xL+1:end) = []; pyc_b2(1:xL+1:end) = [];
% ######################################################################
% Now prepare p update

% Rm calculated at cell centres: want in diagonal form
Rm = spdiags(1./sqrt(((m_vedge2 - m_vedge1)/dx - px).^2 + ((m_hedge2 - m_hedge1)/dx - py).^2 + gamma_c.^2), 0, xL*yL, xL*yL);

Rmx1 = 1./sqrt(((m - mx_ghost1)/dx - px_vedge1).^2 + ((mc_b1 - mc_t1)/dx - py_vedge1).^2 + gamma_c.^2);
Rmx2 = 1./sqrt(((mx_ghost2 - m)/dx - px_vedge2).^2 + ((mc_b2 - mc_t2)/dx - py_vedge2).^2 + gamma_c.^2);
Rmy1 = 1./sqrt(((mc_t2 - mc_t1)/dx - px_hedge1).^2 + ((m - my_ghost1)/dx - py_hedge1).^2 + gamma_c^2);
Rmy2 = 1./sqrt(((mc_b2 - mc_b1)/dx - px_hedge2).^2 + ((my_ghost2 - m)/dx - py_hedge2).^2 + gamma_c.^2);


% Rpx on vertical edges
Rpx1 = sqrt(((px - px_ghost1)/dx).^2 + ((pyc_b1 - pyc_t1)/dx).^2 + + ((py - py_ghost1)/dx + (pxc_b1 - pxc_t1)/dx).^2 + gamma_p.^2);
Rpx2 = sqrt(((px_ghost2 - px)/dx).^2 + ((pyc_b2 - pyc_t2)/dx).^2 + + ((py_ghost2 - py)/dx + (pxc_b2 - pxc_b1)/dx).^2 + gamma_p.^2);
% Rpy on horizontal edges
Rpy_t = sqrt(((pxc_t2 - pxc_t1)/dx).^2 + ((py - py_ghost_t)/dx).^2 + + ((pyc_t2 - pyc_t1)/dx + (px - px_ghost_t)/dx).^2 + gamma_p.^2);
Rpy_b = sqrt(((pxc_b2 - pxc_b1)/dx).^2 + ((py_ghost_b - py)/dx).^2 + + ((pyc_b2 - pyc_b1)/dx + (px_ghost_b - px)/dx).^2 + gamma_p.^2);

% % convert to sparse matrix form
% Rpx1 = sparse(1:xL*yL, 1:xL*yL, Rpx1, xL*yL, xL*yL);
% Rpx2 = sparse(1:xL*yL, 1:xL*yL, Rpx2, xL*yL, xL*yL);
% Rpy_t = sparse(1:xL*yL, 1:xL*yL, Rpy_t, xL*yL, xL*yL);
% Rpy_b = sparse(1:xL*yL, 1:xL*yL, Rpy_b, xL*yL, xL*yL);

% now create matrices for CxRpCx,CyRpCy, CyRpCx, CxRpCy 

% CxRpCx
xrx0 = Rpx1/dx.^2; % x(i-1)
xrx1 = - (Rpx2 + Rpx1)/dx.^2; % x(i)
xrx2 = Rpx2/dx.^2; % x(i+1)
% Make ghost points implicit.
xrx1(1:xL:end) = xrx1(1:xL:end) + xrx0(1:xL:end); % apply x = (0) bc implicitly - as x(0) = -x(1)
xrx1(xL:xL:end) = xrx1(xL:xL:end) + xrx2(xL:xL:end); % apply x = xL+1 bc implicitly
xrx0(1:xL:end) = 0; % set start ghost points to zero
xrx2(xL:xL:end) = 0; % set end ghost points to zero
% form fd matrix
cxRcx = sparse(1:xL*yL, 1:xL*yL, xrx1, xL*yL, xL*yL);
cxRcx = spdiags([xrx0, xrx2], [1; -1], cxRcx)';

% CyRpCy
yry0 =  Rpy_t/dx.^2; % y(j-1)
yry1 = - (Rpy_b + Rpy_t)/dx.^2; %(j)
yry2 = Rpy_b/dx.^2; % y(j+1)
% Make ghost points implicit.
yry1(1:xL) = yry1(1:xL) + yry0(1:xL);
yry1(end-xL+1:end) = yry1(end-xL+1:end) + yry2(end-xL+1:end);
% form fd matrix
cyRcy = spdiags([yry1, yry0, yry2], [0, xL, -xL], xL*yL, xL*yL)'; % adjacent terms. behaviour of spdiag used to exclude ghost points


% The mixed Cx Cy terms are made up of the interpolated corner points, and
% therefore take some care to assemble.

% I will use the interpolation p{i+1/2,j+1/2} = 0.25*(p{i,j) + p{i+1,j} +
% p{i, j+1} + p{i+1,j+1}) in the matrices, ehich assumes uniform cell sizes

% [CyRpCx*p]{i,j} = Rpx2*(p{i+1/2, j+1/2} - p{i+1/2, j-1/2})/(dxe*dyc) - Rpx1*(p{i-1/2, j+1/2} - p{i-1/2, j-1/2})/(dxe*dyc) 

% construct ctc_x, ctc_y matrices
% Find contributing terms to cxT*Rcx*Cx and cyTRcyCy, in terms of
% cell-centre parameters. _pq gives x and y coordinates of cell, where p/q=
% x,y dimension p/q = 0,1,2 are i/j -1,+0,+1

yrx_00 = Rpx1/4/dx.^2; % {i-1, j-1)
yrx_22 = Rpx2/4/dx.^2; % (i+1, j+1)
yrx_10 = (- Rpx2 + Rpx1)/4/dx.^2; % (i, j-1)
yrx_12 = (Rpx2 - Rpx1)/4/dx.^2; % (i, j+1)
yrx_02 = -Rpx1/4/dx.^2; % (i-1, j+1)
yrx_20 = -Rpx2/4/dx.^2; % (i+1, j-1)
% Make ghost point corrections

% j boundaries 
% (i+1, j) diagonal 
yrx_bc2 = sparse(xL*yL, 1); 
yrx_bc2(1:xL) = yrx_20(1:xL);
yrx_bc2(end-xL+1:end) = yrx_22(end-xL+1:xL:end-1); % avoids doubel counting (i+1,j+1 point) at corner
% (i-1, j) diagonal
yrx_bc0 = sparse(xL*yL, 1); 
yrx_bc0(2:xL) = yrx_00(2:xL); % avoids double counting (i-1,j-1 point) at corner
yrx_bc0(end-xL+1:end) = yrx_02(end-xL+1:end);
% (i,j) diagonal
yrx_bc1 = sparse(xL*yL, 1); 
yrx_bc1(1:xL) = yrx_10(1:xL);
yrx_bc1(end-xL+1:end) = yrx_12(end-xL+1:end);
% x boundaries
% (j-1, 0); (j-1, end+1) 
yrx_10(1:xL:end) = yrx_10(1:xL:end) + yrx_00(1:xL:end); 
yrx_10(xL:xL:end) = yrx_10(xL:xL:end) + yrx_20(xL:xL:end);
yrx_00(1:xL:end) = 0;
yrx_20(xL:xL:end) = 0;
% (j+1, 0); (j+1, end+1) 
yrx_12(1:xL:end) = yrx_12(1:xL:end) + yrx_02(xL:xL:end); 
yrx_12(xL:xL:end) = yrx_12(xL:xL:end) + yrx_22(xL:xL:end); 
yrx_02(1:xL:end) = 0;
yrx_22(xL:xL:end) = 0;

% now constuct matrix
cyRcx = spdiags([yrx_10, yrx_12, yrx_00, yrx_22, yrx_02, yrx_20, yrx_bc1, yrx_bc0, yrx_bc2], [xL, -xL, 1 + xL, -1 - xL, -xL + 1, xL - 1, 0, 1, -1], xL*yL, xL*yL)';
                
% CxRpCy
xry_22 = Rpy_b/4/dx.^2;
xry_00 = Rpy_t/4/dx.^2;
xry_21 = (Rpy_b - Rpy_t)/4/dx.^2;
xry_01 = (-Rpy_b + Rpy_t)/4/dx.^2;
xry_20 = -Rpy_t/4/dx^2;
xry_02 = -Rpy_b/4/dx^2;

% j boundaries 
% (i, j+1) diagonal 
xry_bc2 = sparse(xL*yL, 1); 
xry_bc2(1:xL) = xry_02(1:xL);
xry_bc2(end-xL+1:end) = xry_22(end-xL+1:xL:end-1); % avoids doubel counting (i+1,j+1 point) at corner
% (i, j-1) diagonal
xry_bc0 = sparse(xL*yL, 1); 
xry_bc0(2:xL) = xry_00(2:xL); % avoids double counting (i-1,j-1 point) at corner
xry_bc0(end-xL+1:end) = xry_20(end-xL+1:end);
% (i,j) diagonal
xry_bc1 = sparse(xL*yL, 1); 
xry_bc1(1:xL) = xry_01(1:xL);
xry_bc1(end-xL+1:end) = xry_21(end-xL+1:end);
% x boundaries
% (j-1, 0); (j-1, end+1) 
xry_01(1:xL:end) = xry_01(1:xL:end) + xry_00(1:xL:end); 
xry_01(xL:xL:end) = xry_01(xL:xL:end) + xry_02(xL:xL:end);
xry_00(1:xL:end) = 0;
xry_02(xL:xL:end) = 0;
% (j+1, 0); (j+1, end+1) 
xry_21(1:xL:end) = xry_21(1:xL:end) + xry_20(xL:xL:end); 
xry_21(xL:xL:end) = xry_21(xL:xL:end) + xry_22(xL:xL:end); 
xry_20(1:xL:end) = 0;
xry_22(xL:xL:end) = 0;

cxRcy = spdiags([xry_00, xry_22, xry_21, xry_01, xry_20, xry_02, xry_bc1, xry_bc0, xry_bc2], [1+xL, -1-xL, -1, 1, -1+xL, 1-xL, 0, xL, -xL], xL*yL, xL*yL)';
            

% invert signs
% cxRcx = -cxRcx;
% cxRcy = -cxRcy;
% cyRcx = -cyRcx;
% cyRcy = -cyRcy;

% Now solve for p
a11 = Rm + tgv_lagrn*(cxRcx + 0.5*cyRcy);
a12 = tgv_lagrn*0.5*cyRcx;
a21 = tgv_lagrn*0.5*cxRcy;
a22 = Rm + tgv_lagrn*(cyRcy + 0.5*cxRcx);
b1 = (Rmx2.*m_vedge2 - Rmx1.*m_vedge1)/dx;
b2 = (Rmy2.*m_hedge2 - Rmy1.*m_hedge1)/dx;

A = [a11, a12; a21, a22];
b = [b1; b2];

p2 = A\b;
end