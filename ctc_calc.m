% [ctc_x, ctc_y] = ctc_calc(res, px, py, param_x, param_y)
%  Calculates terms in inverse problem using half point central diferences
function [ctc_x, ctc_y, ctrc_px, ctrc_py] = ctc_calc(mesh, input, px, py, gamma, itr)

% key values
m = log10(mesh.res_param1);
tmp_x = unique(mesh.param_x);
tmp_y = unique(mesh.param_y);
xL = length(tmp_x);
yL = length(tmp_y);
dx = (tmp_x(2) - tmp_x(1)); % cell width of x edge cells
dy = tmp_y(2) - tmp_y(1);
tmp_xe = mesh.xe;
tmp_ye = mesh.ye;

% testing
% m = repmat([1:xL]',yL,1);

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

% Now form the corner point vectors needed for ctc. 1/2 is at (i-0.5:end-0.5 / i+0.5:end+0.5). (t/b) (top/bottom) are (j=-0.5:end-0.5) and (j=0.5:end+0.5).
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


if itr == 1 || input.inv_flag == -2 % L2 iteration to start the inversion without Rc terms
    Rcx1 = ones(size(m));
    Rcy_t = Rcx1;
    Rcx2 = Rcx1;
    Rcy_b = Rcx1;
else
    Rcx1 = 1./sqrt(( (m - mx_ghost1)./dx - px_vedge1).^2 + ((mc_b1 - mc_t1)./dy - py_vedge1).^2 + gamma.^2);
    Rcy_t = 1./sqrt(((mc_t2 - mc_t1)/dx - px_hedge1).^2 + ((m - my_ghost1)/dy - py_hedge1).^2 + gamma^2);
    Rcx2 = 1./sqrt(((mx_ghost2 - m)/dx - px_vedge2).^2 + ((mc_b2 - mc_t2)/dy - py_vedge2).^2 + gamma.^2);
    Rcy_b = 1./sqrt(((mc_b2 - mc_b1)/dx - px_hedge2).^2 + ((my_ghost2 - m)/dy - py_hedge2).^2 + gamma.^2);
end

% construct ctc_x, ctc_y matrices
% Find contributing terms to cxT*Rcx*Cx and cyTRcyCy.
ctc_xi0 = Rcx1/dx.^2; % x(i-1)
ctc_xi1 = - (Rcx2 + Rcx1)/dx.^2; % x(i)
ctc_xi2 = Rcx2/dx.^2; % x(i+1)
ctc_yi0 =  Rcy_t/dx.^2; % y(j-1)
ctc_yi1 = - (Rcy_b + Rcy_t )/dx.^2; %(j)
ctc_yi2 = Rcy_b/dx.^2; % y(j+1)

% Make ghost points implicit.
ctc_xi1(1:xL:end) = ctc_xi1(1:xL:end) + ctc_xi0(1:xL:end); % apply x = (0) bc implicitly - as x(0) = -x(1)
ctc_xi1(xL:xL:end) = ctc_xi1(xL:xL:end) + ctc_xi2(xL:xL:end); % apply x = xL+1 bc implicitly
ctc_yi1(1:xL) = ctc_yi1(1:xL) + ctc_yi0(1:xL);
ctc_yi1(end-xL+1:end) = ctc_yi1(end-xL+1:end) + ctc_yi2(end-xL+1:end);

ctc_x = sparse(1:xL*yL, 1:xL*yL, ctc_xi1, xL*yL, xL*yL);
ctc_xi0(1:xL:end) = 0; % set start ghost points to zero
ctc_xi2(xL:xL:end) = 0; % set end ghost points to zero
ctc_x = spdiags([ctc_xi0, ctc_xi2], [1; -1], ctc_x)';

ctc_y = sparse(1:xL*yL, 1:xL*yL, ctc_yi1, xL*yL, xL*yL);
ctc_y = spdiags([ctc_yi0, ctc_yi2], [xL, -xL], ctc_y)';

% next calculate cxTRc and cyTRc terms
if itr == 1 || input.inv_flag == -2
    Rcp = ones(size(px));
else
    Rcp = 1./sqrt(((m_vedge2 - m_vedge1)/dx - px).^2 + ((m_hedge2 - m_hedge1)/dy - py).^2 + gamma.^2);
end

ctrc_px = (py_vedge2 - py_vedge1)/dx.*Rcp;
ctrc_py = (py_hedge2 - py_hedge1)/dx.*Rcp;
end



















