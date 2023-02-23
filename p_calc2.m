% Performs half point-difference calculation
% p2 = [px; py]
% m - log10(resistivty), mu - tgv_lagrn
function [p2] = p_calc2(mesh, lagrn ,tgv_lagrn, gamma_c, gamma_p)

px = mesh.px;
py = mesh.py;

% key values
m = log10(mesh.res_param1);
tmp_x = unique(mesh.param_x);
tmp_y = unique(mesh.param_y);
xL = length(tmp_x);
yL = length(tmp_y);
dx = tmp_x(2) - tmp_x(1); % x cell width
dy = tmp_y(2) - tmp_y(1); % y cell width
dz = 1; % factor of increase in y depth with each layer

if test_flag == 1
    %     m = [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4]';
    %     m = [1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]';
    m= [11:15,21:25,31:35,41:45,51:55]';
    px = 0*m;
    py = 0*m;
    %     px = m./10;
    %     py = m./10;
    xL = 5;yL = 5;
    dx = 1; dy=1; tmp_x = 1:4; tmp_y = tmp_x;
end

% Calculate Rm (=Rc) at cell centres

% Rpx on vertical edges

% Rpx on horizontal edges


% Calculate CRpCt terms

% Cx*Rp*Cxt

% Cy*Rp*Cyt

% Cy*Rp*Cxt

% Cx*Rp*Cyt






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