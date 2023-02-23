function [CxtRcCx, CytRcCy, CxtRc, CytRc] = ctc_calc2(mesh, input, px, py, gamma, itr)
% m0 = m1, m(-1) = m2

test_flag = 0;

% key values
m = log10(mesh.res_param1);
tmp_x = unique(mesh.param_x);
tmp_y = unique(mesh.param_y);
xL = length(tmp_x);
yL = length(tmp_y);
dx = tmp_x(2) - tmp_x(1); % x cell width
dy = tmp_y(2) - tmp_y(1); % y cell width
dz = 1; % factor of increase in y depth with each layer

% dx = 1;
% dy=1;
%     m = repmat([1:xL]',yL,1);

% for testing
if test_flag == 1
    xL = 5;yL = 5;
    dx = 1; dy=1; tmp_x = 1:4; tmp_y = tmp_x;
    
%     m = [1:xL*yL]';
    %     m = [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4]';
    %     m = [1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]';
    m= [11:15,21:25,31:35,41:45,51:55]';
%     px = 0*m;
%     py = 0*m;
        px = m./10;
        py = m./10;

end

% first calculate R terms
% cell centred  Rcx
rcx = (circshift(m, -1) - circshift(m, 1))./(2*dx);
rcx(1:xL:end) = (m(2:xL:end) - m(1:xL:end))./(2*dx); % x(0) boundary
rcx(xL:xL:end) = (m(xL:xL:end) - m(xL-1:xL:end))./(2*dx);% x(xL) boundary

rcy = (circshift(m, -xL) - circshift(m, xL))./(2*dy);
rcy(1:xL) = (m(xL+1:2*xL) - m(1:xL))./(2*dy); % y(0) boundary
rcy(end-xL+1:end) = (m(end-xL+1:end) - m(end-2*xL+1:end-xL))./(2*dy); % y(end) boundary

Rc = 1./sqrt((rcx - px).^2 + (rcy - py).^2 + gamma.^2); % rcx/rcy do not include p here

% Rcx: Rc on vertical edges (-1.5:1:xL+1.5, j)

% pad with extra cells 
mve = zeros((xL+1)*yL, 1);
pxve = mve; 
pyve = mve;
for j = 1:yL
   mve((j-1)*(xL+1) + (1:xL+1)) = [m((j-1)*xL + (1:xL)); 0];
   pxve((j-1)*(xL+1) + (1:xL+1)) = [px((j-1)*xL + (1:xL)); 0];
   pyve((j-1)*(xL+1) + (1:xL+1)) = [py((j-1)*xL + (1:xL)); 0];
end

rcx = (mve - circshift(mve, 1))./dx - (pxve + circshift(pxve, 1))./2; % i =1.5:xL-0.5 (interior points)
rcx(1:xL+1:end) = 0; % i = 0.5
rcx(xL+1:xL+1:end) = 0; % i = xL + 0.5

rcy = (circshift(mve, -(xL+1)) + circshift(mve, -(xL+1)+1) - circshift(mve, (xL+1)) - circshift(mve, (xL+1)+1))./(4*dy) - (pyve + circshift(pyve, 1))./2; 
rcy(2:xL) = (m(xL + (1:xL-1)) + m(xL + (2:xL)) - m(1:xL-1) - m(2:xL))./(4*dy) - (py(1:xL-1) + py(2:xL))./2; % j = 1, i = 1.5:xL-0.5
rcy([1,xL+1]) = (m([1, xL] + xL) - m([1, xL]))/(2*dy); % i = [0.5, xL+0.5], j = 1
rcy(xL+2:xL+1:end-xL-1) = (m((xL+1:xL:end-xL) + xL) - m((xL+1:xL:end-xL) - xL))./(2*dy); % i = -0.5, j=2:yL-1
rcy(2*(xL+1):xL+1:end-xL-1) = (m((2*xL:xL:end-xL) + xL) - m((2*xL:xL:end-xL) - xL))./(2*dy); % i = xL+0.5, j = 2:yL-1
rcy(end-xL+1:end-1) = (m(end - xL + (1:xL-1)) + m(end - xL + (2:xL)) - m(end - 2*xL +(1:xL-1)) - m(end - 2*xL + (2:xL)))./(4*dy) - (py(end - xL +(1:xL-1)) + py(end - xL + (2:xL)))./2; % j = yL, i = 0.5:xL-0.5
rcy([end-xL, end]) = (m([end-xL+1, end]) - m([end-xL+1, end] - xL))/(2*dy); % i = [0.5, xL+0.5], j = yL

Rcx = 1./sqrt(rcx.^2 + rcy.^2 + gamma.^2); 

% Rcy: Rc on horizontal edges
rcx = (circshift(m, -1) + circshift(m, -1+xL) - circshift(m, 1) - circshift(m, 1+xL) )./(4*dx) - (px + circshift(px, +xL))./2; % j = 1.5:yL-0.5, i = 2:xL-1. first row is j = 0.5, 
rcx(1:xL) = [(m(2)-m(1))./(2*dx); (m(3:xL) - m(1:xL-2))./(2*dx); (m(xL) - m(xL-1))./(2*dx)]; % j = 0.5 row.
rcx(xL+1:xL:end) = (m(2:xL:end-xL) + m(xL+2:xL:end) - m(1:xL:end-xL) - m(xL+1:xL:end))./(4*dx) - (px(1:xL:end-xL) + px(xL+1:xL:end))./2; % i = 1, j = 1.5:yL-0.5
rcx(2*xL:xL:end) = (m(xL:xL:end-xL) + m(2*xL:xL:end) - m(xL-1:xL:end-xL) - m(2*xL-1:xL:end))./(4*dx) - (px(xL:xL:end-xL) + px(2*xL:xL:end))./2; % i = xL, j = 1.5:yL-0.5
rcx = [rcx; (m(end-xL+2) - m(end-xL+1))./(2*dx); (m(end-xL+3:end) - m(end-xL+1:end-2))./(2*dx); (m(end) - m(end-1))./(2*dx)]; % add j = yL + 0.5 row

rcy = (circshift(m, -xL) - m)./dy - (circshift(py, -xL) + py)./2; % j = 1.5:yL-0.5, i = 1:xL. 
rcy = [zeros(xL,1) ;rcy]; % j = 0.5
rcy(end-xL+1:end) = 0; % j = yL + 0.5

Rcy = 1./sqrt(rcx.^2 + rcy.^2 + gamma.^2); 

if itr == 1 % does l2 iteration instead
    Rcx = ones(size(Rcx));
    Rcy = ones(size(Rcy));
end

% now construct inversion matrices
% 
% % test
if test_flag == 1
    Rc = [1:xL*yL]';
    Rcx = repmat([0.5:1:xL+0.5]', yL, 1);
    for i = 1:yL+1
        Rcy(1+(i-1)*xL:i*xL) = i - 0.5;
    end
end

% Cx^T*Rc
% px(i - 1, j)
xr0 = Rc./(2*dx); 
xr0(1:xL:end) = 0;
% px(i + 1, j)
xr2 = -Rc./(2*dx); 
xr2(xL:xL:end) = 0;
% applies bc at i = 1, xL
xr1 = sparse(xL*yL, 1); 
xr1(1:xL:end) = Rc(1:xL:end)./(2*dx); % i = 1 bc
xr1(xL:xL:end) = -Rc(xL:xL:end)./(2*dx); % i = xL bc

% px2 = px2.*0;
% form Cx^T*Rc matrix
CxtRc = spdiags([xr0, xr1, xr2], [1, 0, -1], xL*yL, xL*yL)';

% Cy^T*Rc
% py(i, j-1) 
yr = Rc./(2*dy); % off diagonal elements are equal magnitude
yr1 = sparse(xL*yL, 1);
yr1(1:xL) = Rc(1:xL)./(2*dy);
yr1(end-xL+1:end) = -Rc(end-xL+1:end)./(2*dy);

CytRc = spdiags([yr, yr1, -yr], [xL, 0, -xL], xL*yL, xL*yL)';


% Cx^T*Rc*Cx
% main diagonal: m(i,j)
xrx2 = circshift(Rcx, -1) + Rcx; 
xrx2(xL+1:xL+1:end) = [];
xrx2(1:xL:end) = 2*Rcx(1:xL+1:end);% bc at i = 1
xrx2(xL:xL:end) = 2*Rcx(xL+1:xL+1:end); % bc at i = xL

% m(i-1,j)
xrx1 = - circshift(Rcx, 1) + Rcx;
xrx1(xL+1:xL+1:end) = [];
xrx1(1:xL:end) = 0; % bc at i = 1
xrx1(2:xL:end) = -2*Rcx(1:xL+1:end) + Rcx(2:xL+1:end); % bc at i = 2
xrx1(xL:xL:end) = - Rcx(xL-1:xL+1:end); % bc at i = xL

% m(i-2, j)
xrx0 = - circshift(Rcx, 1);
xrx0(xL+1:xL+1:end) = [];
xrx0([1:xL:end, 2:xL:end]) = 0; % bc at i = 1, i = 2 

% m(i+1,j)
xrx3 = circshift(Rcx, -1) - circshift(Rcx, -2);
xrx3(xL+1:xL+1:end) = [];
xrx3(xL:xL:end) = 0; % bc at i = xL
xrx3(1:xL:end) =  - Rcx(3:xL+1:end); % bc at i = 1
xrx3(xL-1:xL:end) = -2*Rcx(xL+1:xL+1:end) + Rcx(xL:xL+1:end); % bc at i = xL-1

% m(i+2,j)
xrx4 = -circshift(Rcx, -2);
xrx4(xL+1:xL+1:end) = [];
xrx4([xL:xL:end, xL-1:xL:end]) = 0; % bc at i = xL, xL-1

CxtRcCx = spdiags([xrx0, xrx1, xrx2, xrx3, xrx4],[2, 1, 0, -1, -2], xL*yL, xL*yL)'./(4*dx.^2);

% Cy^T*Rc*Cy
% main diagonal, m(i,j)
yry2 = circshift(Rcy, -xL) + Rcy;
yry2(end-xL+1:end) = []; % remove extra values
yry2(1:xL) = 2*Rcy(1:xL); % bc at j = 1
yry2(end-xL+1:end) = 2*Rcy(end-xL+1:end);% bc at j = yL

% m(i, j-1)
yry1 = Rcy - circshift(Rcy, xL);
yry1(end-xL+1:end) = []; % remove extra values
yry1(xL+1:2*xL) = -2*Rcy(1:xL) + Rcy(2*xL+1:3*xL); % bc at j = 2
yry1(end-xL+1:end) = -Rcy(end-3*xL+1:end-2*xL); % bc at j = yL

% m(i, j-2)
yry0 = - circshift(Rcy, xL);
yry0(end-xL+1:end) = []; % remove extra values - other boundary points  removed by spdiags

% m(i, j+1)
yry3 = circshift(Rcy, -xL) - circshift(Rcy, -2*xL);
yry3(end-xL+1:end) = []; % remove extra values
yry3(end-2*xL+1:end-xL) = Rcy(end-2*xL+1:end-xL) - 2*Rcy(end-xL+1:end); % bc at j = yL-1
yry3(1:xL) = -Rcy(2*xL+1:3*xL);% bc at j = 1

% m(i, j+2)
yry4 = -circshift(Rcy, -2*xL);
yry4(end-xL+1:end) = []; % remove extra values - spdiags removes other boundary points

CytRcCy = spdiags([yry0, yry1, yry2, yry3, yry4],[2*xL, xL, 0, -xL, -2*xL], xL*yL, xL*yL)'./(4*dy.^2);




end