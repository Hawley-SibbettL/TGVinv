function [ctc_x, ctc_y, cxRc, cyRc] = ctc_calc3(mesh, input, px, py, gamma, itr)
% m0 = m2, m(-1) = m3

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
if test_flag == 1 || test_flag == 2
    xL = 5;yL = 5;
    dx = 1; dy=1; tmp_x = 1:4; tmp_y = tmp_x;
    
%     m = [1:xL*yL]';
    %     m = [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4]';
    %     m = [1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]';
    m= [11:15,21:25,31:35,41:45,51:55]';
    px = 0*m;
    py = 0*m;
%         px = m./10;
%         py = m./10;

end

% first calculate R terms
% cell centred  Rcx
rcx = (circshift(m, -1) - circshift(m, 1))./(2*dx) - px;
rcx(1:xL:end) = - px(1:xL:end); % x(0) boundary
rcx(xL:xL:end) =  - px(xL:xL:end);% x(xL) boundary

rcy = (circshift(m, -xL) - circshift(m, xL))./(2*dy) - py;
rcy(1:xL) = - py(1:xL); % y(0) boundary
rcy(end-xL+1:end) = - py(end-xL+1:end); % y(end) boundary

Rc = 1./sqrt(rcx.^2 + rcy.^2 + gamma.^2); % rcx/rcy do not include p here

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
rcx(1:xL+1:end) = (m(1:xL:end) - m(2:xL:end))./dx - (px(1:xL:end) - px(2:xL:end))./2; % i = 0.5
rcx(xL+1:xL+1:end) = (m(xL-1:xL:end) - m(xL:xL:end))./dx - (px(xL:xL:end) - px(xL-1:xL:end))./2; % i = xL + 0.5

rcy = (circshift(mve, -(xL+1)) + circshift(mve, -(xL+1)+1) - circshift(mve, (xL+1)) - circshift(mve, (xL+1)+1))./(4*dy) - (pyve + circshift(pyve, 1))./2; 
rcy(2:xL) = - (py(2:xL) - py(1:xL-1))/2; % j = 1, i = 1.5:xL-0.5
rcy([1,xL+1]) = - (py([1, xL]) - py([2, xL-1]))/2; % i = [0.5, xL+0.5], j = 1
rcy(end-xL+1:end-1) = - (py(end-xL+1:end-1) + py(end-xL+2:end))/2; % j = yL, i = 0.5:xL-0.5
rcy([end-xL, end]) = - (py([end-xL+1, end]) - py([end-xL+2, end-1]))./2; % i = [0.5, xL+0.5], j = yL
% rcy(xL+2:xL+1:end-xL-1) = (m((xL+1:xL:end-xL) + xL) + m((xL+2:xL:end-xL) + xL) - m((xL+1:xL:end-xL) - xL) - m((xL+2:xL:end-xL) - xL))./(4*dy) - (py(xL+1:xL:end-xL) - py(xL+2:xL:end-xL))/2; % i = -0.5, j=2:yL-1
% rcy(2*(xL+1):xL+1:end-xL-1) = (m((2*xL:xL:end-xL) + xL) + m((2*xL-1:xL:end-xL) + xL) - m((2*xL:xL:end-xL) - xL) - m((2*xL-1:xL:end-xL) - xL))./(4*dy) - (py(2*xL:xL:end-xL) - py(2*xL-1:xL:end-xL))/2; % i = xL+0.5, j = 2:yL-1

Rcx = 1./sqrt(rcx.^2 + rcy.^2 + gamma.^2); 

% Rcy: Rc on horizontal edges
rcx = (circshift(m, -1) + circshift(m, -1+xL) - circshift(m, 1) - circshift(m, 1+xL) )./(4*dx) - (px + circshift(px, +xL))./2; % j = 1.5:yL-0.5, i = 2:xL-1. first row is j = 0.5, 
rcx(2:xL-1) = (m(3:xL) + m(xL + (3:xL)) - m(1:xL-2) - m(xL + (1:xL-2)))/(4*dx) - (py(2:xL-1) - py((2:xL-1) + xL))/2;% j = 0.5 row, i = 2:xL-1
rcx = [rcx; - (px(end-xL+1) - px(end-2*xL+1))/2; (m(end-xL+3:end) + m((end-xL+3:end) - xL) - m((end-xL+1:end-2) - xL) - m(end-xL+1:end-2))/(4*dx) - (px(end-xL+2:end-1) - px((end-xL+2:end-1) - xL))/2; - (px(end) - px(end-xL))/2]; % add j = yL + 0.5 row
rcx(1:xL:end) = [- (px(1) - px(1+xL))/2; - (px(1:xL:end-xL) + px(xL+1:xL:end))/2; - (px(xL) - px(2*xL))/2]; % i = 1, j = 0.5:xL+0.5
rcx(xL:xL:end-xL) = [- (px(xL) - px(2*xL))/2; - (px(xL:xL:end-xL) + px(2*xL:xL:end))/2]; % i = xL, j = 0.5:xL-0.5


rcy = (circshift(m, -xL) - m)./dy - (circshift(py, -xL) + py)./2; % j = 1.5:yL-0.5, i = 1:xL. 
rcy = [(m(1:xL) - m((1:xL) + xL))/dy - (py(1:xL) - py((1:xL) + xL))/2 ;rcy]; % j = 0.5
rcy(end-xL+1:end) = (m(end-xL+1:end) - m((end-xL+1:end) - xL))/dy - (py(end-xL+1:end) - py( (end-xL+1:end)- xL))/2; % j = yL + 0.5

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
xr0([1:xL:end, xL:xL:end]) = 0;
% px(i + 1, j)
xr2 = -Rc./(2*dx); 
xr2([1:xL:end, xL:xL:end]) = 0;

% px2 = px2.*0;
% form Cx^T*Rc matrix
cxRc = spdiags([xr0, xr2], [1, -1], xL*yL, xL*yL)';

% Cy^T*Rc
% py(i, j-1) 
yr0 = Rc./(2*dy); % off diagonal elements are equal magnitude
yr0([1:xL,end-xL+1:end]) = 0;
yr2 = -Rc./(2*dy); % off diagonal elements are equal magnitude
yr2([1:xL,end-xL+1:end]) = 0;

cyRc = spdiags([yr0, yr2], [xL, -xL], xL*yL, xL*yL)';

% Cx^T*Rc*Cx
% main diagonal: m(i,j)
xrx2 = circshift(Rcx, -1) + Rcx; 
xrx2(xL+1:xL+1:end) = []; % remove extra points
xrx2(2:xL:end) = Rcx(2:xL+1:end) + Rcx(3:xL+1:end) - Rcx(1:xL+1:end); % bc at i = 2
xrx2(xL-1:xL:end) = Rcx(xL-1:xL+1:end) + Rcx(xL:xL+1:end) - Rcx(xL+1:xL+1:end);% bc at i = xL-1

% m(i-1,j)
xrx1 = - circshift(Rcx, 1) + Rcx;
xrx1(xL+1:xL+1:end) = []; % remove extra points
xrx1(1:xL:end) = 0; % bc at i = 1
xrx1(xL:xL:end) = Rcx(xL+1:xL+1:end) + Rcx(xL:xL+1:end) - 2*Rcx(xL-1:xL+1:end); % bc at i = xL

% m(i-2, j)
xrx0 = - circshift(Rcx, 1);
xrx0(xL+1:xL+1:end) = []; % remove extra points
xrx0([1:xL:end, 2:xL:end]) = 0; % bc at i = 1, i = 2 
xrx0(xL:xL:end) = -2*Rcx(xL-1:xL+1:end); % bc at i = xL

% m(i+1,j)
xrx3 = circshift(Rcx, -1) - circshift(Rcx, -2); 
xrx3(xL+1:xL+1:end) = []; % remove extra points
xrx3(xL:xL:end) = 0; % bc at i = xL
xrx3(1:xL:end) = Rcx(1:xL+1:end) + Rcx(2:xL+1:end) - 2*Rcx(3:xL+1:end); % bc at i = 1

% m(i+2,j)
xrx4 = -circshift(Rcx, -2);
xrx4(xL+1:xL+1:end) = []; % remove extra points
xrx4([xL:xL:end, xL-1:xL:end]) = 0; % bc at i = xL, xL-1
xrx4(1:xL:end) = -2*Rcx(3:xL+1:end); % bc at i = 1

ctc_x = spdiags([xrx0, xrx1, xrx2, xrx3, xrx4],[2, 1, 0, -1, -2], xL*yL, xL*yL)'./(4*dx.^2);

% Cy^T*Rc*Cy
% main diagonal, m(i,j)
yry2 = circshift(Rcy, -xL) + Rcy;
yry2(end-xL+1:end) = []; % remove extra values
yry2(xL+1:2*xL) = Rcy(xL+1:2*xL) + Rcy(2*xL+1:3*xL) - Rcy(1:xL); % bc at j = 2
yry2(end-2*xL+1:end-xL) = Rcy(end-3*xL+1:end-2*xL) + Rcy(end-2*xL+1:end-xL) - Rcy(end-xL+1:end);% bc at j = yL-1

% m(i, j-1)
yry1 = Rcy - circshift(Rcy, xL);
yry1(end-xL+1:end) = []; % remove extra values
yry1(1:xL) = 0; % bc at j = 1
yry1(end-xL+1:end) = -2*Rcy(end-3*xL+1:end-2*xL) + Rcy(end-2*xL+1:end-xL) + Rcy(end-xL+1:end); % bc at j = yL

% m(i, j-2)
yry0 = - circshift(Rcy, xL);
yry0(end-xL+1:end) = []; % remove extra values - other boundary points  removed by spdiags
yry0(1:2*xL) = 0;
yry0(end-xL+1:end) = -2*Rcy(end-3*xL+1:end-2*xL);

% m(i, j+1)
yry3 = circshift(Rcy, -xL) - circshift(Rcy, -2*xL);
yry3(end-xL+1:end) = []; % remove extra values
yry3(1:xL) = Rcy(1:xL) + Rcy(xL+1:2*xL) - 2*Rcy(2*xL+1:3*xL);% bc at j = 1
% bc at j = end

% m(i, j+2)
yry4 = -circshift(Rcy, -2*xL);
yry4(end-xL+1:end) = []; % remove extra values - spdiags removes other boundary points
yry4(end-2*xL+1:end) = 0;
yry4(1:xL) = -2*Rcy(2*xL+1:3*xL);


ctc_y = spdiags([yry0, yry1, yry2, yry3, yry4],[2*xL, xL, 0, -xL, -2*xL], xL*yL, xL*yL)'./(4*dy.^2);




end