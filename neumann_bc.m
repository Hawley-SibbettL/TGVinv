    function [num_param, domain_ind, param_x, param_y, u] = neumann_bc(num_param, param_xIn, param_yIn, uIn, neumann)
% This function adds the requred ghost points into the parameter vector for
% the neuman bc in the inverse problem
% neumann = [neumann constant, normal to boundary, order (1 or 2)]
% u is the parameter field (p or m). For p, u = [px, py]; each component is processed in parallel 

% used to weaken bc
dw_bc_flag = 0;
dw_bc_factor = 10;

neumann_const = neumann(1); % constant value of dm/dx, dm/dy at boundary. 
neumann_flag = neumann(2); %    [1] - constant gradient (for p).    [2] - constant 2nd order log(gradient) (for m).

tmp_x = unique(param_xIn);
tmp_y = unique(param_yIn);

% Get indices of edge and edge-adjacent cells
x1i = find(param_xIn == tmp_x(1));
x2i = find(param_xIn == tmp_x(2)); %x2i = x1i + 1;
xni = find(param_xIn == tmp_x(end));
xnmi = find(param_xIn == tmp_x(end-1)); %xnmi = xni - 1;
y1i = find(param_yIn == tmp_y(1));
y2i = find(param_yIn == tmp_y(2)); %y2i = y1i + length(tmp_y);
yni = find(param_yIn == tmp_y(end));
ynmi = find(param_yIn == tmp_y(end-1)); %ynmi = yni - length(tmp_y);
% corner points
% corners of parameter grid: cc = {c,c} is the corner, c1c = {c-1,c} etc 
cc = [find((param_xIn == tmp_x(1)).*(param_yIn == tmp_y(1))); find((param_xIn == tmp_x(1)).*(param_yIn == tmp_y(end))); find((param_xIn == tmp_x(end)).*(param_yIn == tmp_y(1))); find((param_xIn == tmp_x(end)).*(param_yIn == tmp_y(end)))];
c1c = [find((param_xIn == tmp_x(2)).*param_yIn == tmp_y(1)); find((param_xIn == tmp_x(2)).*param_yIn == tmp_y(end)); find((param_xIn == tmp_x(end-1)).*(param_yIn == tmp_y(1))); find((param_xIn == tmp_x(end-1)).*(param_yIn == tmp_y(end)))];
cc1 = [find((param_xIn == tmp_x(1)).*param_yIn == tmp_y(2)); find((param_xIn == tmp_x(1)).*param_yIn == tmp_y(end-1)); find((param_xIn == tmp_x(end)).*(param_yIn == tmp_y(2))); find((param_xIn == tmp_x(end)).*(param_yIn == tmp_y(end-1)))];
c1c1 = [find((param_xIn == tmp_x(2)).*(param_yIn == tmp_y(2))); find((param_xIn == tmp_x(2)).*(param_yIn == tmp_y(end-1))); find((param_xIn == tmp_x(end-1)).*(param_yIn == tmp_y(2))); find((param_xIn == tmp_x(end-1)).*(param_yIn == tmp_y(end-1)))];


% Using neumann bc, we need to add ghost points before and after the
% inversion domain (using one sided differences means derivatives are
% defined on the cell edges). Ordering will be:
% [ mesh.res_param1; m(x=0); m(y=0); m(x=n+1);  m(y=m+1) ]


% spatial edge  gradients
dx0 = tmp_x(2) - tmp_x(1);
dxn = tmp_x(end) - tmp_x(end - 1);
dy0 = tmp_y(2) - tmp_y(1);
dyn = tmp_y(end) - tmp_y(end - 1);

% downweight the p ghost point influence by increasing the distance
% ###########################################################
if dw_bc_flag == 1 %&& size(uIn, 2) > 1
    dx0 = dw_bc_factor*dx0;
    dxn = dw_bc_factor*dxn;
    dy0 = dw_bc_factor*dy0;
    dyn = dw_bc_factor*dyn;
end
% #############################################################j



dc = repmat([(dx0 + dy0); (dx0 + dyn); (dxn + dy0); (dxn + dyn)],1,size(uIn,2));

tmp_x = [tmp_x(1) - dx0; tmp_x; tmp_x(end) + dxn];
tmp_y = [tmp_y(1) - dy0; tmp_y; tmp_y(end) + dyn];

xL = length(tmp_x);
yL = length(tmp_y);

% reorder parameter vectors
param_x = repmat(tmp_x,length(tmp_y),1);
param_y = zeros(size(param_x));
for i = 1:yL
    param_y( (1 + (i-1)*xL):(i*xL) ) = repmat(tmp_y(i),length(xL),1);
end
domain_ind = find((param_x ~= tmp_x(1)).*(param_x ~= tmp_x(end)).*(param_y ~= tmp_y(1)).*(param_y ~= tmp_y(end)));


% Number of parameters including ghost points
num_param = length(param_x);


if neumann_flag == 1
    % for neumann conditions
%     ux0 = 10.^(log10(uIn(x1i,:)) + neumann_const*dx0);
%     uxn = 10.^(log10(uIn(xni,:)) + neumann_const*dxn);
%     uy0 = 10.^(log10(uIn(y1i,:)) + neumann_const*dy0);
%     uyn = 10.^(log10(uIn(yni,:)) + neumann_const*dyn);
%     uC = 10.^(log10(uIn(cc,:)) + neumann_const*dc);    % corners 
    ux0 = uIn(x1i,:) + neumann_const*dx0;
    uxn = uIn(xni,:) + neumann_const*dxn;
    uy0 = uIn(y1i,:) + neumann_const*dy0;
    uyn = uIn(yni,:) + neumann_const*dyn;
    uC = uIn(cc,:) + neumann_const*dc;    % corners 
elseif neumann_flag == 2 % assumes gradient of u calculated in log space
    % d2m/dx2 = 0.
    % u_{n+1} = u_n(1 + dx_{n+1,n}/dx_{n,n-1}) - u_{n-1}(dx_{n+1,n}/dx_{n,n-1}) +
    % neuman_const*dx_{n+1,n}. But we set dx equal.
    ux0 = 10.^(2*log10(uIn(x1i,:)) - log10(uIn(x2i,:)) + neumann_const*dx0);
    uxn = 10.^(2*log10(uIn(xni,:))- log10(uIn(xnmi,:)) + neumann_const*dxn);
    uy0 = 10.^(2*log10(uIn(y1i,:)) - log10(uIn(y2i,:)) + neumann_const*dy0);
    uyn = 10.^(2*log10(uIn(yni,:)) - log10(uIn(ynmi,:)) + neumann_const*dyn);
    uC = 10.^(4*log10(uIn(cc,:)) - 2*log10(uIn(c1c,:)) - 2*log10(uIn(cc1,:)) + log10(uIn(c1c1,:)) + neumann_const.*dc); % corners
end


ux0 = [uC(1,:); ux0; uC(2,:)];
uxn = [uC(3,:); uxn; uC(4,:)];
uy0 = [uC(1,:); uy0; uC(3,:)];
uyn = [uC(2,:); uyn; uC(4,:)];


% Need to replace corner points in ux's with the true corner
% values.
u = zeros(length(param_x), size(uIn, 2));

u(domain_ind, :) = uIn;

u(param_y == tmp_y(1),:) = uy0;
u(param_y == tmp_y(end),:) = uyn;
u(param_x == tmp_x(1),:) = ux0;
u(param_x == tmp_x(end),:) = uxn;
% needs fix - add in a second layer of boundary points

end