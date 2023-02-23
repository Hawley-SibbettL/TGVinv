% Performs half point-difference calculation
% p2 = [px; py]
% m - log10(resistivty), mu - tgv_lagrn
function [p2] = half_diff_p_calc(mesh, lagrn ,tgv_lagrn, gamma_c, gamma_p)

px = mesh.px;
py = mesh.py;







% Rm calculated at cell centres
Rm = 1./sqrt((mve2 - mve1 - px).^2 + (mhe2 - mhe1 - py).^2 + gamma_c.^2);
Rm = sparse(1:xL*yL, 1:xL*yL, Rm, xL*yL, xL*yL);

% Rm*Cx*m, Rm*Cy*m  
[cx, cy] = c_calc(tmp_x, tmp_y); % Can use fd matrices as points are all at cell centres
b1 = Rm*cy*m;
b2 = Rm*cy*m;

% Rpx
Rpx1 = sqrt(().^2 + ().^2 + + ().^2 + gamma.^2);
% Rpx1 = sqrt((m - mx_ghost1 - px_vedge1).^2 + (mc_b1 - mc_t1 - py_vedge1).^2 + gamma.^2);


a11 = Rm + input.tgv_lagrn*(mesh.cx*Rp*mesh.cx' + 0.5*mesh.cy*Rp*mesh.cy');
a12 = input.tgv_lagrn*0.5*mesh.cy*Rp*mesh.cx';
a21 = input.tgv_lagrn*0.5*mesh.cx*Rp*mesh.cy';
a22 = Rm + input.tgv_lagrn*(mesh.cy*Rp*mesh.cy' + 0.5*mesh.cx*Rp*mesh.cx');
b1 = Rm*mesh.cx*log10(res_param);
b2 = Rm*mesh.cy*log10(res_param);

A = [a11, a12; a21, a22];
b = [b1; b2];

p2 = A\b;
end