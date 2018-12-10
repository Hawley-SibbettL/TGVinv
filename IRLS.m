function mesh = IRLS(input,fem,mesh)
% Calcultes the weighting functions for 2D single-slice IRLS Gauss Newton
% (see Farquarson and Oldenburg 1996)
%
% mesh = IRLS(input,fem,mesh)


gamma = 1e-4;               % cutoff for |x| = 0

if input.inv_flag ~= -3
    mesh.px = zeros(size(mesh.res_param1));
    mesh.py = mesh.px;
end

% Data term weight
mesh.Rd = diag(1./sqrt(input.Wd*(log10(input.real_data) - log10(fem.array_model_data)).^2 + gamma^2));
% mesh.Rd = mesh.Rd + (mesh.Rd>gamma)*(1./gamma);

% Regulariser term weight
mesh.Rc = diag(1./sqrt( (mesh.cx*log10(mesh.res_param1) - mesh.px).^2 + (mesh.cy*log10(mesh.res_param1) - mesh.py).^2 + gamma.^2));
% mesh.Rc = mesh.Rc + (mesh.Rc>gamma)*(1./gamma);

% mesh.ctc is now the weighted version- Rc only appears in weights. Rd
% needs to be used separately later
mesh.ctc = mesh.cx'*mesh.Rc*(mesh.cx) + mesh.cy'*mesh.Rc*(mesh.cy);




% Rd for discrepency term
% discrepency = (log10(input.real_data) - log10(fem.array_model_data));
% Rd = zeros(input.num_mes);

% for i = 1: input.num_mes
%     inv_Rdi = abs(input.Wd(i,:)*discrepency);
%     
%      % avoids dividing by 0
%     if inv_Rdi < gamma
%                 inv_Rdi =  gamma;
%     end
%     
%     Rd(i,i) = 1./inv_Rdi;
% end


% % Rc for smoothing term
% Rcx = zeros(mesh.num_param);
% Rcy = zeros(mesh.num_param);
% 
% for i = 1:mesh.num_param
%     inv_Rcix = abs(mesh.cx(i,:)*log10(mesh.res_param1) - mesh.px(i));
%     inv_Rciy = abs(mesh.cy(i,:)*log10(mesh.res_param1) - mesh.py(i));
% 
%     
% %     avoids dividing by 0
%     if inv_Rcix < gamma
%                 inv_Rcix =  gamma;
%     end
%     if inv_Rciy < gamma
%                 inv_Rciy =  gamma;
%     end
% 
%     
%     Rcx(i,i) = 1./inv_Rcix;
%     Rcy(i,i) = 1./inv_Rciy;
% end        
%     

% mesh.Rcx = Rcx;
% mesh.Rcy = Rcy;
% mesh.ctc = mesh.cx'*Rcx*mesh.cx + mesh.cy'*Rcy*mesh.cy;



end