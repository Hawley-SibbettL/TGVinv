% Function to apply l_curve method to a gauss-newton type inversion (single
% time at present)

function [input,mesh] = l_curve(mesh,fem,input,itr)

old_lagrn = input.lagrn;
Nlc = 100;       % number of (log spaced) values in curve
lagmin = 0.0001; % Minimum lagr multiplier value to be used
lagmax = 10;     % Maximum lagrange multiplier to be used

log_lagr = linspace(log10(lagmin),log10(lagmax),Nlc);
lagr = 10.^log_lagr;
lc_misfit = zeros(1,Nlc);
lc_roughness = zeros(1,Nlc);

% Calculate terms in inversion which are constant with lagr first
if input.inv_flag == -1
    JTJ=fem.array_jacobian.'*input.Wd'*mesh.Rd*input.Wd*fem.array_jacobian;
else
    JTJ=fem.array_jacobian.'*input.Wd*fem.array_jacobian;
end
if input.inv_flag ~= -1 % l2 data constraint
    JtWd_misfit = fem.array_jacobian.'*input.Wd*(log10(input.real_data)-log10(fem.array_model_data));
else % l1 data constraint
    JtWd_misfit = fem.array_jacobian.'*input.Wd'*mesh.Rd*input.Wd*(log10(input.real_data)-log10(fem.array_model_data));
end
ctcm = mesh.ctc*log10(mesh.res_param1);



% For each value of lagrange multiplier, calculate l-curve points
for i = 1:length(lagr)
    
    % Carry out inversion
    [tmp_model dm] = inversion_calc(mesh, JTJ, JtWd_misfit, ctcm, lagr(i));
    
    % Calculate data misfit of l-curve point (l1 or l2 norm as appropriate)
    approx_newdata = log10(fem.array_model_data) + fem.array_jacobian*(dm);
    
    if input.inv_flag == -1     % l1
        lc_misfit(i) = sum(sum( abs(log10(input.real_data) - (approx_newdata)) ));
    else    % l2
        lc_misfit(i) = sum(sum( (log10(input.real_data) - log10(approx_newdata)).^2 ));
    end
    
    % Calculate model roughness of l-curve point (l1 or l norm as nappropriate)
    if input.inv_flag ~=2   % l1
        lc_roughness(i) = (sum(sum( abs((mesh.cx + mesh.cy)*log10(tmp_model)) )));
    else    % l2
        lc_roughness(i) = (sum(sum( ((mesh.cx + mesh.cy)*log10(tmp_model)).^2 )));
    end    
end

for i = 1:size(lc_roughness,1)
    lc_roughness(i,:) = smooth(lc_roughness(i,:),5)';
    lc_misfit = smooth(lc_misfit(i,:),5)';
end

% calculate l-curve curvature
dreg = gradient(log10(lc_roughness));
ddreg = gradient(dreg);
dmis = gradient(log10(lc_misfit));
ddmis = gradient(dmis);



lc_curvature = (dmis.*ddreg - ddmis.*dreg)./((dreg.^2 + dmis.^2).^1.5);



[~, lagr_ind] = max(lc_curvature(i,:));

% % the maximum of the curvature is the the lagr selected by the l-curve
% [~, peak_ind] = findpeaks(lc_curvature);
% % Select peak with smallest maximum spacing between adjacent points on the l-curve in both dim.
% % This is because the points become more spread in the either the x/y 
% % direction further away from the point of maximum curvature
% % Normalised by range value of each axis
% if length(peak_ind) ~=1;  peak_ind(peak_ind == 1) = []; end; % removes ends of curve as not likely to be good (or curve doesn't cover right space)
% if length(peak_ind) ~=1;  peak_ind(peak_ind == Nlc) = []; end;
% [~, max_ind] = min(max(abs(lc_roughness(peak_ind) - lc_roughness(peak_ind-1))/range(lc_roughness),abs(lc_misfit(peak_ind) - lc_misfit(peak_ind-1))/range(lc_misfit)));
% lagr_ind = peak_ind(max_ind);

lc_lagr = lagr(lagr_ind);

% Save results of calculation
% mesh.lc_peak_ind = peak_ind;
mesh.lc_curvature = lc_curvature;
mesh.lc_lagr = lc_lagr;
mesh.lc_roughness = lc_roughness;
mesh.lc_misfit = lc_misfit;
% mesh.lc_curvature_smooth = lc_curvature_smooth;
input.lagrn = lc_lagr;
mesh.lagr_ind = lagr_ind;
if itr == 1
    mesh.lc_trial_lagr = lagr;
end

disp(['l-curve lagr = ', num2str(input.lagrn)]);
% reduction_factor = 1/0.75;
% % Adjust for cooling behaviour
% if input.lagrn < old_lagrn/4;
%     input.lagrn = old_lagrn/reduction_factor;
%     disp(['This is less than old lagrn/',num2str(reduction_factor),' - using old lagrn/',num2str(reduction_factor),' instead (= ',num2str(input.lagrn),')'])
% % elseif input.lagrn > 2*old_lagrn
% %     input.lagrn = old_lagrn/reduction_factor;
% %     disp(['This is greater than old lagrn - using old lagrn/',num2str(reduction_factor) ,'instead (= ',num2str(input.lagrn),')'])
% end
input.lagrn = max(old_lagrn/2,input.lagrn);
if input.lagrn == old_lagrn/2
    disp(['Smaller than old_lagrn/2, replacing with old_lagrn/2 = ', num2str(input.lagrn)])
end

mesh.lc_lagr = input.lagrn;


end



    function [tmp_res_param, dx1] = inversion_calc(mesh, JTJ, JtWd_misfit, ctcm, tmp_lag)
        
        dx1=(JTJ + tmp_lag*mesh.ctc);
        
        tmpmtx=tmp_lag*ctcm;
        
        dm1 = JtWd_misfit - tmpmtx;
        
        dx1=dx1\dm1; % matrix inverse of dx1*dm1 solves minimisation equation
        
        tmp_res_param = 10.^(log10(mesh.res_param2) + dx1); % New model
               
    end
