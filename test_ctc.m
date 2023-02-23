   clearvars
   close all
    load('inv_cntr test workspace2')
%     mesh.res_param1 = repmat([1:length(unique(param_x))]',length(unique(param_y)),1);
%     tmp_res = log10(mesh.res_param1);
%     mesh.px = [1:length(tmp_res)]';
%     mesh.py = px;
gamma_C = 1e-10;
input.diff_type = 3;
bc_flag = 0; % 0 - periodic for m,p ; 1 - 1st order neumann for m, p ; 2 - 2nd order neaumann m, 1st order neuman for p

tmp_lag = 0.05;


xL = length(unique(mesh.param_x));
yL = length(unique(mesh.param_y));

    if input.diff_type == 3
        
        % for now need to expand px, py as now defined on edges
        px = zeros((xL+1)*yL);
        py = zeros(xL*(yL+1));
        
        [ctc, ctx, cty] = ctc_half(mesh, input, gamma_C, itr, bc_flag);
        mesh.ctc = (ctc_x + ctc_y);
        Fm = ctc*log10(tmp_res) - ctx*px - cty*py;
    else
        if itr == 1 || input.inv_flag == 2
            mesh.Rc = eye(length(tmp_res));
        elseif input.inv_flag == -1 || input.inv_flag == -3
            mesh.Rc = diag(1./sqrt( (mesh.cx*log10(tmp_res) - px).^2 + (mesh.cy*log10(tmp_res) - py).^2 + gamma_C.^2));
        end
        ctc_x = mesh.cx'*mesh.Rc*mesh.cx; 
        ctc_y = mesh.cy'*mesh.Rc*mesh.cy;
        mesh.ctc = mesh.cx'*mesh.Rc*mesh.cx + mesh.cy'*mesh.Rc*mesh.cy;
        % Calculate intermediate m terms for removal of ghost points
        gradx = mesh.cx'*mesh.Rc*mesh.cx*log10(tmp_res);
        grady = mesh.cy'*mesh.Rc*mesh.cy*log10(tmp_res);
        cxRc = mesh.cx'*mesh.Rc;
        cyRc = mesh.cy'*mesh.Rc;   
        Fm = ((gradx - cxRc*px) +  (grady - cyRc*py));
    end

    ctc2 = mesh.ctc(domain_ind,:); 
    
    
    dx1 = (JTJ + tmp_lag*ctc2(:, domain_ind));
    dm1 = fem.array_jacobian.'*input.Wd'*mesh.Rd*input.Wd*(log10(input.real_data) - log10(fem.array_model_data)) - tmp_lag*Fm(domain_ind);
    
    dx1 = dx1\dm1;
    
    tmp_res = 10.^(log10(tmp_res(domain_ind)) + dx1); % update tmp_res, eliminate ghost nodes
    
    

    figure
    imagesc(reshape(log10(tmp_res), length(unique(mesh.param_x)),length(unique(mesh.param_y)))')
    colorbar
    title('m_{n+1}')
    figure
    imagesc(reshape(ctc_x*log10(mesh.res_param1), length(unique(mesh.param_x)),length(unique(mesh.param_y)))')
    colorbar
    title('ctcx*m')
    figure
    imagesc(reshape(ctc_y*log10(mesh.res_param1), length(unique(mesh.param_x)),length(unique(mesh.param_y)))')
    colorbar
    title('ctcy*m')
    
    figure
    plot(dx1)
    title('dx')

