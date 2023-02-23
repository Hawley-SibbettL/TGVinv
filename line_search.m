% line_search(input, mesh, fem)
% input.lagrn_ls is vector of vlaues to include in the line search
% input.lagrn_w
% input.line search_flag == 1|2|3 do line searches in lambda|dm|both
% respectively
%
function [mesh_out, input, dx1] = line_search(input, mesh, fem, itr)

% Use input flag to set up type of inversion
if input.line_search_flag == 1 % lambda only, no line search along dm
    w = 1; % weights in line search direction = w.
    lagrn_ls = input.lagrn_ls;
    outer_loop_itr = 8;
elseif input.line_search_flag == 2
    w = input.lagrn_w;
    lagrn_ls = input.lagrn;
    outer_loop_itr = 8;
elseif input.line_search_flag == 3
    w = input.lagrn_w;
    lagrn_ls = input.lagrn_ls;
    outer_loop_itr = 8;
elseif input.line_search_flag == 4 % lambda only, golden section search
    w = 1; % weights in line search direction = w.
    lagrn_ls = input.lagrn_ls;
    outer_loop_itr = 8;
    g_ratio = (1 + sqrt(5))/2;
end

if input.line_search_flag == 4 % golden search
    if input.ls_limit_flag == 0
        % outer points defined
        if itr < 4
            lag1 = 1e-5;
            lag4 = 5e-1;
        else
            lag1 = 1e-6;
            lag4 = 5e-2;
        end
    else
        lag1 = 5e-5;
        lag4 = 5e-1;
    end
    % inner points caluculated via golden ratio (in log space)
    lag2 = 10.^( log10(lag4) - ( log10(lag4) - log10(lag1) ) / g_ratio );
    lag3 = 10.^( log10(lag1) + ( log10(lag4) - log10(lag1) ) / g_ratio );
    
    lagrn_ls = [lag1, lag2, lag3, lag4];
    
    p1 = 1; p2 = 2; p3 = 3; p4 = 4; % index of search points 1-4    
    
else % bisection: may have problems using 3 points to define interval
    % for bisection, replace lagrn ls with intial range:
    if itr < 4
        lagrn_ls = [1e-6, 5e-4, 2e-1];
    else
        lagrn_ls = [1e-7, 5e-5, 5e-2]; % shrink range for later iterations
    end
end

% initialise ls data structure
for i = 1:length(lagrn_ls)
    ls(i).lagrn = lagrn_ls(i);
    ls(i).w = w;
    ls(i).rms = zeros(size(w));
end


% precalculate terms in the inverse problem that are invarient with lambda
[num_param, domain_ind, param_x, param_y, tmp_res, dm1a, Fm, JTJ, ctc2, mesh] = inverse_precalc(input, mesh, fem, itr);
% initialises temporary mesh file
mesh_tmp = mesh;


% Line search is applied: RMS calculated for result of each
for q = 1:outer_loop_itr
    
    if input.line_search_flag == 4 % golden search
        %         lag3 = 10.^( log10(lag4) - ( log10(lag4) - log10(lag1) ) / g_ratio );
        %         lag2 = 10.^( log10(lag1) + ( log10(lag4) - log10(lag1) ) / g_ratio );
        if q == 1
            p = 1:length(lagrn_ls); % all valued
        else
            p = length(lagrn_ls); % new value of lambda to test
            % create ls listing for it
            ls(p).lagrn = lagrn_ls(p);
            ls(p).w = w;
            ls(p).rms = zeros(size(w));
        end
        
        
    else
        if q == 1
            p = 1:length(lagrn_ls);

        else
            p = q + 2;
            % identify new line search bisection point
            lagrn_ls(p) = 10.^mean([log10(lagrn_ls(qmin)), log10(lagrn_ls(qmax))]); % new q value to try
            ls(p).lagrn = lagrn_ls(p);
            ls(p).w = w;
            ls(p).rms = zeros(size(w));
        end
    end
    
    
    for i = p
        
        tmp_lag = lagrn_ls(i);
        fprintf(['line search, lagrn = ', num2str(tmp_lag), ', w = 1 *** '])
        
        % Here solve inverse problem - add result to tmp_mesh ready for forward
        % calculation
        dm1 = dm1a - tmp_lag*Fm(domain_ind);
        dx1 = (JTJ + tmp_lag*ctc2(:, domain_ind));
        %         dx1 = JTJ + tmp_lag*mesh.cx'*mesh.Rc*mesh.cx + mesh.cy'*mesh.Rc*mesh.cy;
        dx1 = dx1\dm1;
        
        % apply forward model to solution for w = 1
        mesh_tmp.res_param1 = 10.^(log10(tmp_res(domain_ind)) + dx1); % update tmp_res, eliminate ghost nodes
        
        % make sure the inverted model is reasonable(-ve res is a sign of an
        % ill-conditioned matrix)
        if any(mesh_tmp.res_param1 < 0)
            ls(i).rms(w == 1) = 10000; % attempt to catch unstable results
            if input.line_search_obj_flag == 1
                ls(i).phi_m(w==1) = [];
                if input.inv_flag == -3
                    ls(i).phi_tgv(w==1) = [];
                end
            end
        else
            h_warn = warning('error', 'MATLAB:nearlySingularMatrix');
            try
                [fem_tmp ,mesh_tmp ,input_tmp] = prop_change(itr+1,input,mesh_tmp,fem);
                fem_tmp = mes_control_fast(itr+1, input_tmp, mesh_tmp, fem, 1);
                % calculate/extract rms, print, then store ready for next loop.
                [~, fem_tmp, ~] = rms(itr + 1, 0, input_tmp, mesh_tmp, fem_tmp);
                ls(i).rms(w == 1) = fem_tmp.rms_sum1;
                if input.line_search_obj_flag == 1
                    ls = obj_calc(ls, input, fem_tmp, mesh_tmp, find(w==1), i, w);
                end
            catch
                ls(i).rms(w == 1) = 10000; % attempt to catch unstable results
                if input.line_search_obj_flag == 1
                    ls(i).phi_m(w==1) = [];
                    if input.inv_flag == -3
                        ls(i).phi_tgv(w==1) = [];
                    end
                end
            end
            warning(h_warn);
        end
        
        
        %     fprintf([' *** solution rms => ', num2str(fem_tmp.rms_sum1 ,'%.3g')])
        
        % here need to calcualte incremental objective functions - calling tgv
        % if necessary (probably make optional as good for paper but not
        % routine).
        
        
        % could add condition here to add w values if the solution is
        % diverging.
        
        
        % Here can apply line search in the direction of lambda if desired (
        % can also be applied at the end to save computational cost). It is
        % an especially good idea to do this here if rms is diverging. The
        % problem with this is that it is multiplicative with n_lagrn in
        % computation time
        if input.line_search_flag == 2 || input.line_search_flag == 3
            for j = 1:length(ls.w(i))
                if ls(i).w(j) ~= 1
                    fprintf(['line search along dm, lagrn = ', num2str(tmp_lag), ', w = ', num2str(ls(i).w(j)),' *** '])
                    % calculate new resistivity along search direction
                    mesh_tmp.res_param1 = 10.^(log10(tmp_res(domain_ind)) + dx1*ls(i).w(j));
                    
                    %                 fprintf([' *** solution rms => ', num2str(fem_tmp.rms_sum1 ,'%.3g')],'\n')
                    % store result
                    if any(mesh_tmp.res_param1 < 0) % catch unstable results
                        ls(i).rms(j) = 10000;
                        if input.line_search_obj_flag == 1
                            ls(i).phi_m(j) = [];
                            if input.inv_flag == -3
                                ls(i).phi_tgv(j) = [];
                            end
                        end
                    else
                        h_warn = warning('error', 'MATLAB:nearlySingularMatrix');
                        try
                            % apply forward problem and calculate rms
                            [fem_tmp ,mesh_tmp ,input_tmp] = prop_change(itr+1,input,mesh_tmp,fem);
                            fem_tmp = mes_control_fast(itr+1, input_tmp, mesh_tmp, fem, 1);
                            [~, fem_tmp, ~] = rms(itr + 1, 1, input_tmp, mesh_tmp, fem_tmp);
                            ls(i).rms(j) = fem_tmp.rms_sum1;
                            if input.line_search_obj_flag == 1
                                ls = obj_calc(ls, input, fem_tmp, mesh, j, i, w);
                            end
                        catch
                            ls(i).rms(j) = 10000;
                            if input.line_search_obj_flag == 1
                                ls(i).phi_m(j) = [];
                                if input.inv_flag == -3
                                    ls(i).phi_tgv(j) = [];
                                end
                            end
                        end
                        warning(h_warn);
                    end
                    
                end
            end
        end
        
    end
    
    
    % save lagrn_ls and rms results
    if input.line_search_flag == 4 % golden search
        
        if input.target_decrease > 0 
            
            target = input.target_decrease*fem.rms_sum1;
            % test where solution lies
            quad12 = ( ( ls(p1).rms - target) > 0 ) ~= ( ( ls(p2).rms - target) > 0 ); 
            quad23 = ( ( ls(p2).rms - target) > 0 ) ~= ( ( ls(p3).rms - target) > 0 ); 
            quad34 = ( ( ls(p3).rms - target) > 0 ) ~= ( ( ls(p4).rms - target) > 0 );
            
            if quad34 == 1 % if a solution lies in upper quadrant, search there
                lower_quad = 0;
            elseif quad23 == 1% solution is in middle quadrant - pick nearest side to solution
                lower_quad = 2*(abs( (ls(p2).rms - target) > 0 ) < abs( (ls(p3).rms - target) > 0 ));
            elseif quad12 == 1 % solution only in lower quadrant
                lower_quad = 1;
            else 
                disp('WARNING, DESIRED TARGET DECREASE OUTSIDE SEARCH BRACKETS') 
                % if target decrese smaller than best result, try to
                % minimise function. Otherwise go with upper quadrant
                if target < min( [ls(p1).rms, ls(p2).rms, ls(p3).rms, ls(p4).rms])
                    lower_quad = ls(p2).rms < ls(p3).rms;
                else
                    lower_quad = 0;    
                end
            end
            
            disp(['target decrease = ', num2str(target), ' lower_quad = ', num2str(lower_quad)])

        else
            lower_quad = ls(p2).rms < ls(p3).rms; % find which quadrant the minimum lies in
        end
        
        
        % find bracket solution lies within 
        if lower_quad  % solution in lower quadrant
            p4 = p3;
            p3 = p2;
            p2 = length(lagrn_ls) + 1;
            lagrn_ls( p2 ) = 10.^( log10(ls(p4).lagrn) - ( log10( ls(p4).lagrn ) - log10( ls(p1).lagrn ) ) / g_ratio  );
        else % solution in upper quadrant
            p1 = p2;
            p2 = p3;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
            p3 = length(lagrn_ls) + 1;
            lagrn_ls( p3 ) = 10.^( log10(ls(p1).lagrn) + ( log10( ls(p4).lagrn ) - log10( ls(p1).lagrn ) ) / g_ratio  );
        end
        disp(['lagrn1 = ', num2str(lagrn_ls(p1)),'; lagrn2 = ', num2str(lagrn_ls(p2)),'; lagrn3 = ', num2str(lagrn_ls(p3)),'; lagrn4 = ', num2str(lagrn_ls(p4))])
        
    else
        if input.target_decrease > 0
            if q == 1
                target_pc = (([ls(1).rms, ls(2).rms, ls(3).rms] - input.target_decrease*fem.rms_sum1));
                upper_test = target_pc(3) > 0 && target_pc(2) < 0;
                lower_test = target_pc(1) > 0 && target_pc(2) < 0;
                
                if itr == 1
                    if upper_test && lower_test % target lies in both brackets
                        % choose higher lambda bracket to select on first itr
                        qmax = 3;
                        qmin = 2;
                    elseif lower_test
                        qmin = 1;
                        qmax = 2;
                    else % pick upper bracket - if it lies outside range, then bisection won't help
                        qmax = 3;
                        qmin = 2;
                    end
                else
                    [target_pc, qsort] = sort(abs([ls(1).rms, ls(2).rms, ls(3).rms] - input.target_decrease*fem.rms_sum1));
                    qmin = min([(qsort(1)), (qsort(2))]);
                    qmax = max([(qsort(1)), (qsort(2))]);
                end
            else
                [target_pc, qsort] = sort(abs([ls(qmin).rms, ls(i).rms, ls(qmax).rms] - input.target_decrease*fem.rms_sum1));
                ql = [qmin, i, qmax];
                if ql(qsort(3)) == i
                    qmin = min([ql(qsort(1)), i]);
                    qmax = max([ql(qsort(1)), i]);
                else
                    qmin = min([ql(qsort(1)), ql(qsort(2))]);
                    qmax = max([ql(qsort(1)), ql(qsort(2))]);
                end
            end
            
        else
            if q == 1
                % if have a target decrease, check to find if it falls in which
                %bracket - if in doubt choose higher one.
                % assume target decrease convexity right from now
                [~, qsort] = sort(abs([ls(1).rms, ls(2).rms, ls(3).rms] - input.target_decrease*fem.rms_sum1));
                qmin = qsort(1);
                qmax = qsort(2);
            else
                [~, qsort] = sort(abs([ls(qmin).rms, ls(i).rms] - input.target_decrease*fem.rms_sum1));
                if qsort(1) == 2
                    qmax = qmin;
                    qmin = i;
                elseif qsort(1) == 1
                    qmax = i;
                end
            end
        end
    end
    
    if input.line_search_flag == 4
        % No convergence criteria for now
    else
        if 100*abs(ls(qmin).rms - ls(qmax).rms)/ls(qmin).rms < input.conv_rate
            break
        end
    end
    
end

% loop creates extra entries.
if input.line_search_flag == 4
    lagrn_ls(end) = [];
end

% next analyse result

% compare with target misfit, find closest result to target or just
% smallest (could just implement smallest rms first)
for i = 1:length(lagrn_ls)
    [min_rms, min_ind] = min(ls(i).rms);
    ls(i).rms_min = min_rms;
    ls(i).w_min = w(min_ind);
end

% sets a target % rms misfit for the next iteration.
[min_rms, min_ind] = min(abs([ls(:).rms_min] - input.target_decrease*fem.rms_sum1));
min_lagrn = ls(min_ind).lagrn;
min_w = ls(min_ind).w_min;

% Apply any bisection searches if desired to further refine the result

% can apply line search in direction of lambda here just for the final
% solution to be more efficient

% recalculate final m to save storing all the intermediate results
tmp_lag = min_lagrn;
dm1 = dm1a - tmp_lag*Fm(domain_ind);
dx1 = (JTJ + tmp_lag*ctc2(:, domain_ind));
dx1 = dx1\dm1;
mesh_out = mesh;
mesh_out.res_param1 = 10.^(log10(tmp_res(domain_ind)) + dx1*min_w); % update tmp_res, eliminate ghost nodes

mesh_out.ls = ls;
input.lagrn = min_lagrn;
fprintf('\n')
disp([' *** line search results: lagrn = ', num2str(min_lagrn ,'%.3g'), '; w = ', num2str(min_w ,'%.3g'), '; rms = ' num2str(min_rms ,'%.3g'),'  ***'])

end





function [num_param, domain_ind, param_x, param_y, tmp_res, dm1a, Fm, JTJ, ctc2, mesh] = inverse_precalc(input, mesh, fem, itr)

% cutoffs for |.| = 0
gamma_D = mesh.gamma_d;
gamma_C = mesh.gamma_c;

if input.ghost_flag == 1
    [num_param, domain_ind, param_x, param_y, tmp_res] = neumann_bc(mesh.num_param, mesh.param_x, mesh.param_y, mesh.res_param1, [0, input.dm_bc]);
else
    num_param = mesh.num_param;
    domain_ind = 1:num_param;
    param_x = mesh.param_x;
    param_y = mesh.param_y;
    tmp_res = mesh.res_param1;
end


if input.inv_flag ~= -3
    px = zeros(size(tmp_res));
    py = px;
elseif input.ghost_flag == 1
    [~, ~, ~, ~, p_tmp] = neumann_bc(mesh.num_param, mesh.param_x, mesh.param_y, [mesh.px, mesh.py], [0, 1]);
    px = p_tmp(:,1);
    py = p_tmp(:,2);
else
    px = mesh.px;
    py = mesh.py;
end

% Regulariser term weight - IRLS part
if (input.inv_flag == -1 && itr ~= 1) || (input.inv_flag == -3 && itr ~= 1)
    mesh.Rc = diag(1./sqrt( (mesh.cx*log10(tmp_res) - px).^2 + (mesh.cy*log10(tmp_res) - py).^2 + gamma_C.^2));
else
    mesh.Rc = eye(length(tmp_res));
end

if input.l1_data_flag == 0
    mesh.Rd = eye(length(input.real_data));
else
    mesh.Rd = diag(1./sqrt((input.Wd*(log10(input.real_data) - log10(fem.array_model_data))).^2 + gamma_D^2));
end

mesh.ctc = mesh.cx'*mesh.Rc*mesh.cx + mesh.cy'*mesh.Rc*mesh.cy;


JTJ = fem.array_jacobian.'*input.Wd'*mesh.Rd*input.Wd*fem.array_jacobian;

Fm = mesh.cx'*mesh.Rc*mesh.cx*(log10(tmp_res) - px) +  mesh.cy'*mesh.Rc*(mesh.cy*log10(tmp_res) - py);

dm1a = fem.array_jacobian.'*input.Wd'*mesh.Rd*input.Wd*(log10(input.real_data) - log10(fem.array_model_data));
ctc2 = mesh.ctc(domain_ind,:);

end

function ls = obj_calc(ls, input, fem, mesh, j, i, w)
% j is index of the value of w being used

% Calculate data term
if input.l1_data_flag == 0
    ls(i).phi_d(w == 1) = sum((input.Wd*(log10(input.real_data) - log10(fem.array_model_data))).^2);
else
    ls(i).phi_d(w == 1) = sum(abs(input.Wd*(log10(input.real_data) - log10(fem.array_model_data))));
end

% need to apply bc for regularisation term
[num_param, domain_ind, param_x, param_y, res_param] = neumann_bc(mesh.num_param, mesh.param_x, mesh.param_y, mesh.res_param1, [0, input.dm_bc]);


if input.inv_flag ~= -3
    px = zeros(size(res_param));
    py = px;
end

grad_mx = mesh.cx*log10(res_param);
grad_my = mesh.cy*log10(res_param);

% tgv p calculation if tgv inversion
if input.inv_flag == -3
    
    input.sp2_itr = 100;
    %initialise p - includes bc ghost point terms
    if input.p_init_flag == 0 || input.p_init_flag == 2
        px = zeros(mesh.num_param,1);
        py = mesh.px;
    elseif input.p_init_flag == 1 % p = grad(m) on domain of m.
        px = mesh.cx*res_param;
        px = px(domain_ind);
        py = mesh.cy*res_param;
        py = py(domain_ind);
    end
    
    gamma_p = mesh.gamma_p;               % cutoff for |x| -> 0
    gamma_c = mesh.gamma_c;
    % Used in performance metrics
    p1 = zeros(2*mesh.num_param,1);
    rms_p = zeros(1,input.sp2_itr+1);
    % for use in p loop
    sp2_itr = 100;
    
    
    for k = 2:(sp2_itr+1)
        % add ghost points to p.
        [~, ~, ~, ~, p_tmp] = neumann_bc(mesh.num_param, mesh.param_x, mesh.param_y, [px, py], [0, 1]);
        px_tmp = p_tmp(:,1);
        py_tmp = p_tmp(:,2);
        
        % calculate weights
        if input.p_init_flag == 2 && k == 2 % l2 1st iteration as initialisation
            Rm = eye(length(grad_mx));
            Rp = eye(length(mesh.cx'*px_tmp));
        else
            Rm = diag(1./sqrt((grad_mx - px_tmp).^2 + (grad_my - py_tmp).^2 + gamma_c^2));
            Rp = diag(1./sqrt((mesh.cx'*px_tmp).^2 + (mesh.cy'*py_tmp).^2 + 0.5*(mesh.cx'*py_tmp + mesh.cy'*px_tmp).^2 + gamma_p.^2));
        end
        
        % Set up matrix equation for px, py and solve
        a11 = Rm + input.tgv_lagrn*(mesh.cx*Rp*mesh.cx' + 0.5*mesh.cy*Rp*mesh.cy');
        a12 = input.tgv_lagrn*0.5*mesh.cy*Rp*mesh.cx';
        a21 = input.tgv_lagrn*0.5*mesh.cx*Rp*mesh.cy';
        a22 = Rm + input.tgv_lagrn*(mesh.cy*Rp*mesh.cy' + 0.5*mesh.cx*Rp*mesh.cx');
        b1 = Rm*mesh.cx*log10(res_param);
        b2 = Rm*mesh.cy*log10(res_param);
        
        A = [a11, a12; a21, a22];
        b = [b1; b2];
        p2 = A\b;
        
        clear A b b1 b2 a11 a12 a21 a22
        
        % separate px, py
        px = p2(1:length(p2)/2);
        py = p2((length(p2)/2+1):length(p2));
        % remove ghost points
        px = px(domain_ind);
        py = py(domain_ind);
        p3 = p2;
        p2 = [px; py];
        
        % calculate convergence and print (may supress in future)
        rms_p(k-1) = sqrt(mean((p2 - p1).^2))/(0.5*sqrt(mean( (mesh.cx*log10(res_param)).^2 + (mesh.cy*log10(res_param)).^2)));
        disp(['tgv iteration ',num2str(k-1),' : rms dp = ',num2str(rms_p(k-1))])
        if ((k > 4) && (rms_p(k-1) < input.p_conv)) %  stop if change smaller than p_conv%
            sp2_i = i;  % store number of iterations needed
            break
        elseif k > 2 && (rms_p(k-1) > rms_p(k-2)) || ~isreal(p2)|| any(p2<0) % stop if diverging
            %             sp2_i = i-1;  % store number of iterations needed
            px = p1(1:length(p1)/2);
            py = p1((length(p1)/2+1):length(p1));
            break
        end
        
        %         if i == (input.sp2_itr+1)
        %             sp2_i = i;
        %         end
        
        p1 = p2;
        
    end
    
    
    % TGV smooth objective term
    tgv_term2 = (mesh.cx'*p3(1:length(p3)/2)).^2 + (mesh.cy'*p3((length(p3)/2+1):length(p3))).^2 + 0.5*(mesh.cx'*p3((length(p3)/2+1):length(p3)) + mesh.cy'*p3(1:length(p3)/2)).^2;
    ls(i).phi_tgv(j) = sum(sqrt(tgv_term2(domain_ind)));
end


grad_x = grad_mx - p3((1:length(p3)/2));
grad_y = grad_my - p3((length(p3)/2+1):length(p3));

if input.inv_flag == -2 % saves l2 regularisation term
    ls(i).phi_m(j) = sum((grad_x(domain_ind)).^2 + (grad_y(domain_ind)).^2);
elseif input.inv_flag == -1 || input.inv_flag == -3 % saves 'l1' penalty
    ls(i).phi_m(j) = sum(sqrt((grad_x(domain_ind)).^2 + (grad_y(domain_ind)).^2));
end


end