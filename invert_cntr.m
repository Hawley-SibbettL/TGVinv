function [mesh,fem,input]=invert_cntr(itr,tmp_flag,ip_cnt,input,mesh,fem)



%    /* ******************************************************** */
%    /* *************    Error weighting	******************** */
%    /* ******************************************************** */
%
%      /* perform error weighting only when inv_flags are 3 or 4 */



% wd=abs((log10(input.real_data))-log10((fem.array_model_data)));
% %
% %
% %
% %
% % for i=1:input.num_mes
% %     if 100*abs(input.real_data(i)-fem.array_model_data(i))/input.real_data(i) <fem.rms_sum1
% %         wd(i)=fem.rms_sum1/100;
% %     end
% % end
% %
% input.Wd=diag(1./wd);

% Unknown flags - seem to represent error weighting but can't find their
% initialisation
if input.inv_flag==3 || input.inv_flag==4
    
    sum1=0;
    sum2=0;
    
    %/* find initially sums */
    for i=1:input.num_mes
        sum1=sum1+ (input.Wd(i,i).^0.5)   *abs(log10(input.real_data(i))-log10(fem.array_model_data(i)));
        sum2=sum2+ (input.Wd(i,i).^0.25) *(abs(log10(input.real_data(i))-log10(fem.array_model_data(i))).^0.5);
    end
    hit=0;
    w_trial=zeros(input.num_mes);
    for i=1:input.num_mes
        w_trial(i,i)=(input.Wd(i,i)^0.5) / abs(log10(input.real_data(i))-log10(fem.array_model_data(i)));
        w_trial(i,i)=w_trial(i,i)*(sum1/sum2);
        if w_trial(i,i)>input.Wd(i,i)
            w_trial(i,i)=input.Wd(i,i);
            hit=hit+1;
        end
    end
    
    % /* find L1 ratio */
    sum1=0; sum2=0;
    for i=1:input.num_mes
        sum1=sum1+input.Wd(i,i)*abs(log10(input.real_data(i))-log10(fem.array_model_data(i)));
        sum2=sum2+w_trial(i,i)*abs(log10(input.real_data(i))-log10(fem.array_model_data(i)));
    end
    l1=sum1/sum2;
    
    %/* accept changes only when L1 ratio >1 */
    if(l1>1)
        
        input.Wd=w_trial;
        
    end
    
end




% Save original lagrange multiplier for IP
if (ip_cnt==1 && itr==1) ;input.original_lagrn=input.lagrn; end
if (ip_cnt==2 && itr==1) ;input.lagrn=input.original_lagrn; end

% Reduce the lagrange multiplier every itr
%  if (itr~=1  && tmp_flag==1); input.lagrn=input.lagrn/input.lagrn_reduction; end %&& itr<5
%
% Quasi Newton update Jacobian
if(itr>=2 && input.jacobian_flag==1) ;fem=quasi_newton(input,mesh,fem); end


%
%  %				  /* find jacobian scale */
%  sum=0;
%  	for j=1:mesh.num_param
% 	  for i=1:mesh.num_param
% 		  tmp11=JTJ(i,j);
% 		  sum=sum+tmp11*tmp11;
%       end
%     end
%
%  jacscale=sqrt(sum);



%--------------------ACB+dx1----------------------------------------------
% dx1 = JTJ + labda*CTC - resistiviy in log units. LHS of inversion equ.
% Should work for l1 too, as JTJ, ctc are updated (as long as acb not
% ticked)
JTJ=fem.array_jacobian.'*input.Wd*fem.array_jacobian;


%   if input.acb_flag==0
tmp_lag=input.lagrn;
%       dx1=(JTJ + tmp_lag*mesh.ctc);
%       ctc=mesh.ctc;
%   elseif (input.acb_flag==1)
%       [ctc,L1]=acb(input,mesh,fem);
%       dx1=(JTJ+ctc);
%       % keep ACB
%       fem.L1=L1;
%   end



% Here calcualtes Gauss-Newton updates from my TGV modifications

% Modes: TV (-1), TGV (-3), l2-l2 (2), l1-l2 (-2)
% note that l2-l2 overwrites original GN functionality, so that Wd is used
% consistently across all my inversions

if input.inv_flag == -1 || input.inv_flag == -3 || input.inv_flag == -2 || input.inv_flag == 2

    
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
end


% line search happens here
if input.line_search_flag ~= 0
    [mesh, input, dx1] = line_search(input, mesh, fem, itr);
    tmp_res = mesh.res_param1;
else
    
    % data terms always the same
    if input.l1_data_flag == 0
        mesh.Rd = speye(length(input.real_data));
    else
        mesh.Rd = spdiags(1./sqrt( ( input.Wd*( log10(input.real_data) - log10(fem.array_model_data) ) ).^2 + mesh.gamma_d^2 ), 0, mesh.num_param, mesh.num_param);
    end
    
    JTJ=fem.array_jacobian.'*input.Wd'*mesh.Rd*input.Wd*fem.array_jacobian;
    
    
    % Regulariser term weight - IRLS part
    if input.diff_type == 3
        [ctc_x, ctc_y, cxtrc, cytrc] = ctc_calc3(mesh, input, px, py, gamma_C, itr);
%         [ctc_x, ctc_y, cxtrc, cytrc] = ctc_calc4(mesh, input, px, py, gamma_C, 1);
        mesh.ctc = ctc_x + ctc_y;
        Fm = (ctc_x*log10(tmp_res) - cxtrc*py) + (ctc_y*log10(tmp_res) - cytrc*py);
    else
        if input.inv_flag == 2 || ( itr == 1 && input.m_init_flag == 2 )  
            mesh.Rc = speye(length(tmp_res));
        elseif input.inv_flag == -1 
            mesh.Rc = spdiags(1./sqrt( (mesh.cx*log10(tmp_res)).^2 + (mesh.cy*log10(tmp_res)).^2 + mesh.gamma_c.^2), 0, mesh.num_param, mesh.num_param);
        elseif input.inv_flag == -3
            mesh.Rc = spdiags(1./sqrt( (mesh.cx*log10(tmp_res) - px).^2 + (mesh.cy*log10(tmp_res) - py).^2 + mesh.gamma_c.^2), 0, mesh.num_param, mesh.num_param);
            
            %             if input.diff_type == 2
            %                 mesh.Rc = diag(1./sqrt( ((mesh.cx*log10(tmp_res) - px) + (mesh.cy*log10(tmp_res) - py)).^2 + gamma_C.^2));
            %             end
        end
        mesh.ctc = mesh.cx'*mesh.Rc*mesh.cx + mesh.cy'*mesh.Rc*mesh.cy;
        % Calculate intermediate m terms for removal of ghost points
%         gradx = mesh.cx'*mesh.Rc*mesh.cx*log10(tmp_res);
%         grady = mesh.cy'*mesh.Rc*mesh.cy*log10(tmp_res);
%         cxRc = mesh.cx'*mesh.Rc;
%         cyRc = mesh.cy'*mesh.Rc;   
%         Fm = ((gradx - cxRc*px) +  (grady - cyRc*py));
        Fm = mesh.cx'*mesh.Rc*(mesh.cx*log10(tmp_res) - px) +  mesh.cy'*mesh.Rc*(mesh.cy*log10(tmp_res) - py);
    end

    ctc2 = mesh.ctc(domain_ind,:); 
    
    
    dx1 = (JTJ + tmp_lag*ctc2(:, domain_ind));
    dm1 = fem.array_jacobian.'*input.Wd'*mesh.Rd*input.Wd*(log10(input.real_data) - log10(fem.array_model_data)) - tmp_lag*Fm(domain_ind);
    
    dx1 = dx1\dm1;
    
    tmp_res = 10.^(log10(tmp_res(domain_ind)) + dx1); % update tmp_res, eliminate ghost nodes
    
end




% updates res_param

for i=1:mesh.num_param
    % GN; GN diff; L-M; Occam dif
    if input.inv_flag==6 || input.inv_flag==0 ||input.inv_flag==5
        % mesh_param2 holds parameters for last itr. If rms terminates
        % inversion, nothing will be stored.
        mesh.res_param1(i)=10^(log10(mesh.res_param2(i)) + dx1(i)); % Stored as resistivity
    elseif input.inv_flag == -1 || input.inv_flag == -3 || input.inv_flag == -2 || input.inv_flag == 2
        mesh.res_param1 = tmp_res;
    else % TL, Occam
        mesh.res_param1(i)=10^(dx1(i));
    end
    %     if imag(mesh.res_param1(i))>0 ; mesh.res_param1(i)=complex(real(mesh.res_param1(i)),-0.01); end
    if input.limit_res==1
        if mesh.res_param1(i)>input.max_res; mesh.res_param1(i)=input.max_res; end
        if mesh.res_param1(i)<input.min_res; mesh.res_param1(i)=input.min_res; end
    end
    
end



% keep resolution matrix
fem.resolution_matrix=dx1\(JTJ);



tt=sprintf('**  ITERATION =>  %d  **\n',itr);
tt=cellstr(tt);
drawnow;



end















