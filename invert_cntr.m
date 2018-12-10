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
 if (itr~=1  && tmp_flag==1); input.lagrn=input.lagrn/input.lagrn_reduction; end %&& itr<5
 
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


  if input.acb_flag==0
      tmp_lag=input.lagrn;
      dx1=(JTJ + tmp_lag*mesh.ctc);
      ctc=mesh.ctc;
  elseif (input.acb_flag==1)
      [ctc,L1]=acb(input,mesh,fem);
      dx1=(JTJ+ctc);
      % keep ACB
      fem.L1=L1;
  end

 
  


 %---------------C*dx------------------------------------------------------ 
 % tmptx stores part of RHS of inversion equation for GN -lambda*F*m 
 % tmptx2 only for difference/occam
 
 % prevents crashing
      tmpmtx2=zeros(input.num_mes,1);
           tmpmtx=zeros(mesh.num_param,1);


 
 
 % Gauss Newton ) 
 if input.inv_flag == 2 
     if input.acb_flag==0
        tmpmtx=tmp_lag*mesh.ctc*log10(mesh.res_param1);
     else
        tmpmtx=ctc*log10(mesh.res_param1);  % Why no lagrange multiplier? = CTCm (end term in inversion) must be in acb workings
     end
     tmpmtx2=zeros(input.num_mes,1);
 end
 
% Not needed due to coupled x y terms 
%  % TGV s.p. 1 - as l1 GN but separate ctc terms for x and y gradients
%  if input.inv_flag == -3
%      tmpmtx = tmp_lag*( mesh.cx'*mesh.Rcx*(mesh.cx*log10(mesh.res_param1) - mesh.px) + mesh.cy'*mesh.Rcy*(mesh.cy*log10(mesh.res_param1) - mesh.py) );
%      tmpmtx2=zeros(input.num_mes,1);
%  end
 
 
 % Levenburg-Marquadt or Occam difference - both equal 0. Neither have
 % smoothing term
  if input.inv_flag==0 || input.inv_flag==5 
     tmpmtx=zeros(mesh.num_param,1);
     tmpmtx2=zeros(input.num_mes,1);
  end
 
 % GN difference - background used
 if input.inv_flag==6 
     
     if input.acb_flag==0
         if itr>1 
             tmpmtx=tmp_lag*ctc*(log10(mesh.res_param1) - log10(mesh.bgr_param));
         else
              tmpmtx=tmp_lag*ctc*(log10(mesh.res_param1));
         end
     else
        if itr>1
            tmpmtx=ctc*(log10(mesh.res_param1) - log10(mesh.bgr_param));
        else
            tmpmtx=ctc*(log10(mesh.res_param1));
        end
     end     
     
     tmpmtx2=zeros(input.num_mes,1);
 end
 
 % Occam
 if input.inv_flag==1 
     tmpmtx=zeros(mesh.num_param,1);
     tmpmtx2=fem.array_jacobian*log10(mesh.res_param1);
 end
%--------------------------------------------------------------------------

 
 
%---------------------Jt*dy------------------------------------------------
% Calculates RHS of inversion equation, dm1

% Most inversions
 if input.inv_flag~=5 && input.inv_flag~=6 
     dm1=fem.array_jacobian.'*input.Wd*(log10(input.real_data)-log10(fem.array_model_data) +tmpmtx2 ) - tmpmtx;
 % l1 data constraint
  else % difference inversions
     if itr>1 
         dm1=fem.array_jacobian.'*input.Wd*(log10(input.real_data)-log10(input.real_data0)-log10(fem.array_model_data)+log10(input.array_model_data_bgr)) - tmpmtx;
     else
         dm1=fem.array_jacobian.'*input.Wd*(log10(input.real_data)-log10(fem.array_model_data)) - tmpmtx;
     end
 end
 

% save('cntr_workspace')

dx1=dx1\dm1; % matrix inverse of dx1*dm1 solves minimisation equation
              % new dx1 is change in model state 
fem.dx1=dx1; % Keep this in case of Quasi Newton update


% Here calcualtes Gauss-Newton updates from my TGV modifications

% Modes: TV (-1), TGV (-3), l2-l2 (2), l1-l2 (-2) 
% note that l2-l2 overwrites original GN functionality, so that Wd is used
% consistently across all my inversions

if input.inv_flag == -1 || input.inv_flag == -3 || input.inv_flag == -2 || input.inv_flag == 2
    gamma_D = 5e-2;               % cutoff for |x| = 0
    gamma_C = 1e-3; 
    
    if input.inv_flag ~= -3
        mesh.px = zeros(size(mesh.res_param1));
        mesh.py = mesh.px;
    end
    
   
    tmp_res = mesh.res_param2;
    
    % Set i = 1 for normal operation!! i > 1 for separate IRLS and GN
    % iterations
    for i = 1    
            % Regulariser term weight - IRLS part
            if input.inv_flag == -1 || input.inv_flag == -3
                mesh.Rc = diag(1./sqrt( (mesh.cx*log10(tmp_res) - mesh.px).^2 + (mesh.cy*log10(tmp_res) - mesh.py).^2 + gamma_C.^2));
            else
                mesh.Rc = eye(length(mesh.res_param1));
            end
            
            if input.inv_flag == 2
                mesh.Rd = eye(length(input.real_data));
            else
                mesh.Rd = diag(1./sqrt((input.Wd*(log10(input.real_data) - log10(fem.array_model_data))).^2 + gamma_D^2));
            end
            
            mesh.ctc = mesh.cx'*mesh.Rc*(mesh.cx) + mesh.cy'*mesh.Rc*(mesh.cy);
            
            JTJ=fem.array_jacobian.'*input.Wd'*mesh.Rd*input.Wd*fem.array_jacobian;
 
            dm1 = fem.array_jacobian.'*input.Wd'*mesh.Rd*input.Wd*(log10(input.real_data)-log10(fem.array_model_data)) - tmp_lag*(mesh.cx'*mesh.Rc*(mesh.cx*log10(tmp_res)-mesh.px) +  mesh.cy'*mesh.Rc*(mesh.cy*log10(tmp_res)-mesh.py));

            dx1 = (JTJ + tmp_lag*mesh.ctc);
            
            dx1 = dx1\dm1;
            
            tmp_res = 10.^(log10(tmp_res) + dx1); % update tmp_res 
            
%             net_res(i) = sum(abs(tmp_res));

    end
    
%     figure(99)
%     plot(net_res)
    
    
    
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











