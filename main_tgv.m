% t_fwd is time for forward problem, t_inv is time solving matrix equation
function [input, mesh, final, fem,ex,t_fwd,t_inv] = main_tgv(input,mesh,tgv_itr,fem)
ex=0;
final=[];
t_inv = 0;

if tgv_itr == 1
    
    %     if matlabpool('size') == 0
    %         matlabpool open
    %     end
    
    % mesh.Rd = 1;
    
    % fem=[]; % initialize fem as empty matrix
    
    
    if input.wd_flag == 1
        if input.data_type == 2
            input.stdev_error = input.stdev_error.*mesh.geofac;
        end
        input.stdev_error(input.stdev_error == 0) = 0.5*min(input.stdev_error(input.stdev_error>0)); % ensures no infinities
        input.Wd=diag(1./input.stdev_error);
    elseif input.wd_flag == 2
        input.Wd = spdiags(input.real_data*input.noise_dev, 0, length(input.real_data), length(input.real_data));
    else
        input.Wd = diag(ones(length(input.real_data),1));
    end
    
    
    if isfield(input,'lc_flag') == 0  % creates an lc flag if not present
        input.lc_flag = 0;
        disp('...added lc_flag = 0..')
    elseif input.lc_flag == 1 && ~(input.inv_flag==-2 || input.inv_flag==-1 || input.inv_flag==2)
        disp( ' L-curve method not available for current inversion technique, using cooling-off method instead')
        input.lc_flag = 0;
    end
    
    disp(['input.inv_flag = ', num2str(input.inv_flag)])
    if input.lc_flag == 0 ; disp(['input.lagrn = ', num2str(input.lagrn)]); end
    
    
    % Calculate k and g
    [mesh.k,mesh.g]=rwaven(mesh.probe_spacing,mesh.probe_spacing*mesh.max_n);
    % rdis=double([mesh.probe_spacing:mesh.probe_spacing:mesh.max_n*mesh.probe_spacing]);
    % k_init= logspace(-2,0.3,5);
    % [mesh.k,mesh.g,obj,err]=xu_inversion(rdis,k_init);
end

% global textt

% set(textt,'String','****** INVERSION  STARTS *****'); drawnow update;
clear global big part
clear global textt



if input.time_lapse_flag==0
    for ip_cnt=1:input.ip_num
        for itr = tgv_itr
            tic
            % skips first section after first sp1 itr to avoid
            % double-calculating
            % fem, convergence checks
            fem=mes_control_fast(itr,input,mesh,fem,0);
            t_fwd = toc;
            
            [ex,fem,mesh]=rms(itr,ip_cnt,input,mesh,fem); % rms_tgv

        end
        % New itr of inversion
        tic
        
        
        if input.m_init_flag ~= 2 && itr == 1 % skip itr 1 inversion and insert initial model
            mesh.res_param1 = mesh.m_init;
            final=info_save(0,itr,ip_cnt,input,mesh,fem,final);
            mesh.ls = 0;
        else
            t_inv = toc;
            
            
            if input.lc_flag == 1
                [input,mesh] = l_curve(mesh,fem,input,itr);
            elseif input.line_search_flag == 0
                [input]=update_lagran(itr,ip_cnt,1,input);
            end
            final=info_save(0,itr,ip_cnt,input,mesh,fem,final);
            if (ex==1 || ex==2) ;
                if itr ~= input.itn+1;
                    mesh.res_param1 = mesh.res_param2;
                end
                final.itr=itr-1;
                break;
            end         % Exit as per RMS results          
            [mesh,fem,input]=invert_cntr(itr,1,ip_cnt,input,mesh,fem);
            t_inv = toc;
        end
        [fem,mesh,input]=prop_change(itr,input,mesh,fem);
        %             save('itr_results','mesh','input')
        auto_contour(1,itr,ip_cnt,input,mesh,fem); colormap cool;
        
        
        if itr==input.itn+1; final.itr=itr;end
        if input.ip_flag==1
            [input,mesh]=ip_calc(input,mesh);
        else
            break;
        end
    end
    
elseif input.time_lapse_flag==1
    for ip_cnt=1:input.ip_num
        for itr=1:input.itn
            [fem,mesh]=d4_prepare_data(itr,input,mesh,fem);
            [ex,fem,mesh]=rms_4d(itr,ip_cnt,input,mesh,fem);
            final=info_save(0,itr,ip_cnt,input,mesh,fem,final);
            if (ex==1 || ex==2) ;final.itr=itr-1;break; end
            [input]=update_lagran(itr,ip_cnt,1,input);
            [mesh]=kim_inversion2(input,mesh,fem);
            auto_contour_4d(itr,input,mesh,fem);
        end
        if itr==input.itn+1; final.itr=itr;end
        if input.ip_flag==1
            [input,mesh]=ip_calc(input,mesh);
        else
            break;
        end
    end
end

% Here save outputs...
info_save(3,itr,ip_cnt,input,mesh,fem,final);
% save('inv_results.mat','final');


end