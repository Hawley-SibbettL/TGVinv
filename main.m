function main(input,mesh)
ex=0;

% if matlabpool('size') == 0
%        matlabpool open
% end

% mesh.Rd = 1;


global textt
fem=[]; % initialize fem as empty matrix
final=[];

mesh=smooth_mtx_surface4(input,mesh);set(textt,'String','********* SMOOTH MATRIX ********');drawnow update;

if input.wd_flag == 1
    if input.data_type == 2
        input.stdev_error = input.stdev_error.*mesh.geofac;
    end
    input.stdev_error(input.stdev_error == 0) = 0.5*min(input.stdev_error(input.stdev_error>0)); % ensures no infinities
    input.Wd=diag(1./input.stdev_error);
elseif input.wd_flag == 2
    input.Wd = spdiags(input.real_data*input.noise_dev, 0, length(input.real_data), length(input.real_data));
else
    input.Wd = speye(length(input.real_data));
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

set(textt,'String','****** INVERSION  STARTS *****'); drawnow update;

if input.time_lapse_flag==0
    for ip_cnt=1:input.ip_num
        for itr=1:input.itn+1
            % fem, convergence checks
            tic
            fem = mes_control_fast(itr,input,mesh,fem,1);
            toc
            [ex,fem,mesh]=rms(itr,ip_cnt,input,mesh,fem);

            % New itr of inversion
            %             if (input.inv_flag == -1) || (input.inv_flag == -2);
            %                 mesh = IRLS(input, fem, mesh);
            % %                 disp('IRLS')
            %             end
            
            if input.m_init_flag ~= 2 && itr == 1 % skip itr 1 inversion and insert initial model                
                mesh.res_param1 = mesh.m_init;
                mesh.Rc = speye(length(mesh.res_param1));
            else % normal inversion
                if input.lc_flag == 1 %&& itr > 1
                    [input,mesh] = l_curve(mesh,fem,input,itr);
                elseif input.line_search_flag == 0
                    [input]=update_lagran(itr,ip_cnt,1,input);
                end
                final=info_save(0,itr,ip_cnt,input,mesh,fem,final);
                if (ex==1 || ex==2) ;final.itr=itr-1; break; end
                [mesh,fem,input]=invert_cntr(itr,0,ip_cnt,input,mesh,fem);
            end
            [fem,mesh,input]=prop_change(itr,input,mesh,fem);
            auto_contour(1,itr,ip_cnt,input,mesh,fem); colormap cool;
        end
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
if isfield(input,'batch_save') == 0 
    save('inv_results.mat','final');
else
    curr_dir = pwd;    
    if ~exist(input.out_path, 'dir')
        mkdir(input.out_path);
    end
    cd(input.out_path);
    save(input.batch_save,'final')
    cd(curr_dir);
end

end