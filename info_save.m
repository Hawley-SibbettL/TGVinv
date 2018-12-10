    function final=info_save(stf,itr,ip_cnt,input,mesh,fem,final)
       
if stf==0
    final.num_param=mesh.num_param;
    final.num_of_bor=mesh.num_of_bor;
    final.param_x=mesh.param_x;
    final.param_y=mesh.param_y;
    % final.borehole_no_surface_flag=input.borehole_no_surface_flag;
    final.itn=itr;
    final.lagrn(itr)=input.lagrn;
    final.inv_flag=input.inv_flag;
    final.bgr_res_flag=input.bgr_res_flag;
    % final.par_nm=input.par_nm;
    final.acb_flag=input.acb_flag;
    final.mes_in=input.mes_in;
    final.num_mes=input.num_mes;
    final.time_lapse_flag=input.time_lapse_flag;
    final.ip_flag=input.ip_flag;
    final.sip_flag=input.sip_flag;
    final.dc_flag=input.dc_flag;
 
    
    if input.lc_flag == 1 && itr > 1
        final.lc_curvature(itr-1,:) = mesh.lc_curvature;            
        final.lc_misfit(itr-1,:) = mesh.lc_misfit;
        final.lc_roughness(itr-1,:) = mesh.lc_roughness;
        final.lagr_ind(itr-1) = mesh.lagr_ind;
%         final.lc_curvature_smooth(itr-1,:) = mesh.lc_curvature_smooth;
        if itr == 2;  final.lc_trial_lagr =  mesh.lc_trial_lagr; end
        final.lc_lagr(itr-1) = mesh.lc_lagr;
%         final.lc_peak_ind{itr-1} = mesh.lc_peak_ind';
    end
    
    
    if input.bgr_res_flag~=0 &&input.time_lapse_flag==0  && itr>1;final.bgr_res_param=mesh.bgr_param;end

    if input.time_lapse_flag==1 ;final.num_files=input.num_files; end
    
    if itr>2 && input.acb_flag==1 
        final.ACB(:,itr-1)=fem.L1;
    end

    if itr>2 && input.time_lapse_flag==0 ;final.resolution_1=fem.resolution_matrix; end
    
    % Keep models
    final.res_param1(:,itr)=mesh.res_param2;
    final.acb(:,itr)=1;
    final.array_model_data=fem.array_model_data2;
    final.RMS(itr)=fem.rms_sum2;


    if input.ip_flag==1&& ip_cnt==1
        final.res_model=mesh.res_param2;
        final.res_final_data=final.array_model_data;
        final.all_res_model=final.res_param1;
        final.res_rms=final.RMS;
    end


    if input.ip_flag==1&& ip_cnt==2    
        final.chargeb=1000*(mesh.res_param2-mesh.res_final_param)./(mesh.res_param2);
        final.ip_model=final.res_param1;
        final.ip_rms(itr)=fem.rms_sum2;
    end

    if input.time_lapse_flag==1
        final.d4_res_param1(:,:,itr)=mesh.d4_res_param2;
    end
    
    if itr == 1
       final.half_space_jac = fem.array_jacobian;
    end
end






if stf==3
%     write param and measurement	file
        if input.time_lapse_flag==0 %  I need on file for difference or backgroud inversion. All other files will be in mat form.
            datout=fopen('datout.inv','wt');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if input.bgr_res_flag==0
                fprintf(datout,'INV_TYPE= %d DATAFILE= %s RMS= %f ITR= %d LGRN= %f NUM_PARAM= %d NUM_MEAS= %d\n',input.inv_flag,input.mes_in,fem.rms_sum2,itr-1,input.lagrn,mesh.num_param,input.num_mes);
                if input.dc_flag==1
                    for i=1:mesh.num_param
                        fprintf(datout,'%f %f %f\n',mesh.param_x(i),-mesh.param_y(i),mesh.res_param2(i));
                    end
                    for i=1:input.num_mes
                        fprintf(datout,'%d %f %f %f %f %f %f %f %f %f %f\n',1,input.ax(i),input.az(i),input.bx(i),input.bz(i),input.mx(i),input.mz(i),input.nx(i),input.nz(i),...
                        input.real_data(i),fem.array_model_data(i));
                    end
                elseif input.ip_flag==1
                    % Not yet supported
                elseif input.sip_flag==1
                    for i=1:mesh.num_param
                        fprintf(datout,'%f %f %f %f\n',mesh.param_x(i),-mesh.param_y(i),abs(mesh.res_param2(i)),1000*atan2(imag(mesh.res_param2(i)),real(mesh.res_param2(i))));
                    end
                    for i=1:input.num_mes
                        fprintf(datout,'%d %f %f %f %f %f %f %f %f %f %f %f %f\n',1,input.ax(i),input.az(i),input.bx(i),input.bz(i),input.mx(i),input.mz(i),input.nx(i),input.nz(i),...
                        real(input.real_data(i)),imag(input.real_data(i)),real(fem.array_model_data(i)),imag(fem.array_model_data(i)));
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fclose(datout);
        end
          
end

    end