% Script control for a total generalised variation inversion
clearvars
close all

% Input data files and settings here

% Filename and pathname of file
source_files = {'clifton.inv'}; % cell array of filenames (strings) to invert
pathname = 'C:\Users\Luke\Documents\MATLAB\ip4Di_TGV_PAPER\TGV data\2d_quarry_data';
out_path = 'C:\Users\Luke\Documents\MATLAB\ip4Di_TGV_PAPER\TGV results\no_ghost\w75_real';
append_txt = ''; % used to add additional text tFo the output filename

% default file format is IP4DI files (.d)- see IP4DI manual50
% can also take .dat filetypes as generated from res2dinv
res2d_flag = 1; % 0 - .d file; 1 - .dat file, general format. 2 - .dat file dipole-dipole format
wd_flag = 0; % 0 - no data errors in inversion. 2 - use errors from input file
data_type = 1; % 1 = apparent resistivity, 2 = resistance
res2d_headerlines = 12; % manually input number of headerlines in .dat file. [synth/imager = 6/12]
p_init_flag = 2; % choice of how to initialise p in sp2. [0 is p = 0] [1 is p = grad(m)] [2 is result of an l2 iteration strting from p = 2]
p_conv = 0.001; % 0.01 before
mu = {2}; % TGV trade off parameter. Cell array allows multiple values
lagrn_in = [0.05]; % damping parameter (lambda)
itn = 10; % number of gauss newton iterations in inversion

% Advanced options
ref_flag = 0;  % defines cell width. 0 - 1 cell per electrode. 1 - 2 cells per electrode (default)
manual_flag = 1; % 0 - default. 1 - mesh spacing and depth can be customised by modifying depth_calculator.m
extend_mesh_flag = 1; % 0 - default. 1 - extends mesh outside
line_search_flag = 3; % 0 - cooling off scheme for lambda. 1 - line search along lambda direction. 2 - line search along dm direction. 3 - line search along both
lagrn_ls = [0.0005, 0.001, 0.003, 0.005, 0.007, 0.01, 0.03, 0.05, 0.07, 0.11, 0.15, 0.2]; % lagrn values for line search
lagrn_w = [0.75]; % w_lagrn values along line search (fraction of lambda to travel in search direction)
line_search_obj_flag = 0; % calculates and saves objective functions during line search (if selected)
dm_bc = 1; % order of dm neumann bc
l1_data_flag = 0; % 0 - l2 data term, 1 - l1 data term
diff_type = 1;% 1 - forward-backward. 2 - central. 3 - central, symmetrical (m_{i+1/2} - m_{i-1/2})
n_ext = 8; % number of extra edge cells. 2 usually.
ghost_flag = 0; % = 1 to use ghost points
%% Inversion process begins
curr_dir = pwd;

for k = 1:length(lagrn_in)
    lagrn = lagrn_in(k);
    
    for q = 1:length(mu)
        
        for r = 1:length(source_files)
            clearvars -except q r i source_files ghost_flag n_ext res2d_headerlines pathname out_path curr_dir lagr_txt lagrn mu ref_flag damp_surface_flag append_txt res2d_flag wd_flag data_flag itn manual_flag extend_mesh_flag data_type p_init_flag p_conv k lagrn_in lagrn line_search_flag dm_bc lagrn_ls lagrn_w l1_data_flag line_search_obj_flag diff_type
            close all
            
            
            
            % --------- Initialisation, do not change ------------------------
            % Includes read data, meshing, as well as legacy parameter
            % assignments from IP4DI base code
            
            % Experimental options, advise against use
            input.lc_flag = 0; % l curve method to determine damping parameter
            damp_surface_flag = 0; % attempt to damp
            input.topography_flag = 0; % My flag for reading topography straight from file.
            input.diff_type = diff_type;
            input.n_ext = n_ext; 
            
            input.tgv_lagrn = mu{q};
            input.lagrn = lagrn;
            
            input.res2d_headerlines = res2d_headerlines;
            input.damp_surface_flag = damp_surface_flag;
            input.out_path = out_path;
            
            % Filename and pathname of file- .d for single .dat TL
            filename = source_files{r};
            %     pathname = source_path(i);
            [filename_path, short_filename,file_ext] = fileparts(filename);
            
            if ~isempty(filename_path)
                pathname2 = fullfile(pathname,filename_path);
            else
                pathname2 = pathname;
            end
            
            if line_search_flag == 0
                lagr_txt = [' lag',num2str(1000*input.lagrn,'%d'),'_tgv_lag',num2str(10*input.tgv_lagrn),append_txt];
            else
                lagr_txt = [' lagls',num2str(length(lagrn_ls)),'_tgv_lag',num2str(10*input.tgv_lagrn),append_txt];
            end
            input.batch_save = ['TGV ',short_filename,lagr_txt];
            
            mkdir(input.out_path,input.batch_save); % output folder
            
            % ---Read Data Parameters---------
            input.wd_flag = wd_flag;
            input.time_lapse_flag = 0; % TGV Time lapse not implemented
            input.data_type = data_type;    % 1 = apparent resistivity, 2 = resistance
            input.res2d_flag = res2d_flag;
            input.fd_weight = 1;    % = 1 for true fd weighting for cell separation. always use
            input.ext_mesh_flag = extend_mesh_flag;
            input.tgv_itr = itn;
            input.line_search_flag = line_search_flag;
            input.lagrn_ls = lagrn_ls;
            input.lagrn_w = lagrn_w;
            input.dm_bc = dm_bc;
            input.l1_data_flag = l1_data_flag;
            input.p_init_flag = p_init_flag;
            input.p_conv = p_conv;
            input.line_search_obj_flag = line_search_obj_flag;
            input.ghost_flag = ghost_flag;
            % Bypasses borehole detection currently.
            if ref_flag == 1;
                input.refine_mesh_flag = 1;  % = 1 for cell widths equal to half electride spacing
                manual_flag = 1;       % flag for manual depth (need to custimise depth function)
            else
                input.refine_mesh_flag = 0;
            end
            
            % --- Fixed Parameters ------------
            
            % I am always using DC data
            input.sip_flag=0;
            input.ip_flag=0;
            input.dc_flag=1;
            
            % Background file load settings
            input.bgr_res_flag = 0; % 0 if no background file, 1 otherwise
            input.pathname = '[]';
            input.filename = '[]'; % *.mod file for difference *.mat for timelapse. Why??
            input.par_nm = fullfile(input.pathname, input.filename);
            
            % --- Carry out read data and meshing -----
            
            input.mes_in = fullfile(pathname2, [short_filename, file_ext]);
            
            %  try
            if input.time_lapse_flag==0
                [input]=read_data(input);
                disp('Data read...')
            else
                [input]=read_4d_data(input);
                disp('data read....')
            end
            
            [input,mesh]=create_mesh3_script(input,manual_flag);
            disp('..mesh done.')
            
            input.inv_flag = -3; % -3, flag for TGV inversion
            
            % initialise p
            mesh.px = zeros(mesh.num_param,1);
            mesh.py = mesh.px;
            
            
            % ---- Fixed input/mesh parameters ----
            
            % boundary conditions for FEM (always hom. dirichlet for internal boundary)
            input.bc_option = 2; % 1=air/surface is dirichlet; 2=mixed (air/surface neumann)
            input.acb_flag = 0; % Not using active constraint balancing
            input.lagrn_min = 0.005; % Minimum lagrange multiplier value (stops cooling off)
            input.lagrn_max = 1; % Max lagrange multiplier (acb only)
            input.limit_res = 0; % =1 if want to limit max/min res values, 0 otherwise
            input.min_res = 1; % minimum allowed res (default = 0)
            input.max_res=100000000; % Max allowed res (default = 100000000)
            input.jacobian_flag = 2; % 2 for full Jacobian; 1 for Quasi-Newton
            input.atc_flag = 0; % = 1 for active time constraints only.
            input.conv_rate = 1; % very strict convergence checker to prevent termination in main
            
            % temp variables to measure convergence
            rms1 = 0;
            rms2 = 0;
            
            % Parameters set by inversion_parameters.m on buttonpress
            input.data_type=1;
            input.elec_array_type=1;
            input.conv_rate=2;
            input.current=1;
            input.strike_dist=0;
            input.loke_colors=0;
            
            % Need a textbox textt to write in to let main work
            figure(1)
            texthand = uicontrol('Style','Text','Units','Normalized','Position',[0.1,0.1,0.8,0.8]);
            global textt big_part
            textt = texthand;
            figure(2)
            big_part = gca;
            
            % Computes gradient operators
            mesh=smooth_mtx_surface4(input,mesh);
            
            fem = struct();
            %         mesh=smooth_mtx_surface4(input,mesh);
            
            % ---------- TGV outer loop begins ------------------------------
            
            %         starting_res = mesh.res_param1;
            %     lagrn = 0.15;
            itr = 1;
            
            
            while itr <= input.tgv_itr
                
                
                
                % -------- SUBPROBLEM 1 --------
                % Fix p, update resistivty model
                
                
                % Done in main_tgv instead
%                      input = update_lagran(itr,0,1,input);
                
                
                input.itn = 1; % Perform only one main_tgv itr
                input.lc_flag = 0;
                
                % Calls single Gauss-Newton iteration
                [input, mesh, final, fem, ex, t_fwd, t_inv] = main_tgv(input,mesh,itr,fem);
                
                
                disp(['Forward problem time = ',num2str(t_fwd),'s ; sp1 matrix equation time = ',num2str(t_inv),'s'])
                
                cd(fullfile(input.out_path,input.batch_save));
                save(['tgv_it_', num2str(itr), '_sp1'],'final')
                cd(curr_dir);
                

                
                disp(['*** Iteration ', num2str(itr),' ***'])
                
                
                % -------- SUBPROBLEM 2 --------
                % Fix model, update p
                tic
                if itr == input.tgv_itr; break; end
                
                if input.ghost_flag == 1
                % add ghost points to apply bc to m
                [num_param, domain_ind, param_x, param_y, res_param] = neumann_bc(mesh.num_param, mesh.param_x, mesh.param_y, mesh.res_param1, [0, input.dm_bc]);
                else
                    num_param = mesh.num_param;
                    domain_ind = 1:num_param;
                    param_x = mesh.param_x;
                    param_y = mesh.param_y;
                    res_param = mesh.res_param1;
                end 
                
                grad_mx = mesh.cx*log10(res_param);
                grad_my = mesh.cy*log10(res_param);
                input.sp2_itr = 100;
                
                %initialise p - includes bc ghost point terms
                if input.p_init_flag == 0 || input.p_init_flag == 2
                    mesh.px = zeros(mesh.num_param,1);
                    mesh.py = mesh.px;
                elseif input.p_init_flag == 1 % p = grad(m) on domain of m.
                    mesh.px = mesh.cx*res_param;
                    mesh.px = mesh.px(domain_ind);
                    mesh.py = mesh.cy*res_param;
                    mesh.py = mesh.py(domain_ind);
                end
                
                
                % Used in main calculation
                gamma_p = 1e-15;               % cutoff for |x| -> 0
                gamma_c = 1e-15;
                % Used in performance metrics
                p1 = zeros(2*mesh.num_param,1);
                rms_p = zeros(1,input.sp2_itr+1);
                
                % main sp2 loop
                if itr ~= input.tgv_itr
                    for i = 2:(input.sp2_itr+1)
                        
                        
                        if input.ghost_flag == 1
                        % add ghost points to p. Ax p,m have same domain, cx, cy are
                        % the same, but different bc (first order homogeneous neumann).
                        [~, ~, ~, ~, p_tmp] = neumann_bc(mesh.num_param, mesh.param_x, mesh.param_y, [mesh.px, mesh.py], [0, 1]);
                        px_tmp = p_tmp(:,1);
                        py_tmp = p_tmp(:,2);
                        else
                            px_tmp = mesh.px;
                            py_tmp = mesh.py;
                        end

                        
                        % calculate weights
                        if input.p_init_flag == 2 && i == 2 % l2 1st iteration as initialisation
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
                        
                        % Domain_ind first
                        a11 = a11(domain_ind, domain_ind);
                        a12 = a12(domain_ind, domain_ind);
                        a21 = a21(domain_ind, domain_ind);
                        a22 = a22(domain_ind, domain_ind);
                        b1 = b1(domain_ind);
                        b2 = b2(domain_ind);       
                        
                        %                     A = [a11, a12; a21, a22];
                        %                     b = [b1; b2];
                        A = [a11, a12; a21, a22];
                        b = [b1; b2];
                        %                     condition = cond(A);
                        %                     disp(['itr: ', num2str(i-1), ' | cond(A): ', num2str(condition,'%.2g')])
                        p2 = A\b;
                        %                     p2 = p2(domain_ind);
                        clear A b b1 b2 a11 a12 a21 a22
                        
                        % spearate px, py
                        mesh.px = p2(1:length(p2)/2);
                        mesh.py = p2((length(p2)/2+1):length(p2));
                        % eliminate ghost points
%                         mesh.px = mesh.px(domain_ind);
%                         mesh.py = mesh.py(domain_ind);
                        p2 = [mesh.px; mesh.py];
                        
                        % Store results
                        px(:,i) = mesh.px;
                        py(:,i) = mesh.py;
                        
                        % root mean square difference between last and current
                        % solution normalise by rms p value
                        %                     rms_p(i) = sqrt(mean((p2 - p1).^2))/(0.5*sqrt(mean(p2.^2 + p1.^2));
                        
                        
                        
                        rms_p(i-1) = sqrt(mean((p2 - p1).^2))/(0.5*sqrt(mean( (mesh.cx*log10(res_param)).^2 + (mesh.cy*log10(res_param)).^2)));
                        %         disp(['rms dp = ',num2str(p_rms(i))])
                        %                 rms_change(i) = abs((p_rms(i) - p_rms(i-1)))./(p_rms(i) + p_rms(i-1)); % percentage rms change per itr
                        
                        disp(['tgv iteration ',num2str(i-1),' : rms dp = ',num2str(rms_p(i-1))])
                        
                        
                        if ((i > 4) && (rms_p(i-1) < input.p_conv)) %  stop if change smaller than p_conv after first few itr
                            sp2_i = i;  % store number of iterations needed
                            break
                        elseif (i > 2 && (rms_p(i-1) > rms_p(i-2))) % stop if diverging
                            sp2_i = i-1;  % store number of iterations needed
                            mesh.px = p1(1:length(p1)/2);
                            mesh.py = p1((length(p1)/2+1):length(p1));
                            break
                        end
                        
                        
                        if i == (input.sp2_itr+1)
                            sp2_i = i;
                        end
                        
                        p1 = p2;

                    end
                end
                
                
                t_sp2 = toc;
                disp(['Subproblem 2 time = ',num2str(t_sp2),'s (for ',num2str(sp2_i),' iterations)'])
                
                rms_p(i:end) = [];
                
                % Plot convergence
                figure(10)
                hold off
                plot(rms_p)
                title('rms change in p')
                xlabel('sp 2 iteration')
                ylabel('rms dp')
                
                
                % Save to final
                final.px = mesh.px;
                final.py = mesh.py;
                final.t_fwd(itr) = t_fwd;
                final.t_inv(itr) = t_inv;
                final.t_sp2 = t_sp2;
                final.sp2_itr = sp2_i;
                final.rms_p = rms_p;
                final.cx = mesh.cx;
                final.cy = mesh.cy;
                final.tgv_lagrn = input.tgv_lagrn;
                
                save(fullfile(input.out_path,input.batch_save,['tgv_it_', num2str(itr), '_sp2']),'final')
                itr = itr + 1;
                clear p_tmp px_tmp py_tmp p1 p2 rms_p Rm Rp param_x param_y res_param grad_mx grad_my
                
            end
            
        end
        
        
        
    end
    
end



