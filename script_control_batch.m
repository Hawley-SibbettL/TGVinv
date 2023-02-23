% Script control for l2/TV inversion
clear all
close all
% Input data files and settings here

% Filename and pathname of file
source_files = {'gb2_long2_1pc.dat'}; % cell array of filenames (strings) to invert
 % cell array of filenames (strings) to invert
pathname = ['D:\TGV_revision\thesis\mesh_define\res2dmod files']; %['D:\data\2d_quarry_data'];
out_path = ['D:\TGV_revision\thesis\periodic\golden_test\l1 solution init\reg mesh\wd repeat\strong_cool 1'];
% source_files = {'clifton.inv', 'ingham.inv', 'finningley.inv'}; % cell array of filenames (strings) to invert
% pathname = ['D:\data\2d_quarry_data'];
% out_path = ['D:\TGV_revision\thesis\periodic\'];
append_txt = ''; % append string to filename

res2d_flag = 2; % 0=IP4DI input, 1 = Loke, 2 = loke dd, 3 = resipy
wd_flag = 2; % 0 - no data errors in inversion. 1 - use errors directly from input file 2 - use % error given by 'noise dev'
noise_dev = 0.01; % percentage error applied to all measurements if used
data_type = 1; % 1 = apparent resistivity, 2 = resistance
res2d_headerlines = 12; % number of headerlines to discount in res2d file. [synth/imager = 6/12]
m_init_flag = 4; % Initialisation of m: [1 - interpolated pseudosection, loaded in from res2dmod files][2 - l2 iteration from homogeneous half space] [3 - vertically faulted half space] [4 read in solution from previous inversion. currently hardcoded]

% inv_flag_list: vector of flags specifying inversion types.
% -1 = l1 data term, l1 (TV) regularisation
% -2 = l1 data term, l2 regularisation
% 2 = l2 data term, l2 regularisation
inv_flag_list = [-1, 2];% currently disbled l1 data term %[-1,-2, 2];
lagrn_in = [0.05]; % initial damping parameter, lambda
itn = 14;% number of iterations in inversion

% Advanced Options
ref_flag = 0;  % defines cell width. 0 - 1 cell per electrode. 1 - 2 cells per electrode (default)
manual_flag = 1; % 0 - default. 1 - mesh spacing and depth can be customised by modifying depth_calculator.m
ext_mesh_flag = 1; % 0 - default. 1 - extends mesh outside
line_search_flag = 0; % 0 - cooling off scheme for lambda. 1 - line search along lambda direction. 2 - line search along dm direction. 3 - line search along both 4 - golden section search
ls_limit_flag = 1; % 0 - lambda can be small; 1 - lambda limited to (0.01?) minimum a la Loke
lagrn_ls = [];%[0.0001, 0.0005, 0.001, 0.003, 0.005, 0.007, 0.01, 0.03, 0.05, 0.07, 0.11, 0.15, 0.2]; % lagrn values for line search, [] is bisection
target_decrease = 0; % target decrease in rms during line search (0 for minimiser)
lagrn_w = [1];
line_search_obj_flag = 0; % calculates and saves objective functions during line search
l1_data_flag = 0; % 0 - l2 data term, 1 - l1 data term
diff_type = 1;% 1 - forward-backward. 2 - central. 3 - central half points (m_{i+1/2} - m_{i-1/2}), requires uniform mesh 4 - half point attempt 2
n_ext = 2; % number of extra cells to add outside x range. 2 is standard
ghost_flag = 0; % 1 - use neumann_bc ghost points method, 0 - implicit neumann bc
dm_bc = 1; % order of neumann bc in m if ghost points are used
periodic_bc = 1; % = 1 for periodic bc, moderates diff_type =1 only
diff_weight_flag = 0; % = 1 to weight finited differences by cell spacing
% Output filename details
inv_text = {'_l1', '_l2','_l22'};

%% Inversion process begins
for k = 1:length(lagrn_in)
    lagrn = lagrn_in(k);
    for i = 1:length(source_files)
        for j = 1:length(inv_flag_list) % Different inv parameters
            
            
            if length(lagrn) > 1
                input.lagrn = lagrn(i);
            else
                input.lagrn = lagrn;
            end
            

           lagr_txt = {['_ls0','lag',append_txt]};
            
            
            % --- Experimental options, advise against use. (default = 0) ---
            input.lc_flag = 0; % l curve method to determine damping parameter
            damp_flag = 0; % attempt to damp
            input.topography_flag = 0; % My flag for reading topography straight from file.

            % ---------- Set parameters from ip4di unpackaging, do not change
            
            input.damp_surface_flag = damp_flag;
            input.ext_mesh_flag = ext_mesh_flag; % 0 - default. 1 - extends mesh outside
            input.line_search_flag = line_search_flag;
            input.lagrn_ls = lagrn_ls;
            input.lagrn_w = lagrn_w;
            input.dm_bc = dm_bc;
            input.line_search_obj_flag = line_search_obj_flag;
            input.diff_type = diff_type;
            input.ghost_flag = ghost_flag;
            input.target_decrease = target_decrease;
            input.periodic_bc = periodic_bc;
            input.noise_dev = noise_dev;
            input.n_ext = n_ext;
            input.diff_weight_flag = diff_weight_flag;
            input.m_init_flag = m_init_flag;
            input.ls_limit_flag = ls_limit_flag;
            
            if ref_flag == 1
                manual_flag = 1;
                input.refine_mesh_flag = 1;  % = 1 for cell widths equal to half electride spacing
                
            else
                input.refine_mesh_flag = 0;
            end
            
            
            % --------------------- Fixed Parameters ----------------------------------
            
            % I am always using DC data
            input.sip_flag=0;
            input.ip_flag=0;
            input.dc_flag=1;
            input.time_lapse_flag = 0;    % Time lapse inversion not implemented for TGV
            
            % Background file load settings
            input.bgr_res_flag = 0; % 0 if no background file, 1 otherwise
            input.pathname = '[]';
            input.filename = '[]'; % *.mod file for difference *.mat for timelapse. Why??
            input.par_nm = fullfile(input.pathname, input.filename);
            input.fd_weight = 1;    % Default = 1. Set to equal 0 turns off fd weightings.
            input.l1_data_flag = l1_data_flag;
            input.p_init_flag = [];
            input.p_conv = [];
            % ------------------- Carry out read data and meshing ---------------------
            input.res2d_flag = res2d_flag;
            input.res2d_headerlines = res2d_headerlines;
            input.wd_flag = wd_flag;
            input.data_type = data_type;
            filename = cell2mat(source_files(i));
            
            [filename_path, short_filename,file_ext] = fileparts(filename);
            % Allows path to be given in filename
            if ~isempty(filename_path)
                pathname2 = fullfile(pathname,filename_path);
            else
                pathname2 = pathname;
            end
            input.batch_save = [short_filename,cell2mat(inv_text(j)),cell2mat(lagr_txt)];
            input.inv_flag = inv_flag_list(j);
            input.out_path = out_path;
            input.mes_in = fullfile(pathname2, [short_filename, file_ext]);
            
            if input.time_lapse_flag==0
                [input]=read_data(input);
                disp('Data read...')
            else
                [input]=read_4d_data(input);
                disp('data read....')
            end
            
            [input,mesh]=create_mesh3_script(input,manual_flag);
            disp('..mesh done.')
            
            % Need axis handle to be 'big part'. Part of bad unpacking of IP4DI GUI
            global big_part texxt; big_part = gca;
            
            
            %% -------------------- Inversion Parameters ------------------------------
            % Inversion flag determines type of inversion:
            % SINGLE TIME:
            % -2=IRLS GN (l2 data); -1=IRLS GN (l1 data); 0=Levenburgh-Marquadt;
            % 1=Occam; 2=GN;
            % TIME LAPSE:
            % 7=4D
            
            input.lagrn_reduction = 2; % Factor by which lagrange multiplier is divided each itr
            input.itn = itn; % Number of iterations to perform
            %     input.inv_flag = -1;
            % --------------------- Fixed Parameters ----------------------------------
            % (shouldn't need to change often for my applications)
            
            % boundary conditions for FEM (always hom. dirichlet for internal boundary)
            mesh.bc_option = 2; % 1=air/surface is dirichlet; 2=mixed (air/surface neumann)
            input.acb_flag = 0; % Not using active constraint balancing
            input.lagrn_min = 0.005; % Minimum lagrange multiplier value (stops cooling off)
            input.lagrn_max = 1; % Max lagrange multiplier (acb only)
            input.limit_res = 0; % =1 if want to limit max/min res values, 0 otherwise
            input.min_res = 1; % minimum allowed res (default = 0)
            input.max_res=100000000; % Max allowed res (default = 100000000)
            input.jacobian_flag = 2; % 2 for full Jacobian; 1 for Quasi-Newton
            input.atc_flag = 0; % = 1 for active time constraints only.
            
            % Parameters set by inversion_parameters.m on buttonpress
            input.elec_array_type=1;
            input.conv_rate = 1;
            input.current=1;
            input.strike_dist=0;
            input.loke_colors=0;
            
            
            % Need a textbox textt to write in to let main work
            figure(1)
            texthand = uicontrol('Style','Text','Units','Normalized','Position',[0.1,0.1,0.8,0.8]);
            global textt
            textt = texthand;
            figure(2)
            big_part = gca;
            
            % ------------------------ Carry out inversion ----------------------------
            
            mesh.gamma_c = 1e-12;
            mesh.gamma_d = 1e-15;
            
            
            % Alternate initialisations
            if input.m_init_flag == 1
                [~, nm, ~] = fileparts(filename);
                if input.refine_mesh_flag == 1
                    load( fullfile(pathname,'pseudosections',nm) )
                else
                    load( fullfile(pathname,'pseudosections',[nm, '_reg']) )
                end
                mesh.m_init = pseudosection; % pseudosection values
            elseif input.m_init_flag == 4
                % l2 solution
%                 base_file = load('D:\TGV_revision\thesis\periodic\golden_test\reg mesh\gb2_long2_1pc_l2_ls0lag','final');
%                 mesh.m_init = base_file.final.res_param1( :, end );
                % l1 solution
%                 base_file = load('D:\TGV_revision\thesis\periodic\golden_test\l2 itr1 init\reg mesh\step1_cool\gb2_long2_1pc_l1_ls0lag','final');
%                 mesh.m_init = base_file.final.res_param1( :, end ); 
                % homogeneous
                mesh.m_init = log10(mean(input.real_data)*ones(size(mesh.res_param1)));
            elseif input.m_init_flag == 3 % split half space
                mesh.m_init = repmat(mesh.mean_res, [size(mesh.param_x)]);
                mesh.m_init(mesh.param_x > floor(max(mesh.param_x)/2)) = 10.^( log10(mesh.mean_res) + 0.1*log10(mesh.mean_res) ); 
                mesh.m_init(mesh.param_x < floor(max(mesh.param_x)/2)) = 10.^( log10(mesh.mean_res) - 0.1*log10(mesh.mean_res) ); 
            end
            
            main(input,mesh);
            
            
            %         catch
            %             disp('error in inversion... move on to next dataset')
            %         end
            clear input mesh textt big_part texhand fem
            close all
            
            
        end
        
    end
end