function [input]=read_data(input)



% First check for Loke
if input.res2d_flag==0
%%%%%%%%%%%Classic%%%%%%%%%%%%%%%%%%%%
tmp_d=importdata(input.mes_in);
elseif input.res2d_flag==1
%%%%%%LOKE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
try
    disp('tab delimited...?')
    tmp_d=importdata(input.mes_in,'\t',input.res2d_headerlines);
    tmp_d1 = tmp_d.data;
catch
%     try
    disp('space delimited...?')
    tmp_d=importdata(input.mes_in,' ',input.res2d_headerlines);
    tmp_d1 = tmp_d.data;
%     catch
%         error('read data error occured')
%     end
   
end
[nani, nanj] = find(isnan(tmp_d1));   % Nans in incomplete rows
tmp_d1(nani,:) = [];                % Remove incomplete rows
tmp_d=tmp_d1(:,2:end);  

% Ignore topgraphy (otherwise triggers boreholes)
if input.topography_flag == 0
    tmp_d(:,[2,4,6,8]) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif input.res2d_flag == 2    % Loke dipole dipole
    tmp_d1 = importdata(input.mes_in,' ',6); 
    tmp_d1 = tmp_d1.data;
    [nani, nanj] = find(isnan(tmp_d1));   % Nans in incomplete rows
    tmp_d1(nani,:) = [];                % Remove incomplete rows
    x = tmp_d1(:,1);              % Midpoint of each dd structure
    a = tmp_d1(:,2);                    % intra-dipole distance
    n = tmp_d1(:,3);                     % a*n = inter dipole distance
    
    % Now need to go from dipole-dipole format to IP4DI
    tmp_d = zeros(size(tmp_d1,1),9);     % z co-ordinates all zeros 
    tmp_d(:,1) = x - (a.*n)./2;
    tmp_d(:,3) = x - (a.*n)./2 - a;
    tmp_d(:,5) = x + (a.*n)./2;
    tmp_d(:,7) = x + (a.*n)./2 + a; 
    tmp_d(:,9) = abs(tmp_d1(:,4));
    if size(tmp_d1,2) >9
        tmp_d(10:size(tmp_d1,2)) = tmpd1(10:end);
    end
    % round to avoid problems
%     tmp_d = round(10000*tmp_d)./10000;
elseif input.res2d_flag == 3 % resipy input - data stored inside named folder
    [in_path, in_file, in_ext] = fileparts(input.mes_in);
    r2_folder = fullfile(in_path, in_file);
    
    el_pos = importdata(fullfile(r2_folder, 'electrodes.dat'));
    el_pos(:, 2) = []; % remove y dimension as 2d
    
    tmp_d1 = importdata(fullfile(r2_folder, 'R2_forward.dat'),' ',1);
    tmp_d1 = tmp_d1.data;
    tmp_d = [el_pos(tmp_d1(:,2),:), el_pos(tmp_d1(:,3),:), el_pos(tmp_d1(:,4),:), el_pos(tmp_d1(:,5),:), tmp_d1(:,7)];
end




% This is the common part
input.ax=tmp_d(:,1);
input.az=abs(tmp_d(:,2));
input.bx=(tmp_d(:,3));
input.bz=abs(tmp_d(:,4));
input.mx=tmp_d(:,5);
input.mz=abs(tmp_d(:,6));
input.nx=tmp_d(:,7);
input.nz=abs(tmp_d(:,8));


input.real_data=tmp_d(:,9);


input.num_mes=length(input.ax);



%%%%%%%%%%%%%%%%%%%%%%%%%
%%%IP DATA or SIP%%%%%%
try
    input.ip_data=tmp_d(:,10); 
    
    
    if input.sip_flag==1
        disp('SIP data found');
        input.ip_num=1;
        input.real_data=complex(input.real_data,input.ip_data);
        % if are in mag and phase
%         real_part=input.real_data.*(cos(input.ip_data./1000));
%         imag_part=input.real_data.*(sin(input.ip_data./1000));
%         input.real_data=complex(real_part,imag_part);
         
    elseif input.ip_flag==1
        disp('IP data found');
        input.ip_num=2;
        input.ip_data=input.ip_data./1000;% Loke has it in msec
    else
        input.stdev_error=input.ip_data; % Currently on fow DC
        input.ip_num=1;
    end
catch
    disp('No IP or SIP data found');
    input.ip_num=1;
    input.stdev_error=ones(input.num_mes,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Here check for IP or SIP and STDEV
% if input.sip_flag==1 || input.ip_flag==1
%     try
%         input.stdev_error=tmp_d(:,11); 
%     catch
%         input.stdev_error=ones(input.num_mes,1);
%     end
% end


clear tmp_d;







        

% Try to find number of electrodes used




%--------------------------------------------------------------------------
    % search for how many electrodes
    Error=0;
    b_el=0;
    n_el=0;
    input.array_type=1;
    ax=length(unique(input.ax));
    ay=length(unique(input.az));
    az=length(unique(input.az));
    
    bx=length(unique(input.bx));
    by=length(unique(input.bz));
    bz=length(unique(input.bz));

    
    mx=length(unique(input.mx));
    my=length(unique(input.mz));
    mz=length(unique(input.mz));

    
    nx=length(unique(input.nx));
    ny=length(unique(input.nz));
    nz=length(unique(input.nz));
    
    if ax==1 && ay==1 && az==1
        disp('Error. A electrode must be there ALWAYS');
        Error=1;
    end
    
    if mx==1 && my==1 && mz==1
        disp('Error. M electrode must be there ALWAYS');
        Error=1;
    end
    
    if bx==1 && by==1 && bz==1
        b_el=1;
        Error=0;
    end
    
    if nx==1 && ny==1 && nz==1
        n_el=1;        
        Error=0;
    end
    
    
   if b_el==1 && n_el==0
       disp('Pole Dipole array?');
       input.elec_array_type=2;
   elseif b_el==1 && n_el==1
       disp('Pole Pole array?');
       input.elec_array_type=3;
   elseif b_el==0 && n_el==1
       disp('Not tested array');
       pause
   else
       disp('4 electrode array?');
       input.elec_array_type=1;
   end



if Error==1
    pause
end

      
input.m_array_type(1:input.num_mes,1)=input.elec_array_type;


  
end