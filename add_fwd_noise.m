% add noise: adds gaussian noise to a .d datafile
percent_noise = 3 / 100;

input = struct('mes_in','0','res2d_flag',0);
input.mes_in = 'forward.d';

% First check for Loke
if input.res2d_flag==0
    %%%%%%%%%%%Classic%%%%%%%%%%%%%%%%%%%%
    tmp_d=importdata(input.mes_in);
elseif input.res2d_flag==1
    %%%%%%LOKE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    tmp_d=importdata(input.mes_in,' ',9);
    tmp_d=tmp_d.data(1:length(tmp_d.data)-1,2:10);
    ind1=find(tmp_d(:,9)>0);
    tmp_d=tmp_d(ind1,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add the noise
noisy_data = input.real_data + percent_noise*input.real_data.*randn(size(input.real_data));


% Save data again 
filename = 'noisy_forward.d';

fid = fopen(filename,'w');

for i = 1:length(input.ax)
    fprintf(fid,'%f %f %f %f %f %f %f %f %f \n',input.ax(i),input.az(i),input.bx(i),input.bz(i),input.mx(i),input.mz(i),input.nx(i),input.nz(i),noisy_data(i));
end

fclose(fid);

