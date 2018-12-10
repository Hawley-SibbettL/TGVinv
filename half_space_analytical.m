% Analytical half space solution (for I = 1)
function analytical_hs = half_space_analytical(final,rho)
% First read in electrode data
input = struct('mes_in',final.mes_in,'res2d_flag',0);

input = readindata(input);


% Calculate data (assumes no topography as half space)
analytical_hs = 0.*input.real_data;

AN = abs(input.ax - input.nx);
AM = abs(input.ax - input.mx);
BN = abs(input.bx - input.nx);
BM = abs(input.bx - input.mx);

analytical_hs = (rho/(2*pi)) * (AM.^(-1) - AN.^(-1) - BM.^(-1) + BN.^(-1));

% Now plot against inverted results (ordered in order of magnitude for
% convenieance)
figure
[~, ind] = sort(analytical_hs(:,end));   % Picks last iteration!!!

plot(analytical_hs(ind),'k-^','markersize',3)
plot(final.res_param1(ind),'r-v','markersize',3)


end




% Reads in datafile
function [input] = readindata(input)

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

end