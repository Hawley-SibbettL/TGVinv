% clearvars
% load('jac_files')

% jactmp=0;
% sum1=0;
% jam=0;
% jan=0;
% jbm=0;
% jbn=0;
% 
% % sum1 = zeros(input.num_mes,mesh.num_elements);
% jam = zeros(input.num_mes,max(max(mesh.param)));
% jan = jam;
% jbm = jan;
% jbn = jan;
% tic
% % mes = 1;
% % m = 1;
% % k = 1;
% 
% 
% for mes=1:input.num_mes 
% %     disp(num2str(mes))
%     for m=1:mesh.num_param
%         for k=1:mesh.no_elm_per_param(m)
% 
% 
% %          if (mesh.param(m,k)~=0)
% 
% % Takes each element inside the current parameter, identifies
% % node coordinates and calculates area
% l=mesh.param(m,k);
% 
% % a's are the indicies of the nodes in each element 
% a1=mesh.icon(1,l);    %I
% a2=mesh.icon(2,l);    %J
% a3=mesh.icon(3,l);    %L
% 
% b_node(1)=mesh.y_node_coord(a2)-mesh.y_node_coord(a3);
% b_node(2)=mesh.y_node_coord(a3)-mesh.y_node_coord(a1);
% b_node(3)=mesh.y_node_coord(a1)-mesh.y_node_coord(a2);
% 
% c_node(1)=mesh.x_node_coord(a3)-mesh.x_node_coord(a2);
% c_node(2)=mesh.x_node_coord(a1)-mesh.x_node_coord(a3);
% c_node(3)=mesh.x_node_coord(a2)-mesh.x_node_coord(a1);
% 
% area = (b_node(2)*c_node(3)-b_node(3)*c_node(2))/2;
% 
% % First loop (am)
% i=mesh.pa(mes);
% j=mesh.pm(mes);
% 
% % aaa is electrodes x nodes (potential @ each node for a source at each
% % electride pos???)
% s1=aaa(i,a1)+aaa(i,a2)+aaa(i,a3);
% s2=aaa(j,a1)+aaa(j,a2)+aaa(j,a3);
% sum1 = 1*(aaa(i,a1)*aaa(j,a1)+aaa(i,a2)*aaa(j,a2)+aaa(i,a3)*aaa(j,a3))+s1*s2;
% 
% sum1 = area*sum1*mesh.k(kit)^2/(12);
% 
% x_flux1=aaa(i,a1)*b_node(1)+aaa(i,a2)*b_node(2)+aaa(i,a3)*b_node(3);
% x_flux2=aaa(j,a1)*b_node(1)+aaa(j,a2)*b_node(2)+aaa(j,a3)*b_node(3);
% 
% y_flux1=aaa(i,a1)*c_node(1)+aaa(i,a2)*c_node(2)+aaa(i,a3)*c_node(3);
% y_flux2=aaa(j,a1)*c_node(1)+aaa(j,a2)*c_node(2)+aaa(j,a3)*c_node(3);
% 
% jam(mes,l)=(x_flux1*x_flux2+y_flux1*y_flux2)/(4*area)+sum1;
% 
% i=mesh.pa(mes);
% j=mesh.pn(mes);
% 
% s1=aaa(i,a1)+aaa(i,a2)+aaa(i,a3);
% s2=aaa(j,a1)+aaa(j,a2)+aaa(j,a3);
% sum1=1*(aaa(i,a1)*aaa(j,a1)+aaa(i,a2)*aaa(j,a2)+aaa(i,a3)*aaa(j,a3))+s1*s2;
% 
% sum1=area*sum1*mesh.k(kit)^2/(12);
% 
% x_flux1=aaa(i,a1)*b_node(1)+aaa(i,a2)*b_node(2)+aaa(i,a3)*b_node(3);
% x_flux2=aaa(j,a1)*b_node(1)+aaa(j,a2)*b_node(2)+aaa(j,a3)*b_node(3);
% 
% y_flux1=aaa(i,a1)*c_node(1)+aaa(i,a2)*c_node(2)+aaa(i,a3)*c_node(3);
% y_flux2=aaa(j,a1)*c_node(1)+aaa(j,a2)*c_node(2)+aaa(j,a3)*c_node(3);
% 
% jan(mes,l)=(x_flux1*x_flux2+y_flux1*y_flux2)/(4*area)+sum1;
% 
% i=mesh.pb(mes);
% j=mesh.pm(mes);
% 
% s1=aaa(i,a1)+aaa(i,a2)+aaa(i,a3);
% s2=aaa(j,a1)+aaa(j,a2)+aaa(j,a3);
% sum1=1*(aaa(i,a1)*aaa(j,a1)+aaa(i,a2)*aaa(j,a2)+aaa(i,a3)*aaa(j,a3))+s1*s2;
% 
% sum1=area*sum1*mesh.k(kit)^2/(12);
% 
% x_flux1=aaa(i,a1)*b_node(1)+aaa(i,a2)*b_node(2)+aaa(i,a3)*b_node(3);
% x_flux2=aaa(j,a1)*b_node(1)+aaa(j,a2)*b_node(2)+aaa(j,a3)*b_node(3);
% 
% y_flux1=aaa(i,a1)*c_node(1)+aaa(i,a2)*c_node(2)+aaa(i,a3)*c_node(3);
% y_flux2=aaa(j,a1)*c_node(1)+aaa(j,a2)*c_node(2)+aaa(j,a3)*c_node(3);
% 
% jbm(mes,l)=(x_flux1*x_flux2+y_flux1*y_flux2)/(4*area)+sum1;
% 
% i=mesh.pb(mes);
% j=mesh.pn(mes);
% 
% s1=aaa(i,a1)+aaa(i,a2)+aaa(i,a3);
% s2=aaa(j,a1)+aaa(j,a2)+aaa(j,a3);
% sum1=1*(aaa(i,a1)*aaa(j,a1)+aaa(i,a2)*aaa(j,a2)+aaa(i,a3)*aaa(j,a3))+s1*s2;
% 
% sum1=area*sum1*mesh.k(kit)^2/(12);
% 
% x_flux1=aaa(i,a1)*b_node(1)+aaa(i,a2)*b_node(2)+aaa(i,a3)*b_node(3);
% x_flux2=aaa(j,a1)*b_node(1)+aaa(j,a2)*b_node(2)+aaa(j,a3)*b_node(3);
% 
% y_flux1=aaa(i,a1)*c_node(1)+aaa(i,a2)*c_node(2)+aaa(i,a3)*c_node(3);
% y_flux2=aaa(j,a1)*c_node(1)+aaa(j,a2)*c_node(2)+aaa(j,a3)*c_node(3);
% 
% jbn(mes,l)=(x_flux1*x_flux2+y_flux1*y_flux2)/(4*area)+sum1;
% 
%         end
%     end
% end
% 
% t1 = toc;
% % ans1 = sum1;
% % jam1 = jam;
% % jan1 = jan;
% % jbn1=jbn;
% % jbm1 = jbm;
% disp(['t1 = ',num2str(t1)])
% 
% jsum = (jam-jbm-jan+jbn);
% 
% clear jan jam jbn jbm
% 
% 
% jsum = bsxfun(@rdivide, jsum, mesh.prop(1:max_grid_ind).^2);
% 
% % So far jacobian has been calculated per element, but we want per
% % parameter - so sum over elements in each parameter
% trans_jac1 = zeros(input.num_mes,mesh.num_param);
% for i = 1:mesh.num_param
%     trans_jac1(:,i) = sum(jsum(:,nonzeros(mesh.param(i,:))),2);
% end
% 
% trans_jac_old = trans_jac1;
% save('jac_test','trans_jac1')

%% 
clearvars

load('jac_files')
% This version of the jacobian calculation aims to be the vectorised
% representation which uses the least memory, and hence is a trade off
% between speed and memory. Currently it requires 2 copies of the largest 
% matrix, requiring 2*(input.num_mes*max_grid_ind)*8bytes of RAM
% 
% Note - when there are mixes of 4,3,2 point measurements may become less
% efficient, esp. if largest type doesn't dominate

tic
% maximum index of finite element inside the parameter grid 
% (excludes elements outside the parameter grid)
max_grid_ind = max(max(mesh.param));

b_node1 = mesh.y_node_coord(mesh.icon(2,1:max_grid_ind)) - mesh.y_node_coord(mesh.icon(3,1:max_grid_ind));
b_node2 = mesh.y_node_coord(mesh.icon(3,1:max_grid_ind)) - mesh.y_node_coord(mesh.icon(1,1:max_grid_ind));
b_node3 = mesh.y_node_coord(mesh.icon(1,1:max_grid_ind)) - mesh.y_node_coord(mesh.icon(2,1:max_grid_ind));

c_node1 = mesh.x_node_coord(mesh.icon(3,1:max_grid_ind)) - mesh.x_node_coord(mesh.icon(2,1:max_grid_ind));
c_node2 = mesh.x_node_coord(mesh.icon(1,1:max_grid_ind)) - mesh.x_node_coord(mesh.icon(3,1:max_grid_ind));
c_node3 = mesh.x_node_coord(mesh.icon(2,1:max_grid_ind)) - mesh.x_node_coord(mesh.icon(1,1:max_grid_ind));

% Each element of area is the area of each element in order given by icon
area = (b_node2.*c_node3 - b_node3.*c_node2)/2;

% Calculate the values for sources at at each position (if required) 
% and sum for the net potential
jsum = zeros(input.num_mes,max_grid_ind);

if (sum(mesh.pa.*mesh.pm)~=0)    
    ab_elec = mesh.pa; mn_elec = mesh.pm;
    sum1 = (aaa(ab_elec,mesh.icon(1,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(3,1:max_grid_ind))) + (aaa(ab_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind))) .* (aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(mn_elec,mesh.icon(3,1:max_grid_ind)));
    sum1 = bsxfun(@times,area',sum1).*mesh.k(kit)^2/(12);
    jsum = jsum + (bsxfun(@rdivide,(bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),b_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),b_node3')) + (bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),c_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),c_node3')),4*area') + sum1);   
end

if (sum(mesh.pa.*mesh.pn)~=0)    
    ab_elec = mesh.pa; mn_elec = mesh.pn;
    sum1 = (aaa(ab_elec,mesh.icon(1,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(3,1:max_grid_ind))) + (aaa(ab_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind))) .* (aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(mn_elec,mesh.icon(3,1:max_grid_ind)));
    sum1 = bsxfun(@times,area',sum1).*mesh.k(kit)^2/(12);
    jsum = jsum - (bsxfun(@rdivide,(bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),b_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),b_node3')) + (bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),c_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),c_node3')),4*area') + sum1);   
end


if (sum(mesh.pb.*mesh.pm)~=0)   
    ab_elec = mesh.pb; mn_elec = mesh.pm;
    sum1 = (aaa(ab_elec,mesh.icon(1,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(3,1:max_grid_ind))) + (aaa(ab_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind))) .* (aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(mn_elec,mesh.icon(3,1:max_grid_ind)));
    sum1 = bsxfun(@times,area',sum1).*mesh.k(kit)^2/(12);
    jsum = jsum - (bsxfun(@rdivide,(bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),b_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),b_node3')) + (bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),c_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),c_node3')),4*area') + sum1);
end

if (sum(mesh.pb.*mesh.pn)~=0)
    ab_elec = mesh.pb; mn_elec = mesh.pn;
    sum1 = (aaa(ab_elec,mesh.icon(1,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(3,1:max_grid_ind))) + (aaa(ab_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind))) .* (aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(mn_elec,mesh.icon(3,1:max_grid_ind)));
    sum1 = bsxfun(@times,area',sum1).*mesh.k(kit)^2/(12);
    jsum = jsum + (bsxfun(@rdivide,(bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),b_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),b_node3')) + (bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),c_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),c_node3')),4*area') + sum1);
end

clear sum1

jsum = bsxfun(@rdivide, jsum, mesh.prop(1:max_grid_ind).^2);

% So far jacobian has been calculated per element, but we want per
% parameter - so sum over elements in each parameter
trans_jac1 = zeros(input.num_mes,mesh.num_param);
for i = 1:mesh.num_param
    trans_jac1(:,i) = sum(jsum(:,nonzeros(mesh.param(i,:))),2);
end

t2 = toc;
disp(['t2 = ',num2str(t2)])

trans_jac_slow_vec = trans_jac1;
save('jac_test','trans_jac_slow_vec','-append')
%% Avoid repeated calc - fast_vec
clearvars
load('jac_files')

% This verion of the jacobean calcultion is the vectorised form which aims
% maximises speed at the cost of memory. It currently requires 11 copies of 
% the largest matrix, requiring 11*(input.num_mes*max_grid_ind)*8bytes 
% of RAM.
%
% Note - when there are mixes of 4,3,2 point measurements may become less
% efficient, esp. if largest type doesn't dominate

tic

% maximum index of finite element inside the parameter grid 
% (excludes elements outside the parameter grid)
max_grid_ind = max(max(mesh.param));

b_node1 = mesh.y_node_coord(mesh.icon(2,1:max_grid_ind)) - mesh.y_node_coord(mesh.icon(3,1:max_grid_ind));
b_node2 = mesh.y_node_coord(mesh.icon(3,1:max_grid_ind)) - mesh.y_node_coord(mesh.icon(1,1:max_grid_ind));
b_node3 = mesh.y_node_coord(mesh.icon(1,1:max_grid_ind)) - mesh.y_node_coord(mesh.icon(2,1:max_grid_ind));

c_node1 = mesh.x_node_coord(mesh.icon(3,1:max_grid_ind)) - mesh.x_node_coord(mesh.icon(2,1:max_grid_ind));
c_node2 = mesh.x_node_coord(mesh.icon(1,1:max_grid_ind)) - mesh.x_node_coord(mesh.icon(3,1:max_grid_ind));
c_node3 = mesh.x_node_coord(mesh.icon(2,1:max_grid_ind)) - mesh.x_node_coord(mesh.icon(1,1:max_grid_ind));

% Each element of area is the area of each element in order given by icon
area = (b_node2.*c_node3 - b_node3.*c_node2)/2;

% Precalculates values and then calculates jsum. Order is chosen so as to
% minimise memory usage.

% Calculate the values for sources at at each position (if required) 
% and sum for the net potential
jsum = zeros(input.num_mes,max_grid_ind);

% a precalculations
if sum(mesh.pa) ~= 0
    s1_a = aaa(mesh.pa,mesh.icon(1,1:max_grid_ind)) + aaa(mesh.pa,mesh.icon(2,1:max_grid_ind)) + aaa(mesh.pa,mesh.icon(3,1:max_grid_ind));
    xflux_a = bsxfun(@times,aaa(mesh.pa,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(mesh.pa,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(mesh.pa,mesh.icon(3,1:max_grid_ind)),b_node3');
    yflux_a = jsum + bsxfun(@times,aaa(mesh.pa,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(mesh.pa,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(mesh.pa,mesh.icon(3,1:max_grid_ind)),c_node3');
end

% m precalculations - named mn as n will overwrite later ofr efficency 
if sum(mesh.pm) ~= 0
    s2_mn = aaa(mesh.pm,mesh.icon(1,1:max_grid_ind)) + aaa(mesh.pm,mesh.icon(2,1:max_grid_ind)) + aaa(mesh.pm,mesh.icon(3,1:max_grid_ind));
    xflux_mn = bsxfun(@times,aaa(mesh.pm,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(mesh.pm,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(mesh.pm,mesh.icon(3,1:max_grid_ind)),b_node3');
    yflux_mn = bsxfun(@times,aaa(mesh.pm,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(mesh.pm,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(mesh.pm,mesh.icon(3,1:max_grid_ind)),c_node3');
end

% am calculation
if sum(mesh.pa.*mesh.pm)~=0
    sum1 = (aaa(mesh.pa,mesh.icon(1,1:max_grid_ind)).*aaa(mesh.pm,mesh.icon(1,1:max_grid_ind))+aaa(mesh.pa,mesh.icon(2,1:max_grid_ind)).*aaa(mesh.pm,mesh.icon(2,1:max_grid_ind))+aaa(mesh.pa,mesh.icon(3,1:max_grid_ind)).*aaa(mesh.pm,mesh.icon(3,1:max_grid_ind))) + s1_a.*s2_mn;
    sum1 = bsxfun(@times,area',sum1).*mesh.k(kit)^2/(12);
    jsum = (bsxfun(@rdivide,(xflux_a.*xflux_mn + yflux_a.*yflux_mn),4*area') + sum1);
end

% b precalculations
if sum(mesh.pb) ~=0
    s1_b = aaa(mesh.pb,mesh.icon(1,1:max_grid_ind)) + aaa(mesh.pb,mesh.icon(2,1:max_grid_ind)) + aaa(mesh.pb,mesh.icon(3,1:max_grid_ind));
    xflux_b = bsxfun(@times,aaa(mesh.pb,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(mesh.pb,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(mesh.pb,mesh.icon(3,1:max_grid_ind)),b_node3');
    yflux_b = bsxfun(@times,aaa(mesh.pb,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(mesh.pb,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(mesh.pb,mesh.icon(3,1:max_grid_ind)),c_node3');
end

% bm calculations
if sum(mesh.pb.*mesh.pm) ~= 0
    sum1 = (aaa(mesh.pb,mesh.icon(1,1:max_grid_ind)).*aaa(mesh.pm,mesh.icon(1,1:max_grid_ind))+aaa(mesh.pb,mesh.icon(2,1:max_grid_ind)).*aaa(mesh.pm,mesh.icon(2,1:max_grid_ind))+aaa(mesh.pb,mesh.icon(3,1:max_grid_ind)).*aaa(mesh.pm,mesh.icon(3,1:max_grid_ind))) + s1_b.*s2_mn;
    sum1 = bsxfun(@times,area',sum1).*mesh.k(kit)^2/(12);
    jsum = jsum - (bsxfun(@rdivide,(xflux_b.*xflux_mn + yflux_b.*yflux_mn),4*area') + sum1);
end

% n precalculations - called mn as overwrites m for efficency
if sum(mesh.pn) ~= 0
    s2_mn = aaa(mesh.pn,mesh.icon(1,1:max_grid_ind)) + aaa(mesh.pn,mesh.icon(2,1:max_grid_ind)) + aaa(mesh.pn,mesh.icon(3,1:max_grid_ind));
    xflux_mn = bsxfun(@times,aaa(mesh.pn,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(mesh.pn,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(mesh.pn,mesh.icon(3,1:max_grid_ind)),b_node3');
    yflux_mn = bsxfun(@times,aaa(mesh.pn,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(mesh.pn,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(mesh.pn,mesh.icon(3,1:max_grid_ind)),c_node3');
end

% an calculations
if sum(mesh.pa.*mesh.pn) ~= 0
    sum1 = (aaa(mesh.pa,mesh.icon(1,1:max_grid_ind)).*aaa(mesh.pn,mesh.icon(1,1:max_grid_ind))+aaa(mesh.pa,mesh.icon(2,1:max_grid_ind)).*aaa(mesh.pn,mesh.icon(2,1:max_grid_ind))+aaa(mesh.pa,mesh.icon(3,1:max_grid_ind)).*aaa(mesh.pn,mesh.icon(3,1:max_grid_ind))) + s1_a.*s2_mn;
    clear s1_a
    sum1 = bsxfun(@times,area',sum1).*mesh.k(kit)^2/(12);
    jsum = jsum - (bsxfun(@rdivide,(xflux_a.*xflux_mn + yflux_a.*yflux_mn),4*area') + sum1);
end
clear s1_a xflux_a yflux_a

% bn calculations
if sum(mesh.pb.*mesh.pn)~=0
    sum1 = (aaa(mesh.pb,mesh.icon(1,1:max_grid_ind)).*aaa(mesh.pn,mesh.icon(1,1:max_grid_ind))+aaa(mesh.pb,mesh.icon(2,1:max_grid_ind)).*aaa(mesh.pn,mesh.icon(2,1:max_grid_ind))+aaa(mesh.pb,mesh.icon(3,1:max_grid_ind)).*aaa(mesh.pn,mesh.icon(3,1:max_grid_ind))) + s1_b.*s2_mn;
    clear s1_b s2_mn
    sum1 = bsxfun(@times,area',sum1).*mesh.k(kit)^2/(12);
    jsum = jsum + (bsxfun(@rdivide,(xflux_b.*xflux_mn + yflux_b.*yflux_mn),4*area') + sum1);
end


% Don't clear everything as should be working within limits of RAM.
% Clearing s matrices means will have space fro trans_jac1, plus staying 1
% matrix below peak usage if possible
clear s1_b s2_mn

jsum = bsxfun(@rdivide, jsum, mesh.prop(1:max_grid_ind).^2);

% So far jacobian has been calculated per element, but we want per
% parameter - so sum over elements in each parameter
trans_jac1 = zeros(input.num_mes,mesh.num_param);
for i = 1:mesh.num_param
    trans_jac1(:,i) = sum(jsum(:,nonzeros(mesh.param(i,:))),2);
end

t3 = toc;
disp(['t3 = ',num2str(t3)])

trans_jac_fast_vec = trans_jac1;
save('jac_test','trans_jac_fast_vec','-append')

%%
clearvars
load('jac_files')

jactmp=0;
sum1=0;
jam=0;
jan=0;
jbm=0;
jbn=0;

% aaa=fem.tmp_aaa;

trans_jac1=zeros(input.num_mes,mesh.num_param);

% save('jac_files')

tic
%/*  printf(" JAC_CONTROL..\n");*/
for mes=1:input.num_mes
    %     disp(num2str(mes))
    for m=1:mesh.num_param
        for k=1:mesh.no_elm_per_param(m)
            %          if (mesh.param(m,k)~=0)
            
            % Takes each element inside the current parameter, identifies
            % node coordinates and calculates area
            l=mesh.param(m,k);
            
            a1=mesh.icon(1,l);    %I
            a2=mesh.icon(2,l);    %J
            a3=mesh.icon(3,l);    %L
            
            b_node(1)=mesh.y_node_coord(a2)-mesh.y_node_coord(a3);
            b_node(2)=mesh.y_node_coord(a3)-mesh.y_node_coord(a1);
            b_node(3)=mesh.y_node_coord(a1)-mesh.y_node_coord(a2);
            
            c_node(1)=mesh.x_node_coord(a3)-mesh.x_node_coord(a2);
            c_node(2)=mesh.x_node_coord(a1)-mesh.x_node_coord(a3);
            c_node(3)=mesh.x_node_coord(a2)-mesh.x_node_coord(a1);
            
            area=(b_node(2)*c_node(3)-b_node(3)*c_node(2))/2;
            
            
            if ((mesh.pa(mes)*mesh.pm(mes))~=0)
                
                i=mesh.pa(mes);
                j=mesh.pm(mes);
                
                s1=aaa(i,a1)+aaa(i,a2)+aaa(i,a3);
                s2=aaa(j,a1)+aaa(j,a2)+aaa(j,a3);
                sum1=1*(aaa(i,a1)*aaa(j,a1)+aaa(i,a2)*aaa(j,a2)+aaa(i,a3)*aaa(j,a3))+s1*s2;
                
                sum1=area*sum1*mesh.k(kit)^2/(12);
                
                x_flux1=aaa(i,a1)*b_node(1)+aaa(i,a2)*b_node(2)+aaa(i,a3)*b_node(3);
                x_flux2=aaa(j,a1)*b_node(1)+aaa(j,a2)*b_node(2)+aaa(j,a3)*b_node(3);
                
                y_flux1=aaa(i,a1)*c_node(1)+aaa(i,a2)*c_node(2)+aaa(i,a3)*c_node(3);
                y_flux2=aaa(j,a1)*c_node(1)+aaa(j,a2)*c_node(2)+aaa(j,a3)*c_node(3);
                
                jam=(x_flux1*x_flux2+y_flux1*y_flux2)/(4*area)+sum1;
                
            end
            
            
            if((mesh.pa(mes)*mesh.pn(mes))~=0)
                
                i=mesh.pa(mes);
                j=mesh.pn(mes);
                
                
                s1=aaa(i,a1)+aaa(i,a2)+aaa(i,a3);
                s2=aaa(j,a1)+aaa(j,a2)+aaa(j,a3);
                sum1=1*(aaa(i,a1)*aaa(j,a1)+aaa(i,a2)*aaa(j,a2)+aaa(i,a3)*aaa(j,a3))+s1*s2;
                
                sum1=area*sum1*mesh.k(kit)^2/(12);
                
                x_flux1=aaa(i,a1)*b_node(1)+aaa(i,a2)*b_node(2)+aaa(i,a3)*b_node(3);
                x_flux2=aaa(j,a1)*b_node(1)+aaa(j,a2)*b_node(2)+aaa(j,a3)*b_node(3);
                
                y_flux1=aaa(i,a1)*c_node(1)+aaa(i,a2)*c_node(2)+aaa(i,a3)*c_node(3);
                y_flux2=aaa(j,a1)*c_node(1)+aaa(j,a2)*c_node(2)+aaa(j,a3)*c_node(3);
                
                jan=(x_flux1*x_flux2+y_flux1*y_flux2)/(4*area)+sum1;
                
                
            end
            
            if((mesh.pb(mes)*mesh.pm(mes))~=0)
                
                i=mesh.pb(mes);
                j=mesh.pm(mes);
                
                
                s1=aaa(i,a1)+aaa(i,a2)+aaa(i,a3);
                s2=aaa(j,a1)+aaa(j,a2)+aaa(j,a3);
                sum1=1*(aaa(i,a1)*aaa(j,a1)+aaa(i,a2)*aaa(j,a2)+aaa(i,a3)*aaa(j,a3))+s1*s2;
                
                sum1=area*sum1*mesh.k(kit)^2/(12);
                
                x_flux1=aaa(i,a1)*b_node(1)+aaa(i,a2)*b_node(2)+aaa(i,a3)*b_node(3);
                x_flux2=aaa(j,a1)*b_node(1)+aaa(j,a2)*b_node(2)+aaa(j,a3)*b_node(3);
                
                y_flux1=aaa(i,a1)*c_node(1)+aaa(i,a2)*c_node(2)+aaa(i,a3)*c_node(3);
                y_flux2=aaa(j,a1)*c_node(1)+aaa(j,a2)*c_node(2)+aaa(j,a3)*c_node(3);
                
                jbm=(x_flux1*x_flux2+y_flux1*y_flux2)/(4*area)+sum1;
                
            end
            
            if((mesh.pb(mes)*mesh.pn(mes))~=0)
                
                i=mesh.pb(mes);
                j=mesh.pn(mes);
                
                s1=aaa(i,a1)+aaa(i,a2)+aaa(i,a3);
                s2=aaa(j,a1)+aaa(j,a2)+aaa(j,a3);
                sum1=1*(aaa(i,a1)*aaa(j,a1)+aaa(i,a2)*aaa(j,a2)+aaa(i,a3)*aaa(j,a3))+s1*s2;
                
                sum1=area*sum1*mesh.k(kit)^2/(12);
                
                x_flux1=aaa(i,a1)*b_node(1)+aaa(i,a2)*b_node(2)+aaa(i,a3)*b_node(3);
                x_flux2=aaa(j,a1)*b_node(1)+aaa(j,a2)*b_node(2)+aaa(j,a3)*b_node(3);
                
                y_flux1=aaa(i,a1)*c_node(1)+aaa(i,a2)*c_node(2)+aaa(i,a3)*c_node(3);
                y_flux2=aaa(j,a1)*c_node(1)+aaa(j,a2)*c_node(2)+aaa(j,a3)*c_node(3);
                
                jbn=(x_flux1*x_flux2+y_flux1*y_flux2)/(4*area)+sum1;
                
            end
            jactmp=jactmp+(jam-jbm-jan+jbn);
            
            %         end
            
        end
        
        trans_jac1(mes,m)=jactmp/(mesh.prop(l)^2); %4
        jactmp=0;
    end
end
t1 = toc;
disp(['t1 = ',num2str(t1)])

trans_jac_old = trans_jac1;
save('jac_test','trans_jac1','-append')

%% DEFINITELY WORKS, DO NOT TOUCH!!!!
clearvars
load('jac_files')
tic
% This is the last element index contained in the parameter grid (not
% extending to inf)
max_grid_ind = max(max(mesh.param));

% find coordonates of element nodes projected onto unit element
b_node1 = mesh.y_node_coord(mesh.icon(2,1:max_grid_ind)) - mesh.y_node_coord(mesh.icon(3,1:max_grid_ind));
b_node2 = mesh.y_node_coord(mesh.icon(3,1:max_grid_ind)) - mesh.y_node_coord(mesh.icon(1,1:max_grid_ind));
b_node3 = mesh.y_node_coord(mesh.icon(1,1:max_grid_ind)) - mesh.y_node_coord(mesh.icon(2,1:max_grid_ind));

c_node1 = mesh.x_node_coord(mesh.icon(3,1:max_grid_ind)) - mesh.x_node_coord(mesh.icon(2,1:max_grid_ind));
c_node2 = mesh.x_node_coord(mesh.icon(1,1:max_grid_ind)) - mesh.x_node_coord(mesh.icon(3,1:max_grid_ind));
c_node3 = mesh.x_node_coord(mesh.icon(2,1:max_grid_ind)) - mesh.x_node_coord(mesh.icon(1,1:max_grid_ind));

% Gives the areas of all elements as ordered in mesh.icon 
area = (b_node2.*c_node3 - b_node3.*c_node2)/2;

% Perform separate calculations  
ab_elec = mesh.pa; mn_elec = mesh.pm;
sum1 = (aaa(ab_elec,mesh.icon(1,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(3,1:max_grid_ind))) + (aaa(ab_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind))) .* (aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(mn_elec,mesh.icon(3,1:max_grid_ind)));
sum1 = bsxfun(@times,area',sum1).*mesh.k(kit)^2/(12);
jam = bsxfun(@rdivide,(bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),b_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),b_node3')) + (bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),c_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),c_node3')),4*area') + sum1;

ab_elec = mesh.pa; mn_elec = mesh.pn;
sum1 = (aaa(ab_elec,mesh.icon(1,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(3,1:max_grid_ind))) + (aaa(ab_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind))) .* (aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(mn_elec,mesh.icon(3,1:max_grid_ind)));
sum1 = bsxfun(@times,area',sum1).*mesh.k(kit)^2/(12);
jan = bsxfun(@rdivide,(bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),b_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),b_node3')) + (bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),c_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),c_node3')),4*area') + sum1;

ab_elec = mesh.pb; mn_elec = mesh.pm;
sum1 = (aaa(ab_elec,mesh.icon(1,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(3,1:max_grid_ind))) + (aaa(ab_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind))) .* (aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(mn_elec,mesh.icon(3,1:max_grid_ind)));
sum1 = bsxfun(@times,area',sum1).*mesh.k(kit)^2/(12);
jbm = bsxfun(@rdivide,(bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),b_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),b_node3')) + (bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),c_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),c_node3')),4*area') + sum1;

ab_elec = mesh.pb; mn_elec = mesh.pn;
sum1 = (aaa(ab_elec,mesh.icon(1,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind)).*aaa(mn_elec,mesh.icon(3,1:max_grid_ind))) + (aaa(ab_elec,mesh.icon(1,1:max_grid_ind))+aaa(ab_elec,mesh.icon(2,1:max_grid_ind))+aaa(ab_elec,mesh.icon(3,1:max_grid_ind))) .* (aaa(mn_elec,mesh.icon(1,1:max_grid_ind))+aaa(mn_elec,mesh.icon(2,1:max_grid_ind))+aaa(mn_elec,mesh.icon(3,1:max_grid_ind)));
sum1 = bsxfun(@times,area',sum1).*mesh.k(kit)^2/(12);
jbn = bsxfun(@rdivide,(bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),b_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),b_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),b_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),b_node3')) + (bsxfun(@times,aaa(ab_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(ab_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(ab_elec,mesh.icon(3,1:max_grid_ind)),c_node3')).*(bsxfun(@times,aaa(mn_elec,mesh.icon(1,1:max_grid_ind)),c_node1') + bsxfun(@times,aaa(mn_elec,mesh.icon(2,1:max_grid_ind)),c_node2') + bsxfun(@times,aaa(mn_elec,mesh.icon(3,1:max_grid_ind)),c_node3')),4*area') + sum1;

% combine potential measurements, divide by resistivity (chargability)
trans_jac1_elm = bsxfun(@rdivide, (jam-jbm-jan+jbn), mesh.prop(1:max_grid_ind).^2);

% Sum elements inside each parameter to get parameter dependent jacobian
trans_jac1 = zeros(input.num_mes,mesh.num_param);
for i = 1:mesh.num_param
    trans_jac1(:,i) = sum(trans_jac1_elm(:,nonzeros(mesh.param(i,:))),2);
end
toc
trans_jac_new1 = trans_jac1;
save('jac_test','trans_jac_new1','-append')

