% This function aims to calculate the size of the penalty terms in the
% objective function for a model
% 
% It also has the option to generate a dataset from an input model and mesh
% flags are used to specify which functions are desired
%
% Inputs : mesh - generated using forward_modelling_2 using an input data
% file containing the measurements. 
%
% Model - this is the resistivity model to be used. Currently generated by
% reading mesh into res2dmod_generator and saving the interpolated output

function objective_calc()
clearvars
close all

% Read in model
filename = 'new_diag_block3';
user = 'bgsvisluke2';
obj_flag = 1; % = 1 if want to calculate objective functions
fwd_flag = 0; % = 1 if want to calculate forward problem and save data
p_plot = 1; % Value of mu to plot

filepath = ['C:\Users\',user,'\OneDrive - The University of Nottingham\UoN Box Migration\Data\New TGV test models\'];
load_model = load([filepath,filename]);
model = load_model.ip4di_im;
% model = 10.^repmat(1.3 + (0.5./13)*(1:13), 45, 1)'; % for loglin_v
param = reshape(model(1:13,:)',[],1); % same from as in ip4di. Cuts bottom layer off to match.
xdim = size(model,2);
ydim = size(model,1);

% Read in mesh, add model to mesh structure
mesh = []; input = [];
load('C:\Users\bgsvisluke2\Box Sync\Data\Sim2\fwd_mesh_mod.mat');
% load('gauss150_mesh');
mesh.res_param1(:) = param; mesh.init_res_param = param;  mesh.res_param2(:) = param; input.bgr_res_param = param;


% Calculate gradient matrix
[mesh.cx, mesh.cy] = calc_grad(mesh);


% Calculate the objective function terms
% TV term: |grad(m)|
grad2 = (mesh.cx*log10(mesh.res_param1)).^2 + (mesh.cy*log10(mesh.res_param1)).^2;
TV_reg = sum(sqrt(grad2(:)));
L2_reg = sqrt(sum(grad2(:)));

% TGV term - need to calculate p from minimisation problem once.
tgv_lagrn = [0.1:0.1:3];
px = zeros(length(mesh.res_param1),length(tgv_lagrn));
py = px;

if obj_flag == 1
    
    for i = 1:length(tgv_lagrn);
        input.tgv_lagrn = tgv_lagrn(i);
        mesh = TGV_sp2(input,mesh);
        TGV_t1_im = (sqrt((mesh.cx*log10(mesh.res_param1) - mesh.px).^2 + (mesh.cy*log10(mesh.res_param1) - mesh.py).^2));
        TGV_t2_im = (sqrt( (mesh.cx'*mesh.px).^2 + (mesh.cy'*mesh.py).^2 +  0.5*(mesh.cy'*mesh.py + mesh.cx'*mesh.px).^2));
        TGV_t1(i) = sum(TGV_t1_im(:));
        TGV_t2(i) = sum(TGV_t2_im(:));
        
        px(:,i) = mesh.px;
        py(:,i) = mesh.py;
        
        
        figure(5)
        surf(unique(mesh.param_x),unique(mesh.param_y),(reshape(TGV_t1_im,xdim,[])'),'edgecolor','none')
        view([0,0,1])
        title(['TGV \lambda_2 = ',num2str(input.tgv_lagrn)])
        set(gca,'ydir','reverse')
        colorbar
        axis image
        pause(0.01)
        if input.tgv_lagrn == 2
            TGV_t1_im_2 = reshape(TGV_t1_im,xdim,[])';
        end
    end
    
    
    TGV_reg = TGV_t1 + input.tgv_lagrn.*TGV_t2;
    p_im = reshape(sqrt(px(:,tgv_lagrn == p_plot).^2 + py(:,tgv_lagrn == p_plot).^2),xdim,[])';
    grad_im = reshape(sqrt(grad2),xdim,[])';
    
    % Print results
    disp([filename, ' objective terms (\lambda_2 = 2)'])
    disp(['L2 objective = ',num2str(L2_reg)])
    disp(['TV objective = ',num2str(TV_reg)])
    disp(['TGV objective = ',num2str(sum(TGV_reg.*(tgv_lagrn == 2)))])
    
    figure(1)
    plot(tgv_lagrn,TGV_reg)
    hold on
    plot(tgv_lagrn,TGV_t1)
    plot(tgv_lagrn,TGV_t2)
    plot(tgv_lagrn,TV_reg.*ones(1,length(tgv_lagrn)))
    plot(tgv_lagrn,L2_reg.*ones(1,length(tgv_lagrn)))
    legend('total TGV objective','TGV |\nabla m - p|_{l1}','TGV 0.5|\nabla p + \nabla p^T|_{l1}','|\nabla m|_{L1}','|\nabla m|_{L2}')
    title('TGV objective function')
    xlabel('\lambda_2')
    ylabel(['Objective function : ',filename])
    
    figure(2)
    surf(unique(mesh.param_x),unique(mesh.param_y),p_im,'edgecolor','none')
    view([0,0,1])
    title('p')
    set(gca,'ydir','reverse')
    colorbar
    axis image
    
    figure(3)
    surf(unique(mesh.param_x),unique(mesh.param_y),grad_im,'edgecolor','none')
    view([0,0,1])
    title('\nabla m')
    set(gca,'ydir','reverse')
    colorbar
    axis image
    
    
    figure(5)
    surf(unique(mesh.param_x),unique(mesh.param_y),TGV_t1_im_2,'edgecolor','none')
    view([0,0,1])
    title('TGV t1 \lambda_2 = 2')
    set(gca,'ydir','reverse')
    colorbar
    axis image
    
end

figure(4)
surf(unique(mesh.param_x),unique(mesh.param_y),reshape((mesh.res_param1),xdim,[])','edgecolor','none')
view([0,0,1])
title('Model')
set(gca,'ydir','reverse')
colorbar
axis image

if fwd_flag == 1
    % Calculate forward problem too - ideally should save data for inversion
    fem = forward_solver(input,mesh);
    
    fid = fopen([filepath,filename,'.d'],'w');
    fprintf(fid,'%g %g %g %g %g %g %g %g %g\n',[input.ax,input.az,input.bx,input.bz,input.mx,input.mz,input.nx,input.nz,fem.array_model_data]');
    fclose(fid);
end

clear -global
end


% Calculates cx, cy as per smooth matrix4 (with extra options removed)
function [cx, cy] = calc_grad(mesh)
cx = zeros(mesh.num_param);
cy = zeros(mesh.num_param);

tmp_x=unique(mesh.param_x);
tmp_y=unique(mesh.param_y);

% In each row, = 1 where an adjacent parameter is present
for i=1:mesh.num_param
    
    current_x=mesh.param_x(i);
    current_y=mesh.param_y(i);
    ind=find(tmp_x==current_x);
    
    % search all other parameters that have the same y and the x=ind+1
    if ind < (length(tmp_x))
        
        cx(i,i) = -1;
        
        for j=1:1:mesh.num_param
            if mesh.param_y(j)==current_y && mesh.param_x(j)==tmp_x(ind+1)
                
                cx(i,j) = 1;
            end
        end
    end
end

for i=1:1:mesh.num_param
    
    current_x=mesh.param_x(i);
    current_y=mesh.param_y(i);
    ind=find(tmp_y==current_y);
    
    % search all other parameters that have the same x and the y=ind+1
    if ind < (length(tmp_y))
        
        cy(i,i) = -1;
        
        for j=1:1:mesh.num_param
            if mesh.param_y(j)==tmp_y(ind+1) && mesh.param_x(j)==current_x
                
                cy(i,j) = 1;
            end
        end
    end
end



end

function mesh = TGV_sp2(input, mesh)

input.sp2_itr = 100;
% Used in main calculation
mesh.px = zeros(mesh.num_param,1);
mesh.py = mesh.px;
gamma_p = 1e-3;               % cutoff for |x| -> 0
gamma_c = 5e-2;
% display(['\gamma = ',num2str(gamma)])
% Used in performance metrics
p1 = zeros(2*mesh.num_param,1);
rms_p = zeros(1,input.sp2_itr+1);

% main sp2 loop
for i = 2:(input.sp2_itr+1)
    
    % calculater weights
    Rm = diag(1./sqrt((mesh.cx*log10(mesh.res_param1) - mesh.px).^2 + (mesh.cy*log10(mesh.res_param1) - mesh.py).^2 + gamma_c^2));
    Rp = diag(1./sqrt((mesh.cx'*mesh.px).^2 + (mesh.cy'*mesh.py).^2 + 0.5*(mesh.cx'*mesh.py + mesh.cy'*mesh.px).^2 + gamma_p.^2));
    
    
    % Set up matrix equation for px, py and solve
    % a11*px + a12*px = b1;  a21*px + a22*py = b2;
    a11 = Rm + input.tgv_lagrn*(mesh.cx*Rp*mesh.cx' + 0.5*mesh.cy*Rp*mesh.cy');
    a12 = input.tgv_lagrn*0.5*mesh.cy*Rp*mesh.cx';
    a21 = input.tgv_lagrn*0.5*mesh.cx'*Rp*mesh.cy;
    a22 = Rm + input.tgv_lagrn*(mesh.cy*Rp*mesh.cy' + 0.5*mesh.cx*Rp*mesh.cx');
    b1 = Rm*mesh.cx*log10(mesh.res_param1);
    b2 = Rm*mesh.cy*log10(mesh.res_param1);
    A = [a11, a12; a21, a22];
    b = [b1; b2];
    p2 = A\b;
    mesh.px = p2(1:end/2);
    mesh.py = p2(end/2+1:end);
    
    % Store results
    px(:,i) = mesh.px;
    py(:,i) = mesh.py;
    
    
    
    % root mean square difference between last and current
    % solution normalise by rms p value
    rms_p(i) = sqrt(mean((p2 - p1).^2))/(0.5*sqrt(mean(p1.^2 + p2.^2)));
    %         disp(['rms dp = ',num2str(p_rms(i))])
    %                 rms_change(i) = abs((p_rms(i) - p_rms(i-1)))./(p_rms(i) + p_rms(i-1)); % percentage rms change per itr
    p1 = p2;
    
    if ((i > 4) && (rms_p(i) < 0.01)) %  stop if change smaller than 1%
        sp2_i = i;  % store number of iterations needed
        break
    elseif ((rms_p(i) > rms_p(i-1)) && rms_p(i-1) < 0.4 && i > 2) % stop if diverging and
        sp2_i = i;  % store number of iterations needed
        mesh.px = px(:,i-1);
        mesh.py = py(:,i-1);
        break
    end
    
    
    if i == (input.sp2_itr+1)
        sp2_i = i;
    end
end

figure(10)
plot(rms_p)
title('rms change in p')
xlabel('sp 2 iteration')
ylabel('rms dp')

end
