function mesh=smooth_mtx_surface4(input,mesh)

if input.ghost_flag == 1
    [num_param, ~, param_x, param_y, ~] = neumann_bc(mesh.num_param, mesh.param_x, mesh.param_y, mesh.res_param1, [0, 2]);
else
    num_param = mesh.num_param;
    param_x = mesh.param_x;
    param_y = mesh.param_y;
end

% factor by which surface layer damping is increased
surface_damp_strength = 15;

c = zeros(num_param,num_param);
cx = zeros(num_param,num_param);
cy = zeros(num_param,num_param);
dcx = 1;
dcy = 1;


tmp_x = unique(param_x);
tmp_y = unique(param_y);

% change x cell spacing to reflect true mesh - no longer required
mesh.ye = mesh.depth_n;
mesh.xe = (tmp_x(1:end-1) + tmp_x(2:end))/2; % internal edges
dx = tmp_x(3) - tmp_x(2);
% tmp_x = [tmp_x(2) - 1.5*dx; tmp_x(2:end-1); tmp_x(end-1) + 1.5*dx];
% mesh.xe = [tmp_x(1) - dx; mesh.xe; tmp_x(end) + dx]; % add external edges, accounting for double-spaced end cells
% [param_x, ~] = param_generator(tmp_x, tmp_y);
mesh.xe = [mesh.xe(1) - dx; mesh.xe; mesh.xe(end) + dx]; % add extra rows if edge cells are normal sized for diff_type = 3

if input.diff_type == 1 || input.diff_type == 3 % forward/backward differences 
    % original, non sparse formulation
    
    %    % In each row, = 1 where an adjacent parameter is present
    %     for i=1:num_param
    %
    %         current_x=param_x(i);
    %         current_y=param_y(i);
    %         indX = find(tmp_x == current_x);
    %         indY = find(tmp_y == current_y);
    %         % search all other parameters that have the same y and the x=ind+1
    %         if indX < (length(tmp_x))
    %
    %             if input.fd_weight == 1
    %                 dcx = (tmp_x(indX + 1) - tmp_x(indX));
    %             end
    %
    %             if input.damp_surface_flag == 1 && abs(current_y) == min(abs(param_y))
    %                 cx(i,i) = -surface_damp_strength.*1./dcx;
    %             else
    %                 cx(i,i) = -1./dcx;
    %             end
    %
    %             for j=1:num_param
    %
    %
    %                 if (param_y(j)==current_y) && (param_x(j)==tmp_x(indX+1))
    %
    %                     if input.damp_surface_flag == 1 && abs(current_y) == min(abs(param_y))
    %                         cx(i,j) = surface_damp_strength.*1./dcx;
    %                     else
    %                         cx(i,j) = 1./dcx;
    %                         break
    %                     end
    %                 end
    %             end
    %
    %         end
    %
    %     end
    %
    %
    %
    %     for i=1:num_param
    %
    %         current_x=param_x(i);
    %         current_y=param_y(i);
    %         indX = find(tmp_x == current_x);
    %         indY = find(tmp_y == current_y);
    %         % search all other parameters that have the same x and the y=ind+1
    %         if indY < (length(tmp_y))
    %
    %             if input.fd_weight == 1
    %                 dcy = (tmp_y(indY+1) - tmp_y(indY));
    %             end
    %
    %             if input.damp_surface_flag == 1 && abs(current_y) == min(abs(param_y))
    %                 cy(i,i) = 1./dcy;
    %                 %             cy(i,i) = -surface_damp_strength.*1./dcx;
    %             else
    %                 cy(i,i) = -1/dcy;
    %             end
    %
    %             for j=1:num_param
    %                 if param_y(j)==tmp_y(indY+1) && param_x(j)==current_x
    %                     if input.damp_surface_flag == 1 && abs(current_y) == min(abs(param_y))
    %                         cy(i,j) = 1./dcy;
    %                         %                     cy(i,i) = -surface_damp_strength.*1./dcx;
    %                     else
    %                         cy(i,j) = 1/dcy;
    %                         break
    %                     end
    %                 end
    %             end
    %         end
    %     end
    
    
    [cx, cy] = c_calc(tmp_x, tmp_y);
    
    if input.periodic_bc == 1
            [cx, cy] = c_calc_periodic(tmp_x, tmp_y, input.diff_weight_flag);
    end
    
    
    % now calculate cx, cy matrices interpolated onto the corners and edges
    % these include those on the boundary of the image domain, which will
    % have their value determined by homogeneous neumann boundary conditions
    if input.diff_type == 3
        
        % In this case, 'true' underlying operator is C*, where
        %  (Cx*m)_{i}= m_{i+1/2} - m_{i-1/2}/(x_{i+1/2} - x_{i-1/2})
        % A factor of 2 is needed so that Cx*Cx*m = Cx^TCxm, due to the
        % difference in dx
        cx = 2.*cx;
        cy = 2.*cy;
        
        [mesh.cxe, mesh.cye] = c_calc(mesh.xe, mesh.ye);
        
    elseif input.diff_type == 4
       
        
        % calculate cell widths, TAKING INTO ACCOUNT ENLARGED X EDGE CELLS
        tmp_ye = mesh.depth_n; % y edges
        dx = tmp_x(3) - tmp_x(2);
        tmp_x = [tmp_x(2) - 1.5*dx; tmp_x(2:end-1); tmp_x(end-1) + 1.5*dx];
        tmp_xe = (tmp_x(2:end) + tmp_x(1:end-1))./2; % x edges
        % finite difference is weighted by cell widths
        dx = tmp_xe(2:end) - tmp_xe(1:end-1);
        dy = tmp_ye(2:end) - tmp_ye(1:end-1);
        xL = length(tmp_x);
        yL = length(tmp_y);
        
        % recalculate mesh.param_x accurately
        [mesh.param_x, ~] = param_generator(tmp_x, tmp_y);
        
        % now calculate cx values
        cx_diag = [-1./dx(1), - ones(xL-2, 1)./dx(2:end-1) - ones(xL - 2, 1)./dx(2:end-1), -1./dx(end)]; % along diagonal

        
        cx_diag = repmat(cx_i, yL, 1); % repeats for every row in model 
        % This creates a sparse matrix of the (i, i) points
        % sparse(i, j, v, xdim, ydim): (i,j) are indices, v is the element
        % value.
        cx_i = sparse(1:xL*yL, 1:xL*yL, cx_diag, xL*yL, xL*yL);
    



        
        
    end
    
elseif input.diff_type == 2 % central differences

    xL = length(tmp_x);
    yL = length(tmp_y);
    
    % Cx first
    % calculate difference spacings in x and y
    dcx = [2*(tmp_x(2) - tmp_x(1)) ;tmp_x(3:end) - tmp_x(1:end-2); 2*(tmp_x(end) - tmp_x(end-1))];
    dcy = [2*(tmp_y(2) - tmp_y(1)) ;tmp_y(3:end) - tmp_y(1:end-2); 2*(tmp_y(end) - tmp_y(end-1))];
    
    % m0 = m2
%     % create diagonals of cx matrix
%     cx0 = repmat(-ones(xL, 1)./dcx, yL, 1); % diagonal at (i - 1)
%     cx2 = repmat(ones(xL, 1)./dcx, yL, 1); % diagonal at (i + 1)
%     % remove points outside boundaries
%     cx0([1:xL:end, xL:xL:end]) = 0;
%     cx2([1:xL:end, xL:xL:end]) = 0;
%     % add bc points to main diagonal
% 
%     % create cy diagonals
%     cy = zeros(xL*yL, 1);
%     for i = 1:yL
%         cy((i-1)*xL+1:i*xL) = 1./dcy(i);
%     end
%     % initialise cy diagonals
%     cy0 = -cy;  cy2 = cy;
%     % apply boundary conditions
%     cy0([1:xL, end-xL+1:end]) = 0;
%     cy2([1:xL, end-xL+1:end]) = 0;
% 
%     cx = spdiags([cx2, cx0], [-1, 1], xL*yL, xL*yL)';
%     cy = spdiags([cy2, cy0], [-xL, xL], xL*yL, xL*yL)';

    % m0 = m1
    cx0 = repmat(-ones(xL, 1)./dcx, yL, 1); % diagonal at (i - 1)
    cx2 = repmat(ones(xL, 1)./dcx, yL, 1); % diagonal at (i + 1)
    % remove points outside boundaries
    cx0([1:xL:end]) = 0;
    cx2([xL:xL:end]) = 0;
    % add bc points to main diagonal
    cx1 = zeros(xL*yL, 1);
    cx1(1:xL:end) = -1./dcx(1);
    cx1(xL:xL:end) = 1./dcx(end); 
%     cx1(xL:xL:end) = -1./dcx(end);
%     cx0(xL:xL:end) = 1./dcx(end); 
%     
    % create cy diagonals
    cy = zeros(xL*yL, 1);
    for i = 1:yL
        cy((i-1)*xL+1:i*xL) = 1./dcy(i);
    end
    % initialise cy diagonals
    cy0 = -cy; cy1 = cy; cy2 = cy;
    % apply boundary conditions
    cy0([1:xL]) = 0;
    cy2([end-xL+1:end]) = 0;
    cy1(1:xL) = - 1./dcy(1);
    cy1(xL+1:end-xL) = 0;
    cy1(end-xL+1:end) = 1./dcy(end);
%      cy1(end-xL+1:end) = -1./dcy(end); % backward difference
%      cy0(end-xL+1:end) = 1./dcy(end);  % backward difference
    
    cx = spdiags([cx2, cx1, cx0], [-1, 0, 1], xL*yL, xL*yL)';
    cy = spdiags([cy2, cy1, cy0], [-xL, 0, xL], xL*yL, xL*yL)';
    
    cx0(2:xL:end) = -2/dcx(1);
    cx2(xL-1:xL:end) = 2/dcx(1);
    cy0(xL+1:2*xL) = -2/dcy(1);
    cy2(end-2*xL+1:end-xL) = 2/dcy(1);
    ctc0 = 4*ones(xL*yL, 1)./dcx(1); % only correct for square cells
    ctc = spdiags([-cy2, -cx2, ctc0, cx0, cy0], [-xL, -1, 0, 1, xL], xL*yL, xL*yL);
    
%     for i=1:num_param
%         
%         current_x=param_x(i);
%         current_y=param_y(i);
%         indX = find(tmp_x == current_x);
%         indY = find(tmp_y == current_y);
%         % for all domain points
%         if indX < (length(tmp_x)) && indX > 1
%             
%             if input.fd_weight == 1
%                 dcx = (tmp_x(indX + 1) - tmp_x(indX - 1));
%             end
%             
%             for j=1:num_param
%                 
%                 if (param_y(j)==current_y)
%                     
%                     if (param_x(j)==tmp_x(indX-1))
%                         if input.damp_surface_flag == 1 && abs(current_y) == min(abs(param_y))
%                             cx(i,i) = -surface_damp_strength.*1./dcx;
%                         else
%                             cx(i,i) = -1./dcx;
%                         end
%                     elseif(param_x(j)==tmp_x(indX+1))
%                         if input.damp_surface_flag == 1 && abs(current_y) == min(abs(param_y))
%                             cx(i,j) = surface_damp_strength.*1./dcx;
%                         else
%                             cx(i,j) = 1./dcx;
%                             break
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     
%     for i=1:num_param
%         
%         current_x=param_x(i);
%         current_y=param_y(i);
%         indX = find(tmp_x == current_x);
%         indY = find(tmp_y == current_y);
%         
%         
%         
%         % search all other parameters that have the same x and the y=ind+1
%         if indY < (length(tmp_y)) && indY > 1
%             
%             if input.fd_weight == 1
%                 dcy = (tmp_y(indY+1) - tmp_y(indY-1));
%             end
%             
%             for j=1:num_param
%                 
%                 if param_x(j)==current_x
%                     
%                     if param_y(j)==tmp_y(indY-1)
%                         if input.damp_surface_flag == 1 && abs(current_y) == min(abs(param_y))
%                             cy(i,i) = 1./dcy;
%                             %             cy(i,i) = -surface_damp_strength.*1./dcx;
%                         else
%                             cy(i,i) = -1/dcy;
%                         end
%                     elseif param_y(j)==tmp_y(indY+1)
%                         if input.damp_surface_flag == 1 && abs(current_y) == min(abs(param_y))
%                             cy(i,j) = 1./dcy;
%                             %                     cy(i,i) = -surface_damp_strength.*1./dcx;
%                         else
%                             cy(i,j) = 1/dcy;
%                             break
%                         end
%                     end
%                 end
%             end
%         end
%     end
    
    
end

% Normalise cx, cy so that the mean element value across both matrices = 1
% (otherwise this will effect the relative weighting of terms)
% c_norm = sum(sum( [abs(cx), abs(cy)] ))/(4*length(cx));
% cx = cx/mean(cx(cx~=0));
% cy = cy/mean(cy(cy~=0));




ctc1=cx'*cx;
ctc2=cy'*cy;

mesh.ctc=ctc1+ctc2;
mesh.cx = cx;
mesh.cy = cy;

% save('C','cx','cy')


%ctc=eye(num_param,num_param);

%/* adjust weighting at the edges */

% if smooth_type==11
%     for i=1:num_param
%         s1=0;
%         for j=1:num_param
%             if(c(i,j)~=0) s1=s1+1; end
%         end
%         if(s1<=2) c(i,i)=1.15*c(i,j); end
%     end
%
% end




% This matrix has no meaning. I just calculate so I can now the elemeents
% that do not have zero. After this I can calcualte the S matrix (ACB, Kim)
c=cx+cy;
c=mesh.ctc;
% /*Here i calculate the S matrix, which is one when C is non zero and zero otherwise*/
for i=1:num_param
    
    for j=1:num_param
        
        if (c(i,j)~=0)
            S(i,j)=1;
        else
            S(i,j)=0;
        end
    end
end


% Here create time related constarain matrix M
if input.time_lapse_flag==1
    mesh.M=eye(input.num_files*mesh.num_param,input.num_files*mesh.num_param);
    for i=1:input.num_files*mesh.num_param
        for j=1:input.num_files*mesh.num_param
            
            if j-i==mesh.num_param
                mesh.M(i,j)=-1;
            end
            
            if i==j && i>(input.num_files-1)*mesh.num_param && j>(input.num_files-1)*mesh.num_param
                mesh.M(i,j)=0;
            end
        end
    end
    
end


end

