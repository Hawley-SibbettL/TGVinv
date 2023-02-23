% [ctc_x, ctc_y, cxRc, cyRc] = ctc_half(mesh, input, px, py, gamma, itr)
% function to calculate spatial terms for the half differences fd scheme
function [ctc, ctpx, ctpy] = ctc_half(mesh, input, gamma, itr, bc_flag)

xL = length(unique(mesh.param_x));
yL = length(unique(mesh.param_y));
res = log10(mesh.res_param1);

% calculate the individual difference terms used in the outputs

% preallocate
dxf = 0*mesh.param_x;
dxb = dxf; dyb = dxf; dyf = dxf; dy_xp = dxf; dy_xm = dxf; dx_yp = dxf; dx_ym = dxf;

% in x
for j = 1:yL
    rowInd = (j-1)*xL + 1 : j*xL;
    xj = res(rowInd); % row of parameters with y = tmp_y(j)
    
    % now take the x rows at i + 1, and i - 1, applying bc if necessary
    if j ~= 1 && j~= yL
        xjp = res(rowInd + xL); % x at j+1
        xjm = res(rowInd - xL);  % x at j-1
    elseif j == 1
        xjp = res(rowInd + xL); % x at j+1
        
        if bc_flag == 0
            xjm = res(end-xL+1:end);  % x at j = yL
        elseif bc_flag == 1
            
        elseif bc_flag == 2
            
        end
        
    elseif j == yL
        xjm = res(rowInd - xL);  % x at j-1
        
        if bc_flag == 0
            xjp = res(1:xL); % x at j = yL
        elseif bc_flag == 1
            
        elseif bc_flag == 2
            
        end
    end


    % periodic bc by default here
    dxf(rowInd) = circshift(xj, -1) - xj; % forward difference in x
    dxb(rowInd) = xj - circshift(xj, -1); % backward difference in x
    dyf(rowInd) = xjp - xj; % forward difference in x
    dyb(rowInd) = xj - xjm; % backward difference in y
    
    % periodic bc by default here
    dy_xp(rowInd) = ( xjp + circshift(xjp, -1) - xjm - circshift(xjm, -1) )/2; % y difference at i + 1/2
    dy_xm(rowInd) = ( xjp + circshift(xjp, 1) - xjm - circshift(xjm, -1))/2; % y difference at i - 1/2
    dx_yp(rowInd) = ( circshift(xj, -1) + circshift(xjp, -1) - circshift(xj, 1) - circshift(xjp, 1) )/2; % x difference at j + 1/2
    dx_ym(rowInd) = ( circshift(xj, -1) + circshift(xjm, -1) - circshift(xj, 1) - circshift(xjm, 1) )/2; % x difference at j - 1/2
    
    % for other bc
    if bc_flag == 1 % dm/dx = 0, dp/dx = 0
        
    elseif bc_flag == 2 % d^2m/dx^2 = p, dp/dx = o
        
    end
    
end

    % combine into ctc calculations

    % ctc_ij (if _11 is the main diagonal)
    ctc_11 = - 1./sqrt( dxf.^2 + dy_xp.^2 + gamma.^2 ) - 1./sqrt( dxb.^2 + dy_xm.^2 + gamma.^2 ) - 1./( dyf.^2 + dx_yp.^2 + gamma.^2 ) - 1./(dyb.^2 + dx_ym.^2 + gamma.^2 ); % m(i,j)
    ctc_21 = 1./sqrt( dxf.^2 + dy_xp.^2 + gamma.^2 ); % m(i+1, j)
    ctc_01 = 1./sqrt( dxb.^2 + dy_xm.^2 + gamma.^2 );% m(i-1, j)
    ctc_12 = 1./sqrt( dyf.^2 + dx_yp.^2 + gamma.^2 );% m(i, j+1)
    ctc_10 = 1./sqrt( dyb.^2 + dx_ym.^2 + gamma.^2 ); % m(i, j-1)
    
    
    % Note the different indexing for ctpx, ctpy, as they act on the px, py
    % vectors, which are longer due to being on the edges
    
    % ctpx_ij (if _11 is main diagonal is i-1/2, j)
    ctpx_11 = 1./sqrt( dxf.^2 + dy_xp.^2 + gamma.^2 ); % p(i+1/2,j)
    ctpx_01 = 1./sqrt( dxb.^2 + dy_xm.^2 + gamma.^2 ); % p(i-1/2,j)
    
    % ctpy (if _11 is main diagonal is i, j-1/2)
    ctpy_11 = 1./sqrt( dyf.^2 + dx_yp.^2 + gamma.^2 ); % p(i, j+1/2)
    ctpy_10 = 1./sqrt( dyb.^2 + dx_ym.^2 + gamma.^2 ); % p(i, j-1/2)
    

    
    ctc = spdiags([ctc_11, ctc_21, ctc_01, ctc_12, ctc_10],[0, -1, 1, -xL, xL],xL*yL, xL*yL)';
    ctpx = spdiags([ctpx_11, ctpx_01], [0, 1], (xL+1)*yL, xL*yL)';
    ctpy = spdiags([ctpy_11, ctpy_10], [0, xL], xL*(yL+1), xL*yL)';



    
    % enforce bc
    if bc_flag == 0
        % wraparound y cells
        ctc(sub2ind( size(ctc), xL*yL-xL+1:xL*yL, 1:xL )) = ctc_12(end-xL+1:end); 
        ctc(sub2ind( size(ctc), 1:xL, xL*yL-xL+1:xL*yL )) = ctc_10(end-xL+1:end);
       
        
        % wraparound x cells - due ot the structure of the matrix, only the
        % end points are not included automatically
        ctc(end, 1) = ctc_21(end);
        ctc(1, end) = ctc_01(1);
%         ctcpx(end, 1) = ctcpx_11;
        
    elseif bc_flag == 1
        
    elseif bc_flag == 2
        
    end
    


end