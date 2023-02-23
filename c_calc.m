% calculates differences in sparse form
function [cx, cy] = c_calc(tmp_x, tmp_y)

    % Sparse matrix version: uses known structure of parameter ordering
    % (row major structure:
    % [m(x=1,y=1); m(2,1); ...; m(end-1,end); m(end,end)]
    % number of rows/columns
    
    xL = length(tmp_x);
    yL = length(tmp_y);
    
    % Cx first
    % calculate cell spacings
    dcx = tmp_x(2:end) - tmp_x(1:end-1);
    dcy = tmp_y(2:end) - tmp_y(1:end-1);
    
    cx_i = [-ones(xL - 1, 1)./dcx; 0]; % row of elements at x(i), zero for x(xL) for neuman bc
    cx_i = repmat(cx_i, yL, 1); % repeats for every row in model 
    % This creates a sparse matrix of the (i, i) points
    % sparse(i, j, v, xdim, ydim): (i,j) are indices, v is the element
    % value.
    cx_i = sparse(1:xL*yL, 1:xL*yL, cx_i, xL*yL, xL*yL);
    
    cx_ii = [ones(xL - 1, 1)./dcx; 0]; % elements at x(i + 1), zero row at end for bc 
    cx_ii = repmat(cx_ii, yL, 1);
    cx_ii = sparse( (1:(xL*yL - 1)), (1:(xL*yL) - 1) + 1, cx_ii(1:end-1), xL*yL, xL*yL); % note : last row is zeros for bc
    % combine elements into one matrix
    cx = cx_i + cx_ii;
    
    % Now Cy. This time the bc zero rows all appear at the end.
    cy_i = zeros(xL*yL,1);
    for j = 1:yL-1
        cy_i((1 + xL*(j-1)):(xL*j)) = -1./dcy(j);
    end
    cy_ii = -cy_i; % elements at y =(i+1)
    cy_i = sparse(1:xL*yL, 1:xL*yL, cy_i, xL*yL, xL*yL);
    cy_ii = sparse( (1:(xL*yL - xL)), (1:(xL*yL) - xL) + xL, cy_ii(1:end-xL), xL*yL, xL*yL); % last rows are zeros (bc) so are not included explicitly
    
    cy = cy_i + cy_ii;


end

