% [param_x, param_y] = param_generator(tmp_x, tmp_y)
% given the unique x and y grid values, tmp_x and tmp_y, constructs param_x
% and param_y mesh vectors for that mesh. 
%
% % [param_x, param_y, param] = param_generator(tmp_x, tmp_y, param, [ghost_flag_x, ghost_flag_y])
% also extends parameter mesh onto the ghost points. 
% ghost flag_(x/y) = 1 : takes parameters from 0:end-1 in (x/y)
% ghost flag_(x/y) = 2 : takes parameters from 2:end in (x/y)
% ghost_flag_(x/y) = 3 : takes parameters from 0:end in (x/y)
% ghost_flag_(x/y) = 4 : takes parameters from 1:end+1 in (x/y)
% Only have ghost_flag_x < 0 || ghost_flag_y < 0 with current
% implementation
function [param_x, param_y, varargout] = param_generator(tmp_x, tmp_y, varargin)

xL = length(tmp_x);
yL = length(tmp_y);

param_x = repmat(tmp_x, yL, 1);

param_y = zeros(size(param_x));
for i = 1:yL
    param_y( (1 + (i-1)*xL):(i*xL) ) = repmat(tmp_y(i),length(xL),1);
end

if nargin > 3
%     try
        param = varargin{1};
        ghost_flag_x = varargin{2}(1);
        ghost_flag_y = varargin{2}(2);

        param_new = zeros(xL*yL,1);
        param_tmp = param;
        
        % Need to work out how long x is in param
        if ghost_flag_x == 1 && ghost_flag_x == 2
            xL = xL-1;
        elseif ghost_flag_x == 1 || ghost_flag_x == 2
              xL1 = xL-1;
        elseif ghost_flag_x == 3 || ghost_flag_x == 4
            xL1 = xL + mod(length(param)/yL, xL);
        else
            xL1 = xL;
        end

        
        % extend in y first      
        if ghost_flag_y == 3 
            param_tmp(end - xL + 1:end) = []; % first remove y(end)
            param_new = [param(1:xL); param_tmp]; % add y(0)
        elseif ghost_flag_y == 4
            param_tmp(1:xL) = []; % first remove y(0)
            param_new = [param_tmp; param(end - xL + 1:end)]; % adds y(end)
        elseif ghost_flag_y == 1
            param_new = [param(1:xL); param];     
        elseif ghost_flag_y == 2
            param_new = [param; param(end-xL+1:end)];
        end
        
        % then extend in x - need to populate param_new.
        if ghost_flag_x > 0
            ind_ghost = false(length(param_new),1);
        end
        if ghost_flag_x == 3
            ind_ghost(1:xL:xL*yL) = 1;
            param_new(ind_ghost) = param_tmp(1:xL1:end); % adds new bc values
            param_tmp(xL:xL:end) = [];
            param_new(~ind_ghost) = param_tmp; % fills in old values
        elseif ghost_flag_x == 4
            ind_ghost(xL:xL:xL*yL) = 1;
            param_new(ind_ghost) = param_tmp(xL1:xL1:end);
            param_tmp(1:xL:end) = [];
            param_new(~ind_ghost) = param_tmp;
        elseif ghost_flag_x == 1
            ind_ghost(1:xL:xL*yL) = 1;
            param_new(ind_ghost) = param_tmp(1:xL1:end);
            param_new(~ind_ghost) = param_tmp;
        elseif ghost_flag_x == 2
            ind_ghost(xL:xL:xL*yL) = 1;
            param_new(ind_ghost) = param_tmp(xL1:xL1:end);
            param_new(~ind_ghost) = param_tmp;
        end
        
        varargout{1} = param_new; % output new parameter vector with ghost points   
        
%     catch
%         error('check optional inputs to param_generator - are the dimensions right?')
%     end
end
end