function [cx, cy] = p_dif(num_param, param_x, param_y)
% Function calculates the fd matrices for p



% factor by which surface layer damping is increased
surface_damp_strength = 15;

c=zeros(num_param,num_param);
cx=zeros(num_param,num_param);
cy=zeros(num_param,num_param);
dcx = 1;
dcy = 1;


tmp_x = unique(param_x);
tmp_y = unique(param_y);


% In each row, = 1 where an adjacent parameter is present
for i=1:num_param
    
    current_x=param_x(i);
    current_y=param_y(i);
    indX = find(tmp_x == current_x);
    indY = find(tmp_y == current_y);
    % search all other parameters that have the same y and the x=ind+1
    if indX < (length(tmp_x)) 
        
        if input.fd_weight == 1
%             dcx = (tmp_y(indY+1) - tmp_y(indY - 1))/2; 
              dcx = (tmp_x(indX + 1) - tmp_x(indX))/2;
        end
        
       if input.damp_surface_flag == 1 && abs(current_y) == min(abs(param_y)) 
            cx(i,i) = -surface_damp_strength.*1./dcx;
       else
            cx(i,i) = -1./dcx;         
       end
        
        for j=1:num_param
            
          
            if (param_y(j)==current_y) && (param_x(j)==tmp_x(indX+1))
               
                if input.damp_surface_flag == 1 && abs(current_y) == min(abs(param_y)) 
                    cx(i,j) = surface_damp_strength.*1./dcx; 
                else
                    cx(i,j) = 1./dcx;
                    break
                end
            end
        end       
       
    end
    
end

        

for i=1:num_param
    
    current_x=param_x(i);
    current_y=param_y(i);
    indX = find(tmp_x == current_x);
    indY = find(tmp_y == current_y);
    % search all other parameters that have the same x and the y=ind+1
    if indY < (length(tmp_y))
        
        if input.fd_weight == 1
%             dcy = (tmp_x(indX + 1) - tmp_x(indX - 1))/2;
            dcy = tmp_y(indY+1) - tmp_y(indY);
        end        
        
        if input.damp_surface_flag == 1 && abs(current_y) == min(abs(param_y))
            cy(i,i) = 1./dcx;  
%             cy(i,i) = -surface_damp_strength.*1./dcx;  
        else
            cy(i,i) = -1/dcy;
        end
        
        for j=1:num_param            
            if param_y(j)==tmp_y(indY+1) && param_x(j)==current_x
                if input.damp_surface_flag == 1 && abs(current_y) == min(abs(param_y))
                    cy(i,j) = 1./dcx;
%                     cy(i,i) = -surface_damp_strength.*1./dcx;  
                else
                    cy(i,j) = 1/dcy;
                    break
                end
            end
        end 
    end
end
            

