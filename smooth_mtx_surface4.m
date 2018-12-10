function mesh=smooth_mtx_surface4(input,mesh)

% factor by which surface layer damping is increased
surface_damp_strength = 15;

c=zeros(mesh.num_param,mesh.num_param);
cx=zeros(mesh.num_param,mesh.num_param);
cy=zeros(mesh.num_param,mesh.num_param);
dcx = 1;
dcy = 1;



% tmp_x=union(mesh.tmp_param(:,1),mesh.tmp_param(:,1));
% tmp_y=union(mesh.tmp_param(:,2),mesh.tmp_param(:,2));


tmp_x = unique(mesh.param_x);
tmp_y = unique(mesh.param_y);

% Add dummy rows/columns at the ends for length calculations
tmp_x = [2*tmp_x(1) - tmp_x(2); tmp_x; 2*tmp_x(end) - tmp_x(end-1)];
tmp_y = [2*tmp_y(1) - tmp_y(2); tmp_y; 2*tmp_y(end) - tmp_y(end-1)];


% In each row, = 1 where an adjacent parameter is present
for i=1:mesh.num_param
    
    current_x=mesh.param_x(i);
    current_y=mesh.param_y(i);
    indX = find(tmp_x == current_x);
    indY = find(tmp_y == current_y);
    % search all other parameters that have the same y and the x=ind+1
    if indX < (length(tmp_x)-1) 
        
        if input.fd_weight == 1
%             dcx = (tmp_y(indY+1) - tmp_y(indY - 1))/2; 
              dcx = (tmp_x(indX + 1) - tmp_x(indX))/2;
        end
        
       if input.damp_surface_flag == 1 && abs(current_y) == min(abs(mesh.param_y)) 
            cx(i,i) = -surface_damp_strength.*1./dcx;
       else
            cx(i,i) = -1./dcx;         
       end
        
        for j=1:mesh.num_param
            
          
            if (mesh.param_y(j)==current_y) && (mesh.param_x(j)==tmp_x(indX+1))
               
                if input.damp_surface_flag == 1 && abs(current_y) == min(abs(mesh.param_y)) 
                    cx(i,j) = surface_damp_strength.*1./dcx; 
                else
                    cx(i,j) = 1./dcx;
                end
            end
        end       
       
    end
    
end

        

for i=1:mesh.num_param
    
    current_x=mesh.param_x(i);
    current_y=mesh.param_y(i);
    indX = find(tmp_x == current_x);
    indY = find(tmp_y == current_y);
    % search all other parameters that have the same x and the y=ind+1
    if indY < (length(tmp_y)-1)
        
        if input.fd_weight == 1
%             dcy = (tmp_x(indX + 1) - tmp_x(indX - 1))/2;
            dcy = tmp_y(indY+1) - tmp_y(indY);
        end        
        
        if input.damp_surface_flag == 1 && abs(current_y) == min(abs(mesh.param_y))
            cy(i,i) = 1./dcx;  
%             cy(i,i) = -surface_damp_strength.*1./dcx;  
        else
            cy(i,i) = -1/dcy;
        end
        
        for j=1:mesh.num_param            
            if mesh.param_y(j)==tmp_y(indY+1) && mesh.param_x(j)==current_x
                if input.damp_surface_flag == 1 && abs(current_y) == min(abs(mesh.param_y))
                    cy(i,j) = 1./dcx;
%                     cy(i,i) = -surface_damp_strength.*1./dcx;  
                else
                    cy(i,j) = 1/dcy;
                end
            end
        end 
    end
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

save('C','cx','cy')
        

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
for i=1:mesh.num_param
    
    for j=1:mesh.num_param
        
        if (c(i,j)~=0)
            mesh.S(i,j)=1;
        else
            mesh.S(i,j)=0;
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
        
