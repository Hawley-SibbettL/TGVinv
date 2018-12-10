function mesh=smooth_mtx_surface4(input,mesh)

w_edge = 0.5; % Weight of repeated difference value at edge cells


c=zeros(mesh.num_param,mesh.num_param);
cx=zeros(mesh.num_param,mesh.num_param);
cy=zeros(mesh.num_param,mesh.num_param);
cx_O1 = cx; % NEW: 1st order smoothingmatrix
cy_O1 = cx; % NEW: 1st order smoothingmatrix


% tmp_x=union(mesh.tmp_param(:,1),mesh.tmp_param(:,1));
% tmp_y=union(mesh.tmp_param(:,2),mesh.tmp_param(:,2));


tmp_x=unique(mesh.param_x);
tmp_y=unique(mesh.param_y);


% In each row, = 1 where an adjacent parameter is present
for i=1:mesh.num_param
    
    current_x=mesh.param_x(i);
    current_y=mesh.param_y(i);
    ind=find(tmp_x==current_x);
    % search all other parameters that have the same y and the x=ind+1
    if ind < (length(tmp_x)) 
        cx_O1(i,i) = -1;        % 1st Order
        
        for j=1:mesh.num_param
            
            if mesh.param_y(j)==current_y && mesh.param_x(j)==tmp_x(ind+1)
                cx(i,j)=1;
                cx_O1(i,j) = 1;  % 1st Order
                %                cx(i,j)=sqrt(      (mesh.tmp_param(j,6)-mesh.tmp_param(j,5))/ ( mesh.tmp_param(j,1)-mesh.tmp_param(i,1)));
            end
        end
        if ind == (length(tmp_x) - 1)
            cx_O1(i,:) = w_edge.*cx_O1(i,:);
        end
        % Need something for end cells: same a s prev row (half magnitude)
    else
        
        for j=1:mesh.num_param
            if mesh.param_y(j)==current_y && mesh.param_x(j)==tmp_x(length(tmp_x) - 1)
                cx_O1(i,j) = -w_edge;
            elseif mesh.param_y(j)==current_y && mesh.param_x(j)==tmp_x(length(tmp_x))
                cx_O1(i,j) = w_edge;  % 1st Order
            end
        end
        
    end
    
end

for i=1:mesh.num_param
   cx(i,i)=-sum(cx(i,:));    % 2nd order
end


ctc1=cx'*cx;
        

for i=1:mesh.num_param
    
    current_x=mesh.param_x(i);
    current_y=mesh.param_y(i);
    ind=find(tmp_y==current_y);
    % search all other parameters that have the same y and the x=ind+1
    if ind < (length(tmp_y))
        cy_O1(i,i) = -1;
        
        for j=1:mesh.num_param            
            if mesh.param_y(j)==tmp_y(ind+1) && mesh.param_x(j)==current_x
                cy(i,j)=1;
                cy_O1(i,j) = 1;
                %                cy(i,j)=sqrt(      (mesh.tmp_param(j,4)-mesh.tmp_param(j,3))/ ( mesh.param_y(j)-mesh.param_y(i)) );
            end
        end
        
        if ind == (length(tmp_y) - 1)
            cy_O1(i,:) = w_edge.*cy_O1(i,:);
        end
    else
        
        for j=1:mesh.num_param
            if mesh.param_y(j)==tmp_y(length(tmp_y)-1) && mesh.param_x(j)==current_x
                cy_O1(i,j) = -w_edge;
            elseif mesh.param_y(j)==tmp_y(length(tmp_y)) && mesh.param_x(j)==current_x
                cy_O1(i,j) = w_edge;
            end
        end            
    end
    
end

for i=1:mesh.num_param
   cy(i,i)=-sum(cy(i,:));    % 2nd order central difference 
end
            
            
 ctc2=cy'*cy;  
 
 mesh.c_O1 = cx_O1 + cy_O1;
 mesh.cx_O1 = cx_O1;
 mesh.cy_O1 = cy_O1;
 mesh.ctc_O1 = mesh.c_O1'.*mesh.c_O1;

        
mesh.ctc=ctc1+ctc2; 
% mesh.c = cx + cy;
mesh.cx = cx;
mesh.cy = cy;
        

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
        
