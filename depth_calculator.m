% Calculates z, a vector of z values giving the depth to each layer, up to
% a specified maximum depth
% max_depth - maximum depth of mesh (calculated from probe spacing by
% probe spacing - electrode spacing on surface
% refine mesh flag - sets cell width to half the electrode spacing.


function [z] = depth_calculator(max_depth,probe_spacing,refine_mesh_flag)

z1 = [0; 0]; % initialise surface as z = 0 (must be a column vector)
i = 2; % layer index
dz_min = 1; % 1 - constant cell depth. 1.1 - standard 10% increment
custom_depth = max_depth*1.3; % set custom maximum depth, Useful in the case that ip4di 
                   % underestimates model depth

if refine_mesh_flag == 1
    probe_spacing = probe_spacing/2;
else
    % HERE PUT CUSTOM DETAILS
%     % For z refinement only
%         probe_spacing = probe_spacing/2;

    max_depth = custom_depth;

end  
    
    
while z1(i-1) < max_depth
    z1(i) = z1(i-1) + probe_spacing*(dz_min).^(i-1);  
    i = i + 1;
end
z = z1;


end