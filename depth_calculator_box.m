% Calculates z, a vector of y values giving the depth to each layer, up to
% a specified maximum depth
function [z] = depth_calculator(max_depth,probe_spacing,refine_mesh_flag)

% Boxford res2Dinv values
z = [.1386,.2912,.4589,.64355,.84657,1.0698,1.3155,1.5857,1.8830,2.2099,2.569,2.965,3.4004,3.8791,4.4057,4.9850,5.6222,6.323,7.0940,7.942,8.875,9.90118,11.02997];
% % z = [.2642,.64355,1.0698,1.5857,2.2,2.569,2.965,3.4004,3.8791,4.4057,4.9850,5.6222,6.323,7.0940,7.942,8.875,9.90118,11.02997];
z = [0, (z(1:end-1) + z(2:end))/2];

end