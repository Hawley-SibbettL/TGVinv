clear all
close all
% creates a psuedosection and interpolates onto res_param mesh for use as
% initialisation. ASSUMES DIPOLE-DIPOLE DATA ONLY

load('pseudo_init_practice_workspace')

% Calculate dipole centres, a and n values
% na is distance between inner potential and current electrodes
a = input.ax - input.bx;
n = ( input.mx - input.ax ) ./ a;

% z over a ratio for n = 1:8
z_a_ratio = [.416, .697, .962, 1.22, 1.476, 1.73, 1.983, 2.236]';

pseudo_x = ( input.mx + input.ax ) / 2; % dipole centre points
pseudo_z = z_a_ratio( n ) .* a; % dipole depth points

% perform interpolation (requires post 2013 matlab version...)
% Use log space for interpolation to be consistant with inversion
% interpFun = scatteredInterpolant( pseudo_x, pseudo_z, log10( input.real_data ) );% , 'ExtrapolationMethod', 'nearest');
interpFun = scatteredInterpolant( pseudo_x, pseudo_z, log10( input.real_data ), 'linear', 'nearest');%, 'ExtrapolationMethod', 'none');
pseudosection = 10.^interpFun( mesh.param_x, mesh.param_y );

[xgrid, ygrid] = meshgrid( unique(mesh.param_x), unique(mesh.param_y) );

figure
surf( xgrid, ygrid, reshape( log10 (pseudosection ), length( unique( mesh.param_x ) ), [])' )
axis equal
view( [ 0, 0 , 1 ] )
set(gca, 'ydir', 'reverse')