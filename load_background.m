
% Script loads in all the bckground models for use in analysing the results
% from the synthetic models. All for the dipole-dipole 35 electrode array
% with 3m spacing

x_param =  [1.50000000000000;4.50000000000000;7.50000000000000;10.5000000000000;13.5000000000000;16.5000000000000;19.5000000000000;22.5000000000000;25.5000000000000;28.5000000000000;31.5000000000000;34.5000000000000;37.5000000000000;40.5000000000000;43.5000000000000;46.5000000000000;49.5000000000000;52.5000000000000;55.5000000000000;58.5000000000000;61.5000000000000;64.5000000000000;67.5000000000000;70.5000000000000;73.5000000000000;76.5000000000000;79.5000000000000;82.5000000000000;85.5000000000000;88.5000000000000;91.5000000000000;94.5000000000000;97.5000000000000;100.500000000000];  
y_param = [0.750000000000000;2.32500000000000;4.05750000000000;5.96350000000000;8.06000000000000;10.3655000000000;12.9020000000000;15.6925000000000;18.7615000000000;22.1375000000000;25.8515000000000;29.9365000000000;34.4300000000000];
[xgrid, ygrid] = meshgrid(x_param, y_param);


% VFHS_loglin_transition
VFHS = 10.^repmat([ones(1,17), 3*ones(1,17)],13,1);
VFHS_trans = 10.^repmat([ones(1,9), linspace(1,3,16), 3*ones(1,9)],13,1);
VFHS_smoothsharp = VFHS_trans;
VFHS_smoothsharp(:,18:end) = repmat(VFHS_trans(:,11),1,34-17);
VFHS_transblock = 10.^repmat([2*ones(1,9), linspace(2,log10(500),16), log10(500)*ones(1,9)],13,1);
VFHS_2trans = 10.^repmat([2*ones(1,9), linspace(2,2.2,8), linspace(2.2,log10(500),9) log10(500)*ones(1,9)],13,1);
VFHS_2trans(:,18) = [];
VFHS_sqtrans = 10.^repmat([2*ones(1,9),2 + linspace(0,sqrt(log10(500) - 2),16).^2, log10(500)*ones(1,9)],13,1);
VFHS_cubtrans = 10.^repmat([2*ones(1,9),2 + linspace(0,(log10(500) - 2).^(1/3),16).^3, log10(500)*ones(1,9)],13,1);
VFHS_slant =  500*(xgrid./max(max(xgrid)) + 0.5*ygrid./max(max(ygrid)));
VFHS_slant(VFHS_slant>500) = 500;VFHS_slant(VFHS_slant<100) = 100;

gaussian_inclusion = 10*ones(13,34);
mu = [35, 8];   % position of gaussian centre (m)
sigma = 6;      % std, in m. Note - will not be a true gaussian due to block lengthening
A = 15;         % Amplitude of gaussian inclusion above background (ohm m)
gaussian_inclusion = gaussian_inclusion + A*exp(-((xgrid-mu(1)).^2 + (ygrid - mu(2)).^2)/(2*sigma.^2));

% Create a fault surface
depth = 20;
width = 10;
faultline = -[-depth*ones(1,width), -depth + ((1:(length(x_param)-width))*sqrt(depth)/(length(x_param)-width)).^2];
% faultline = faultline/max(faultline);
gf_height = 7;
gf_width = 15;
gf_offset = 10;
gf_centre = x_param(length(x_param)/2);
gf = gf_offset + -gf_height*exp(-((x_param-gf_centre).^2)/(2*gf_width.^2))';


% figure(1)
% plot(gf)

VFHS_slant_fault = VFHS_slant;
VFHS_slant_fault(ygrid>repmat(faultline,length(y_param),1)) = 200;

VFHS_slant_gf = VFHS_slant;
VFHS_slant_gf2 = VFHS_slant;
VFHS_slant_gf(ygrid>repmat(gf,length(y_param),1)) = 200;
VFHS_slant_gf2(ygrid>repmat(gf,length(y_param),1)) = 1000;

VFHS_sqtrans_gf2 = VFHS_sqtrans;
VFHS_sqtrans_gf2(ygrid>repmat(gf,length(y_param),1)) = 1000;

