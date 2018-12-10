
clearvars
% close all


inv_file = 'finningley';
inv_folder = 'Paper/new imager';
username = 'bgsvisluke2';

type = 'l1'; % Inversion type: 'l1'; l12'; 'l22'
lag = '15';
append_txt = '_30m';%
ref_flag = 0;
n = 1:10; % indices to plot

if ref_flag == 1; ref_text = '_ref'; else ref_text = ''; end;


file = ['C:\Users\',username,'\OneDrive - The University of Nottingham\UoN Box Migration\Coupled TGV results\',inv_folder,'\',inv_file,'_',type,'_',lag,'lag',ref_text,append_txt];

load(file,'final')

n_im = length(final.RMS); % number of iterations in file
n(n>n_im) = []; % removes non-existant iterations from n

x = unique(final.param_x);
y = unique(final.param_y);
[xgrid, ygrid] = meshgrid(x,y);

for i = n
    
    res_image = reshape(log10(final.res_param1(:,i)),length(x),length(y))';

    
    
    figure(21)
    subplot(ceil(n_im/2),2,find(n==i))
    surf(x,y,res_image,'edgecolor','none')
    view([0,0,1])
    title(['log(m), itr = ',num2str(i),', RMS = ',num2str(final.RMS(i))],'interpreter','none')
    set(gca,'ydir','reverse')
    colorbar
    axis image
    colormap parula
    
    
    
    
end
    
    
    
    
