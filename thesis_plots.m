%% Periodic bc plots

p_title = {'Blocky Model', 'Smooth Model', 'Gaussian model'}; % titles for p plots
p_file = {'block_rpy1_sharp', 'gauss_rpy1', 'plume_layer_rpy2_sharp'};
p_folder = ['D:\TGV_revision\'];

% objective plots
for i = 2
    
    load(fullfile(p_folder, p_file{i}))

    figure(1)
    subplot_tight(3,3, 3 + i)
    plot(tgv_lagrn,TGV_t1)
    hold on
    plot(tgv_lagrn,TGV_t2)
    plot(tgv_lagrn,TV_reg.*ones(1,length(tgv_lagrn)))
    legend('TGV |\nabla m - p|_{l1}','TGV 0.5|\nabla p + \nabla p^T|_{l1}','|\nabla m|_{L1}')% 'total TGV objective',
    title(p_title{i})
    xlabel('\mu')
    ylabel(['Magnitude of regularisation term'])

end