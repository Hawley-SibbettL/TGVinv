% Script to plot the l-curves used in an inversion

filename = 'inv_results';
load(filename)

%%

lagr = final.lc_trial_lagr;     % Lagr values tested in L-curve
lc_lagr = final.lc_lagr;        % Resulting lagr from L-curve
lc_roughness = final.lc_roughness;
lc_misfit = final.lc_misfit;
lc_curvature = final.lc_curvature;
lagr_ind = final.lagr_ind;
% lc_curvature_smooth = final.lc_curvature_smooth;

n_im = length(lc_lagr);

for i = 1:n_im
    disp(['i = ',num2str(i)])
    
    figure(4)
    subplot(ceil(n_im/2),2,i)
    plot((lc_misfit(i,:)),(lc_roughness(i,:)),'bx')   
    hold on
    plot((lc_misfit(i,1:5)),(lc_roughness(i,1:5)),'r^')
    plot((lc_misfit(i,lagr_ind(i))),(lc_roughness(i,lagr_ind(i))),'go','markersize',10,'linewidth',5)
    hold off
    title(['L-curve. itr = ', num2str(i)]);
    xlabel('|log(d^{obs}) - log(d^{pre})|)')
    ylabel('(|(c_x + c_y)log(\bf{m})|)');
    
    figure(5)
    subplot(ceil(n_im/2),2,i)
    semilogx(lagr,lc_curvature(i,:),'rx-')
    hold on
%     semilogx(lagr,lc_curvature_smooth(i,:),'kx-')
    semilogx((lagr(lagr_ind(i))),(lc_curvature(i,lagr_ind(i))),'go','markersize',10,'linewidth',5)
    hold off
    title(['Curvature. itr = ', num2str(i)]);
    xlabel('lagr')
    ylabel('curvature');
    
    figure(6)
    subplot(ceil(n_im/2),2,i)
    [AX, H1, H2] = plotyy(lagr,lc_misfit(i,:),lagr,lc_roughness(i,:),'semilogx');
    hold(AX(1),'on'); hold(AX(2),'on');
    set(get(AX(1),'Ylabel'),'String','Misfit') ;
    set(get(AX(2),'Ylabel'),'String','Roughness');
    set(AX(1),'ycolor','m') ;
    set(AX(2),'ycolor','c');
    set(H1,'marker','x','markeredgecolor','m');
    set(H2,'marker','x','markeredgecolor','c');
    title(['misfit. itr = ', num2str(i)]);
    xlabel('lagr');
    ylabel('misfit');
    hold(AX(1),'off'); hold(AX(2),'off');

    
%     figure(7)
%     subplot(ceil(n_im/2),2,i)
%     semilogx(lagr,lc_roughness(i,:),'gx')
%     title(['roughness. itr = ', num2str(i)]);
%     xlabel('lagr')
%     ylabel('roughness');
%     
    
end
    