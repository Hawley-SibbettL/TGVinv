
load('p_calc test workspace')



if itr ~= input.tgv_itr
    for i = 2:(input.sp2_itr+1)
        
            [p2] = p_calc(mesh, lagrn, input.tgv_lagrn, gamma_c, gamma_p);
        
        % spearate px, py
        mesh.px = p2(1:length(p2)/2);
        mesh.py = p2((length(p2)/2+1):length(p2));
        % eliminate ghost points
        %                         mesh.px = mesh.px(domain_ind);
        %                         mesh.py = mesh.py(domain_ind);
        p2 = [mesh.px; mesh.py];
        
        % Store results
        px(:,i) = mesh.px;
        py(:,i) = mesh.py;
        
        % root mean square difference between last and current
        % solution normalise by rms p value
        %                     rms_p(i) = sqrt(mean((p2 - p1).^2))/(0.5*sqrt(mean(p2.^2 + p1.^2));
        
        
        
        rms_p(i-1) = sqrt(mean((p2 - p1).^2))/(0.5*sqrt(mean( (mesh.cx*log10(res_param)).^2 + (mesh.cy*log10(res_param)).^2)));
        %         disp(['rms dp = ',num2str(p_rms(i))])
        %                 rms_change(i) = abs((p_rms(i) - p_rms(i-1)))./(p_rms(i) + p_rms(i-1)); % percentage rms change per itr
        
        disp(['tgv iteration ',num2str(i-1),' : rms dp = ',num2str(rms_p(i-1))])
        
        
        if ((i > 4) && (rms_p(i-1) < input.p_conv)) %  stop if change smaller than p_conv after first few itr
            sp2_i = i;  % store number of iterations needed
            break
        elseif (i > 20 && (rms_p(i-1) > rms_p(i-2))) % stop if diverging
            sp2_i = i-1;  % store number of iterations needed
            mesh.px = p1(1:length(p1)/2);
            mesh.py = p1((length(p1)/2+1):length(p1));
            break
        end
        
        
        if i == (input.sp2_itr+1)
            sp2_i = i;
        end
        
        p1 = p2;
        
    end
end


rms_p(i:end) = [];

% Plot convergence
figure(10)
hold off
plot(rms_p)
title('rms change in p')
xlabel('sp 2 iteration')
ylabel('rms dp')