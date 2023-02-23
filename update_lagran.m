function [input]=update_lagran(itr,ip_cnt,tmp_flag,input)
 
 
 
 
   %---------------------LAGRANGIAN VALUE AND REDUCTIONS----------------------   
 % Keep original values for ip inversion  
if (ip_cnt==1 && itr==1) 
    input.original_lagrn_min=input.lagrn_min; % ACB
    input.original_lagrn_max=input.lagrn_max; %ACB
    input.original_lagrn=input.lagrn;   %Classic
end
if (ip_cnt==2 && itr==1) 
    input.lagrn_min=input.original_lagrn_min;
    input.lagrn_max=input.original_lagrn_max;
    input.lagrn=input.original_lagrn;
end

if itr>1 && tmp_flag==1
    

    % Ip4di cool
%         if itr==1 ;input.lagrn_reduction=2;end
%         if itr==2 ;input.lagrn_reduction=1.75; end
%         if itr==3 ;input.lagrn_reduction=1.5; end
%         if itr>4 ;input.lagrn_reduction=1.25; end        
%         if itr>5  ;input.lagrn_reduction=1; end

% % slow cool
% input.lagrn_reduction = 1.1;

% slow_cool2
% if itr < 5
%     input.lagrn_reduction = 1.1;
% else
%       input.lagrn_reduction = 1.5;
% end
        

%my scheme from april: 'TGV cool2'
%             if itr==2 ;input.lagrn_reduction=2;end
%             if itr==3 ;input.lagrn_reduction=2; end
%             if itr==4 ;input.lagrn_reduction=1.75; end
%             if itr>4 && itr<8;input.lagrn_reduction=1.5; end
%             if itr>8  ;input.lagrn_reduction=1.25; end

% TGV_cool3
% if itr==2 ;input.lagrn_reduction=1;end
% if itr==3 ;input.lagrn_reduction= 2; end
% if itr==4 ;input.lagrn_reduction=1.75; end
% if itr>4 && itr<8;input.lagrn_reduction=1.5; end
% if itr>8  ;input.lagrn_reduction=1.25; end

% strong_cool
% input.lagrn_reduction = 2;

% old tgv_cool
%             if itr<8 ;input.lagrn_reduction=1.5; end
%             if itr>8  ;input.lagrn_reduction=1.25; end            

% % gradual damping
%         input.lagrn_reduction=1.25;
%         if itr>6  ;input.lagrn_reduction=1; end

%     if mod(itr,2)==1
%         input.lagrn_reduction=1.5;
%     else
%         input.lagrn_reduction=1;
%     end
%     

% ip4di
% red_factor = [2, 1.75, 1.5, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25];

% strong cool
red_factor = [ 1, 1, 1, repmat(2, 1, (itr))];

% % strong_step2_cool
% if mod( itr, 2)
%     red_factor = [ 1, repmat( 2, 1, (itr))];
% else 
%     red_factor = ones( 1, itr);
% end

% strong_step4_cool
% if mod( itr, 4)
%     red_factor = ones( 1, itr);
% else 
%     red_factor = [ 1, repmat( 2, 1, (itr))];
% end

% strong step3_cool
% if mod( itr, 3)
%     red_factor = ones( 1, itr);
% else 
%     red_factor = [ 1, repmat( 2, 1, (itr))];
% end

% root2_cool
% red_factor = [1,  repmat(sqrt(2), 1, (itr))];

% root2_step2_cool
% if mod( itr, 2)
%     red_factor = [ 1, repmat(sqrt(2), 1, (itr))];
% else 
%     red_factor = ones( 1, itr);
% end

% root3 cool
% red_factor = [ sqrt(3), repmat(2, 1, (itr))];

% strongish cool 
% red_factor = [ 1, repmat(2, 1, 4), repmat(1.5, 1, 20) ];

% accelerated1_cool
% red_factor = [ 1, 2, 2, 5, repmat(1.5, 1, 20)];

% accelerated2_cool
% red_factor = [ 1, 2, 2, 5, 5, 5, 1.5, 1.5, 1.5, repmat(1.05, 1, 20)];

% slow2_cool
% red_factor = [ 1, repmat(1.5, 1, 50)];

% slow3_cool
% red_factor = [ 1, repmat(1.25, 1, 50)];

% step1_cool
% red_factor = [ repmat( 1, 1, 4), 10, 1, 1, 10, repmat(1, 1, 100)];

% step2_cool
% red_factor = [ repmat( 1, 1, 5), 10, 1, 1, 1, 10, repmat(1, 1, 100)];


input.lagrn_reduction = red_factor(itr);    


   

input.lagrn_min=input.lagrn_min/input.lagrn_reduction;
input.lagrn_max=input.lagrn_max/input.lagrn_reduction;
input.lagrn=input.lagrn/input.lagrn_reduction;

lag_thresh = 1e-4;
if input.lagrn < lag_thresh
    input.lagrn = lag_thresh;
end

  
% clifton cooling scheme:
% custom_lagrn = [0.1, 0.2, 0.01, 0.0022, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005]; clifton
% custom cool 1:
% custom_lagrn = [0.2,  0.001, 0.001, 0.0005, 0.0005, 0.0001, 0.0001, 0.00005, 0.00005, 0.00001];
% gb2 cool
% custom_lagrn = [0.2, 0.003686633654376, 0.001073175883198, 3.805753885785665e-04, 1.654065451074255e-04, 9.883109356351321e-05, 6.235148904680434e-05, 6.235148904680434e-05];
% custom_lagrn = [0, 0.001073175883198, 3.805753885785665e-04, 1.654065451074255e-04, 9.883109356351321e-05, 6.235148904680434e-05, 6.235148904680434e-05];

% custom_lagrn = [0.200000000000000, 0.460747924251852, 0.004336510379885, 0.001141434284835, repmat(1e-3, 1, 5)];
% custom_lagrn = [0, 0.001141434284835, repmat(1e-3, 1, 5)];

% not_cool
% custom_lagrn = [0.1, 0.1, 0.05, 0.02, 0.1, 0.15, 0.02, 0.01, 0.005, 0.0005, 0.0001]; 

% input.lagrn = custom_lagrn(itr);    

end 

end