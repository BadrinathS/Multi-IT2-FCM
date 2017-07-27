function  [V_final,center_l,center_r,u] = Multi_IT2_FCM_withchange(data,para)
    
    [no_patterns ,no_features] = size(data);

    if para.cent_option == 1
        V_new = rand(para.c,no_features);
    else
        V_new = [0.5312 0.7970; 0.3551 0.5242];
    end
        
    
    %initializations
    Vl = zeros(para.c,no_features);
    Vr = zeros(para.c,no_features);
    V_old = zeros(para.c,no_features);
    h_m1_m2_idx = zeros(para.c,no_patterns);
    l_m1_m2_idx = zeros(para.c,no_patterns);
    u = zeros(para.c,no_patterns);
    center_l = zeros(91,no_features,para.c);
    center_r = zeros(91,no_features,para.c);
    itr=0;

%     %center path
%     c1 = [];
%     c2 = [];
%     % c3 = [];
% 
%     c1 = [c1;V_new(1,:)];
%     c2 = [c2;V_new(2,:)];
%     % c3 = [c3;V_new(3,:)];

    mult_diff_array = [];

    while(1)
        itr = itr + 1;
        old_diff = norm(V_new - V_old);
        mult_diff_array = [mult_diff_array old_diff];
        V_old = V_new;

    %     step ------> 1
    %     Estimating FCM membership using m1
        [U_m1,~] = cal_U(data,V_old,para.m1,para.c);


    %     Estimating FCM membership using m2
        [U_m2,~] = cal_U(data,V_old,para.m2,para.c);


        % step ------> 2
        %Determine Upper and lower membership using U_m1 and U_m2
        for i=1:para.c
            for j=1:no_patterns
                if( U_m1(i,j) > U_m2(i,j) )
                    Uh(i,j) = U_m1(i,j);
                    Ul(i,j) = U_m2(i,j);
                    h_m1_m2_idx(i,j) = para.m1;
                    l_m1_m2_idx(i,j) = para.m2;
                else
                    Ul(i,j) = U_m1(i,j);
                    Uh(i,j) = U_m2(i,j);
                    h_m1_m2_idx (i,j) = para.m2;
                    l_m1_m2_idx(i,j) = para.m1;
                end
            end
        end

    %     U = (Uh + Ul) / 2;

        % step ------> 3
        %Type reduction and Defuzzification using KM algorithm
        for j=1:para.c
    %         [V_new(j,:),~,~] = Multi_KM_withchange(data,Ul(j,:),Uh(j,:),h_m1_m2_idx(j,:),l_m1_m2_idx(j,:));
            [V_new(j,:),~,~,u(j,:),cent_l,cent_r] = Multi_KM(data,Ul(j,:),Uh(j,:),para.m);
            center_l(:,:,j) = cent_l;
            center_r(:,:,j) = cent_r;
        end

%         c1 = [c1;V_new(1,:)];
%         c2 = [c2;V_new(2,:)];
%     %     c3 = [c3;V_new(3,:)];

        new_diff = norm(V_new - V_old); 
    %     disp(norm(V_new - V_old));

        % step ------> 4
        if( new_diff > old_diff || itr>100 || new_diff < para.eps )
            disp(itr);
            break;
        end
    end
    
    V_final = V_new;

%     figure(1);
%     % plot(c1(1,1),c1(1,2),'sb','MarkerSize',5,'LineWidth',3);
%     plot(c1(:,1),c1(:,2),'g+-');
%     plot(c1(end,1),c1(end,2),'xr','MarkerSize',5,'LineWidth',1);
% 
%     % plot(c2(1,1),c2(1,2),'sb','MarkerSize',5,'LineWidth',3);
%     plot(c2(:,1),c2(:,2),'g+-');
%     plot(c2(end,1),c2(end,2),'xr','MarkerSize',5,'LineWidth',1);
% 
%     % plot(c3(1,1),c3(1,2),'sb','MarkerSize',5,'LineWidth',3);
%     % plot(c3(:,1),c3(:,2),'g+-');
%     % plot(c3(end,1),c3(end,2),'xr','MarkerSize',5,'LineWidth',1);

    fprintf("The ultimate center is\n");
    disp(V_new);

    % hard_partition(data,u);

%     figure(2);
%     x = [1:1:itr];
%     y = mult_diff_array;
%     plot(x,y,'-ob','MarkerIndices',1:1:length(mult_diff_array),'LineWidth',1);
%     % ylim([0 0.6]);
%     hold on;
%     xlabel("No. of iterations");
%     ylabel("Error in center");
%     title(" IT2 FCM vs Multi IT2 FCM ");

end


function [U,dist] = cal_U(X,V,m,c)
dist = zeros(length(X),c);
U = zeros(c,length(X));
% Calculating distance of jth point from ith center
    for i=1:c
        for j=1:length(X)
            dist(j,i) = norm( X(j,:) - V(i,:) );
        end
    end
%     Calclulating uij
    for i=1:c
        for j=1:length(X)
            if( dist(j,i) == 0 )
                U(i,j) = 1;
                continue;
            end
            
            temp = 0;
            done = 0;
            for k=1:c
                if( dist(j,k) == 0 )
                    U(i,j) = 0;
                    done = 1;
                    break;
                end
                
                temp = temp + ( dist( j,i) / dist(j,k) )^ (2 /(m-1));
            end
            
            if ( ~done == 1 )
                U(i,j) = 1/temp;
            end
        end
    end
            
            
            
end

function hard_partition(X,u)
    cent_test = zeros(length(X),1);
    count = 0;
    for i=1:length(X)
        [~,cent_test(i)] = max(u(:,i));
    end
    
    check1 = cent_test([1:50],1);
    count = count + length(find(check1 == mode(check1)));

    check2 = cent_test([51:100],1);
    temp1 = check2(check2 == mode(check2));
    temp2 = temp1(temp1 ~= mode(check1));
    count = count + length(temp2); 

    check3 = cent_test([101:150],1);
    temp1 = check3(check3 == mode(check3));
    temp2 = temp1(temp1 ~= mode(check2));
    temp3 = temp2(temp2 ~= mode(check1));
    count = count + length(temp3);

    fprintf('Performance of algorithm is %f',count/length(X)*100);

end
    
    





