function  [V_final,u] = with_cent_passed(data,para,num_cluster)

    [no_patterns ,no_features] = size(data);

    %define center properly
    if para.cent_option == 1
        V_new = rand(num_cluster,no_features); %perfectly random center
    elseif para.cent_option == 2
        V_new = [0.5312 0.7970; 0.3551 0.5242]; %manual center
    else
        %Only true for Iris dataset
        V_new = [];
        for i=1:num_cluster
            V_new = [V_new;data(randi([40*(i-1)+1,40*i]),:)];
        end
    end
    
    %initializations
    V_old = zeros(num_cluster,no_features);
    Vl = zeros(num_cluster,no_features);
    Vr = zeros(num_cluster,no_features);
    u_l = zeros(num_cluster,no_patterns);
    u_r = zeros(num_cluster,no_patterns);
    h_m1_m2_idx = zeros(num_cluster,no_patterns);
    l_m1_m2_idx = zeros(num_cluster,no_patterns);
    itr = 0;

%     %center path
%     c1 = [];
%     c2 = [];
% %     c3 = [];
% 
%     c1 = [c1;V_new(1,:)];
%     c2 = [c2;V_new(2,:)];
% %     c3 = [c3;V_new(3,:)];

    diff_array = [];

    while(1)
        itr = itr + 1;
        old_diff = norm(V_new - V_old);
        diff_array = [diff_array old_diff];
        V_old = V_new;

    %     step ------> 1
    %     Estimating FCM membership using m1
        [U_m1,~] = cal_U(data,V_old,para.m1,num_cluster);


    %     Estimating FCM membership using m2
        [U_m2,~] = cal_U(data,V_old,para.m2,num_cluster);


        % step ------> 2
        %Determine Upper and lower membership using U_m1 and U_m2
        for i=1:num_cluster
            for j=1:length(data)
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
        for j=1:num_cluster
            temp_l = zeros(1,no_patterns);
            temp_r = zeros(1,no_patterns);
            for d=1:no_features
    %             [V_new(j,d),Vl(j,d),Vr(j,d),u_left,u_right] = KM_withchange(data(:,d)',Ul(j,:),Uh(j,:),h_m1_m2_idx(j,:),l_m1_m2_idx(j,:));
                [V_new(j,d),Vl(j,d),Vr(j,d),u_left,u_right] = KM(data(:,d)',Ul(j,:),Uh(j,:),para.m);

                temp_l = temp_l + u_left;
                temp_r = temp_r + u_right;
            end
            u_l(j,:) = temp_l ./ no_features;
            u_r(j,:) = temp_r ./ no_features;
        end

%         c1 = [c1;V_new(1,:)];
%         c2 = [c2;V_new(2,:)];
%         c3 = [c3;V_new(3,:)];

        new_diff = norm(V_new - V_old); 
        disp(norm(V_new - V_old));

        % step ------> 4
        if( new_diff > old_diff || itr>100 || new_diff < para.eps )
            break;
        end
    end

    fprintf("Final centers are:");
    disp(V_new);
    V_final = V_new;
    % disp(u_l);
    % disp(u_r);
    u = (u_l + u_r) ./ 2;
    % hard_partition(data,u);

%     figure(2);
%     x = [1:1:itr];
%     y = diff_array;
%     plot(x,y,'-or','MarkerIndices',1:1:length(diff_array),'LineWidth',1);
%     % ylim([0 0.6]);
%     hold on;
% 
%     figure(1);
%     plot(c1(:,1),c1(:,2),'b+-');
%     plot(c1(end,1),c1(end,2),'xb','MarkerSize',5,'LineWidth',1);
% 
%     % plot(c2(1,1),c2(1,2),'sb','MarkerSize',5,'LineWidth',3);
%     plot(c2(:,1),c2(:,2),'b+-');
%     plot(c2(end,1),c2(end,2),'xb','MarkerSize',5,'LineWidth',1);
% 
%     % plot(c3(1,1),c3(1,2),'sb','MarkerSize',5,'LineWidth',3);
% %     plot(c3(:,1),c3(:,2),'b+-');
% %     plot(c3(end,1),c3(end,2),'xb','MarkerSize',5,'LineWidth',1);
%     % hold off;
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
    






