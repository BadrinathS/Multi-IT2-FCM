function  [V_final,center_l,center_r,u] = Multi_GT2_FCM(data,para)
    
    [no_patterns ,no_features] = size(data);

    %define center properly
    if para.cent_option == 1
        V_new = rand(para.c,no_features); %perfectly random center
    elseif para.cent_option == 2
        V_new = [0.5312 0.7970; 0.3551 0.5242]; %manual center
    else
        %Only true for Iris dataset
        V_new = [];
        for i=1:para.c
            V_new = [V_new;data(randi([50*(i-1)+1,50*i]),:)];
        end
    end
    
    %define the appropriate directions
    if para.dir_option == 1  %random 100 directions.
        dir = rand(100,no_features);
    elseif para.dir_option == 2 %all the unit directions
        dir = eye(no_features);
    else
        if no_features == 2 %if no of features are 2 then go for all directions
            dir = [];
            for i=0:pi/180:pi/2
                val = [cos(i) sin(i)];
                dir = [dir;val];
            end
        end
    end
    
%     disp(dir);

    %Define Fuzzifier
    m_min = para.m1;
    m_max = para.m2;
    no_alpha = para.no_alpha;
    alpha = rand(1,no_alpha);
    [m1_arr, m2_arr] = myfuzz1(m_min,m_max,alpha);
%     disp(m1_arr);
%     disp(m2_arr);
    
    
%     %test for IT2
%     m1_arr = [2];
%     m2_arr = [5];
%     alpha = 1;
%     no_alpha = 1;
%     %Test was successful.Both GenT2 & IT2 gave same results.

    %initializations
    itr = 0;
    V_left = zeros(para.c,no_features,no_alpha);
    V_right = zeros(para.c,no_features,no_alpha);
    V_old = zeros(para.c,no_features);
    Vl = zeros(para.c,no_features);
    Vr = zeros(para.c,no_features);
    u_left = zeros(para.c,no_patterns,no_alpha);
    u_right = zeros(para.c,no_patterns,no_alpha);
    Ul = zeros(para.c,no_patterns);
    Uh = zeros(para.c,no_patterns);
    u_l = zeros(para.c,no_patterns);
    u_r = zeros(para.c,no_patterns);
    h_m1_m2_idx = zeros(para.c,no_patterns);
    l_m1_m2_idx = zeros(para.c,no_patterns);

%     %center path
%     c1 = [];
%     c2 = [];
%     c3 = [];
% 
%     c1 = [c1;V_new(1,:)];
%     c2 = [c2;V_new(2,:)];
%     c3 = [c3;V_new(3,:)];


    %plotting
    % scatter(data([1:50],1),data([1:50],2),25,'k');
    % hold on;
    % scatter(data([51:100],1),data([51:100],2),25,'k');
    % scatter(data([101:150],1),data([101:150],2),25,'k');

    while(1)
        itr = itr + 1;
        V_old = V_new;
        temp1 = zeros(para.c,no_features);
        temp2 = zeros(para.c,no_patterns);

        %plotting
    %     plot(V_new(1,1),V_new(1,2),'sb','MarkerSize',5,'LineWidth',3);
    %     plot(V_new(2,1),V_new(2,2),'sr','MarkerSize',5,'LineWidth',3);
    %     plot(V_new(3,1),V_new(3,2),'sg','MarkerSize',5,'LineWidth',3);

        for k=1:no_alpha
            m1 = m1_arr(k);
            m2 = m2_arr(k);
            m = (m1 + m2)/2;
    %         disp(m1);
    %         disp(m2);
            %Do IT2 FCM calculations.

            %Estimating primary membeships
            [U_m1,~] = cal_U(data,V_old,m1,para.c);
            [U_m2,~] = cal_U(data,V_old,m2,para.c);

            %Determine Upper and lower membership using U_m1 and U_m2
            for i=1:para.c
                for j=1:length(data)
                    if( U_m1(i,j) > U_m2(i,j) )
                        Uh(i,j) = U_m1(i,j);
                        Ul(i,j) = U_m2(i,j);
                        h_m1_m2_idx(i,j) = m1;
                        l_m1_m2_idx(i,j) = m2;
                    else
                        Ul(i,j) = U_m1(i,j);
                        Uh(i,j) = U_m2(i,j);
                        h_m1_m2_idx (i,j) = m2;
                        l_m1_m2_idx(i,j) = m1;
                    end
                end
            end

            %Type reduction using KM algorithm
            for j=1:para.c
%                 for d=1:no_features
%                     [~,Vl(j,d),Vr(j,d),~,~] = KM_withchange(data(:,d)',Ul(j,:),Uh(j,:),h_m1_m2_idx(j,:),l_m1_m2_idx(j,:));
                    [~,Vl(j,:),Vr(j,:),~,u_l(j,:),u_r(j,:),~,~] = Multi_KM(data,Ul(j,:),Uh(j,:),m,dir);
%                     [V_final,V_l,V_r,U_final,U_l,U_r,centers_l,centers_r]

%                 end
            end

            V_left(:,:,k) = Vl;
            V_right(:,:,k) = Vr;
            u_left(:,:,k) = u_l;
            u_right(:,:,k) = u_r;
        end

        %Final defuzzified center and memebership values
        for i=1:no_alpha
            temp1 = temp1 + alpha(i) .* (V_left(:,:,i) + V_right(:,:,i));
            temp2 = temp2 + alpha(i) .* (u_left(:,:,i) + u_right(:,:,i));
        end
        V_new = temp1 ./ (2*sum(alpha));
        u = temp2 ./ (2*sum(alpha));

%         c1 = [c1;V_new(1,:)];
%         c2 = [c2;V_new(2,:)];
%         c3 = [c3;V_new(3,:)];


    %     V_new = sum((V_left + V_right),3) / (2*no_alpha);

        disp(norm(V_new - V_old));

        % step ------> 4
        if( norm(V_new - V_old) < para.eps || itr > 100 )
            break;
        end
    end

    fprintf("Final centers are:");
    disp(V_new);
    hard_partition(data,u);

%     % Plotting the center path
%     % plot(c1(1,1),c1(1,2),'sb','MarkerSize',5,'LineWidth',3);
%     plot(c1(:,1),c1(:,2),'c+-');
%     plot(c1(end,1),c1(end,2),'xr','MarkerSize',5,'LineWidth',1);
% 
%     % plot(c2(1,1),c2(1,2),'sb','MarkerSize',5,'LineWidth',3);
%     plot(c2(:,1),c2(:,2),'c+-');
%     plot(c2(end,1),c2(end,2),'xr','MarkerSize',5,'LineWidth',1);
% 
%     % plot(c3(1,1),c3(1,2),'sb','MarkerSize',5,'LineWidth',3);
%     plot(c3(:,1),c3(:,2),'c+-');
%     plot(c3(end,1),c3(end,2),'xr','MarkerSize',5,'LineWidth',1);
end
function [m1_arr ,m2_arr] = myfuzz1(m_min,m_max,alpha)
    x = [m_min; (m_min + m_max)/2; m_max];
    y = [0; 1 ;0];
    f=fit(x,y,'poly2');
    %plotting coefficients
    a = f.p1;
    b = f.p2;
    c = f.p3;
    c = c - alpha;
    m1_arr = (-b - sqrt(b*b - 4*a*c)) / (2 *a);
    m2_arr = (-b + sqrt(b*b - 4*a*c)) / (2 *a);
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

function [] = hard_partition(X,u)
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
    



