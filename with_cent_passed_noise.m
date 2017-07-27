% 0.0819    0.0726    0.0846    0.0765
% 0.0921    0.0800    0.1101    0.1191
% 0.0693    0.0904    0.0289    0.0141
%fetching data from file to matrix
% data_org=csvread('./Data/rand_square.csv');
% data_org(:,[1,2]) = [];
a = 0.2;
b = 0.4;
x1 = (b -a)*rand(150,1) + a;

a = 0.6;
b = 0.8;
x2 = (b -a)*rand(150,1) + a;

a = 0.4;
b = 0.6;
y1 = (b - a)*rand(150,1) + a;

x = [x1;x2];
y = [y1;y1];

noise_x = rand(30,1);

a = 0.1;
b = 0.3;
noise_y1 = (b -a)*rand(15,1) + a;

a = 0.7;
b = 0.9;
noise_y2 = (b -a)*rand(15,1) + a;

noise_y = [noise_y1;noise_y2];

x = [x;noise_x];
y = [y;noise_y];


data_org = [x y];

% data_org=normc(data_org);
% data_org = ( data_org-min(data_org) ) ./ ( max(data_org) - min(data_org) );

% plot(data_org([1:9],1),data_org([1:9],2),'ob');
% xlim([-5 5]);
% ylim([-5 5]);
% hold on;
% plot(data_org([10:end],1),data_org([10:end],2),'or');
% xlim([-5 5]);
% ylim([-5 5]);
% hold on;
% plot(data_org([101:150],1),data_org([101:150],2),'og');
% hold on;

scatter(x1,y1,20,'k');
xlim([0 1]);
ylim([0 1]);
hold on;
scatter(x2,y1,20,'k');
scatter(noise_x,noise_y,10,'k','d');
% scatter(data_org([51:100],1),data_org([51:100],2),25,'k');
% scatter(data_org([101:150],1),data_org([101:150],2),25,'k');


% No. of clusters(c)
c = 2;
eps = 10^(-5);
V_old = zeros(c,size(data_org,2));
V_new = [0.1 0.2;0.7 0.2];
% V_new = rand(c,size(data_org,2));
Vl = zeros(c,size(data_org,2));
Vr = zeros(c,size(data_org,2));
u_l = zeros(c,size(data_org,1));
u_r = zeros(c,size(data_org,1));
itr = 0;


%m,m1 and m2 values
m1=1.1;
m2=10;
% m=(m1+m2)/2;

%center path
c1 = [];
c2 = [];
% c3 = [];

c1 = [c1;V_new(1,:)];
c2 = [c2;V_new(2,:)];
% c3 = [c3;V_new(3,:)];

h_m1_m2_idx = zeros(c,size(data_org,1));
l_m1_m2_idx = zeros(c,size(data_org,1));

while(1)
    itr = itr + 1;
    old_diff = norm(V_new - V_old);
    V_old = V_new;
    
%     plot(V_new(1,1),V_new(1,2),'sb','MarkerSize',5,'LineWidth',3);
%     plot(V_new(2,1),V_new(2,2),'sr','MarkerSize',5,'LineWidth',3);
%     plot(V_new(3,1),V_new(3,2),'sg','MarkerSize',5,'LineWidth',3);
    
%     plot(V_new(1,1),V_new(1,2),'MarkerSize',5,'LineWidth',1);
%     plot(V_new(2,1),V_new(2,2),'MarkerSize',5,'LineWidth',1);
%     plot(V_new(3,1),V_new(3,2),'MarkerSize',5,'LineWidth',1);

%     step ------> 1
%     Estimating FCM membership using m1
    [U_m1,~] = cal_U(data_org,V_old,m1,c);


%     Estimating FCM membership using m2
    [U_m2,~] = cal_U(data_org,V_old,m2,c);


    % step ------> 2
    %Determine Upper and lower membership using U_m1 and U_m2
    for i=1:c
        for j=1:length(data_org)
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
    
%     U = (Uh + Ul) / 2;

    % step ------> 3
    %Type reduction and Defuzzification using KM algorithm
    for j=1:c
        temp_l = zeros(1,size(data_org,1));
        temp_r = zeros(1,size(data_org,1));
        for d=1:size(data_org,2)
            [V_new(j,d),Vl(j,d),Vr(j,d),u_left,u_right] = KM_withchange(data_org(:,d)',Ul(j,:),Uh(j,:),h_m1_m2_idx(j,:),l_m1_m2_idx(j,:));
            temp_l = temp_l + u_left;
            temp_r = temp_r + u_right;
        end
        u_l(j,:) = temp_l ./ size(data_org,2);
        u_r(j,:) = temp_r ./ size(data_org,2);
    end
    
    c1 = [c1;V_new(1,:)];
    c2 = [c2;V_new(2,:)];
%     c3 = [c3;V_new(3,:)];
    
%     flag1 = [1;2;3];
%     gscatter(V_new(:,1),V_new(:,2),flag1,['y' 'y' 'y']);
% 
%     
%     plot(V_new(3,1),V_new(3,2),'xg','MarkerSize',15,'LineWidth',3);
    
%     disp(V_new);
    new_diff = norm(V_new - V_old); 
    disp(norm(V_new - V_old));
    
    % step ------> 4
    if( new_diff > old_diff || itr>100 || new_diff < eps )
        break;
    end
end

fprintf("Final centers are:");
disp(V_new);
% disp(u_l);
% disp(u_r);
u = (u_l + u_r) ./ 2;
% hard_partition(data_org,u);

% plot(V_new(1,1),V_new(1,2),'xb','MarkerSize',10,'LineWidth',3);
% plot(V_new(2,1),V_new(2,2),'xr','MarkerSize',10,'LineWidth',3);

% plot(V_new(1,1),V_new(1,2),'sb','MarkerSize',5,'LineWidth',3);
% plot(V_new(2,1),V_new(2,2),'sr','MarkerSize',5,'LineWidth',3);
% plot(V_new(3,1),V_new(3,2),'sg','MarkerSize',5,'LineWidth',3);

% flag1 = [1;2;3];
% gscatter(V_new(:,1),V_new(:,2),flag1,['y' 'y' 'y']);
% for i=1:length(c1)
    
% plot(c1(1,1),c1(1,2),'sb','MarkerSize',5,'LineWidth',3);
plot(c1(:,1),c1(:,2),'r+-');
plot(c1(end,1),c1(end,2),'xr','MarkerSize',5,'LineWidth',1);

% plot(c2(1,1),c2(1,2),'sb','MarkerSize',5,'LineWidth',3);
plot(c2(:,1),c2(:,2),'r+-');
plot(c2(end,1),c2(end,2),'xr','MarkerSize',5,'LineWidth',1);

% plot(c3(1,1),c3(1,2),'sb','MarkerSize',5,'LineWidth',3);
% plot(c3(:,1),c3(:,2),'r+-');
% plot(c3(end,1),c3(end,2),'xr','MarkerSize',5,'LineWidth',1);
% hold off;


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
    






