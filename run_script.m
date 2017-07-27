%setting IT2 and multi_IT2 parameters
IT2_para.c = 2;
IT2_para.eps = 10^(-5);
IT2_para.m1 = 2;
IT2_para.m2 = 7;
IT2_para.m = 3;
IT2_para.cent_option = 2;

Multi_IT2_para.c = 2;
Multi_IT2_para.eps = 10^(-5);
Multi_IT2_para.m1 = 2;
Multi_IT2_para.m2 = 7;
Multi_IT2_para.m = 3;
Multi_IT2_para.cent_option = 2;
Multi_IT2_para.dir_option = 1;

Multi_GT2_para.c = 3;
Multi_GT2_para.eps = 10^(-5);
Multi_GT2_para.m1 = 2;
Multi_GT2_para.m2 = 7;
Multi_GT2_para.m = 4;
% Multi_GT2_para.no_alpha = 15;
Multi_GT2_para.cent_option = 3;
Multi_GT2_para.dir_option = 1;





data = csvread("./DATA/Iris.csv");
% data(:,[1,2]) = [];
data = ( data-min(data) ) ./ ( max(data) - min(data) );


%prepocessing the data
data1 = data([1:50],:);
data2 = data([51:100],:);
data3 = data([101:150],:);

shuffledArray1 = data1(randperm(50),:);
shuffledArray2 = data2(randperm(50),:);
shuffledArray3 = data3(randperm(50),:);

data = [shuffledArray1;shuffledArray2;shuffledArray3];



% exp2(IT2_para,Multi_IT2_para);
% Multi_GT2_FCM(data,Multi_GT2_para);

function exp2(IT2_para,Multi_IT2_para)
 x1 = normrnd(3,2,500,1);
 y1 = normrnd(8,6,500,1);
 
 x2 = normrnd(20,2,500,1);
 y2 = normrnd(10,6,500,1);
 
 x3 = normrnd(37,2,500,1);
 y3 = normrnd(12,6,500,1);
 
 figure(1);
 scatter(x1,y1,10,'c','filled');
 hold on;
 scatter(x2,y2,10,'g','filled');
 scatter(x3,y3,10,'b','filled');
 
 x = [x1;x2;x3];
 y = [y1;y2;y3];
 data = [x y];
 [center,cent_l,cent_r] = Multi_IT2_FCM_withchange(data,Multi_IT2_para);
 
 poly = zeros(362,size(cent_l,2),Multi_IT2_para.c);
 
 for i=1:Multi_IT2_para.c
        poly(:,:,i) = cat(1,cent_l(:,:,i),cent_r(:,:,i));
 end
 
 
 for i=1:Multi_IT2_para.c
     scatter(center(i,1),center(i,2),8,'r','filled');
     x =  poly(:,1,i);
     y =  poly(:,2,i);
     k = convhull(x,y);
     plot(x(k),y(k),'r-');
 end

 xlabel("x1");
 ylabel("x2");
 title("Estimating centroid region using Multi");
%  dim = [.2 .5 .3 .3];
% str = 'Centroid region is represented by red region';
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
disp(center);



 
hold off;
 
end
function exp1(IT2_para,Multi_IT2_para)
    %creating own dataset
    a = 0.2;
    b = 0.4;
    x1 = (b-a)*rand(100,1) + a;

    a = 0.6;
    b = 0.8;
    x2 = (b-a)*rand(100,1) + a;

    a = 0.4;
    b = 0.5;
    y1 = (b-a)*rand(100,1) + a;


    x = [x1;x2];
    y = [y1;y1];
    data = [x y];

    figure(1);
    scatter(data(:,1),data(:,2),25,'b','filled');
    xlim([0 1]);
    ylim([0 1]);
    hold on;

%     [center1,u1] = with_cent_passed(data,IT2_para);
%     disp(center1);
%     [center2,~,~,u2] = Multi_IT2_FCM_withchange(data,Multi_IT2_para);
%     disp(center2);


%     center = [mean(x1) mean(y1);mean(x2) mean(y1)];
%     scatter(center(1,1),center(1,2),15,'r','filled');
%     scatter(center(2,1),center(2,2),15,'r','filled');

    arr1 = [];
    arr2 = [];

    for i=0:1:10
        x_n = rand(2,1);

        a = 0;
        b = 0.2;
        y1_n = (b-a)*rand(2,1) + a;

        a = 0.8;
        b = 1;
        y2_n = (b-a)*rand(2,1) + a;

%         figure(1);
%         scatter([x_n;x_n],[y1_n;y2_n],10,'g','d','filled');

    %     adding noise to data.
        x = [x ;x_n ;x_n];
        y = [y ;y1_n ;y2_n];

        data = [x y];
        [cent1,u1] = with_cent_passed(data,IT2_para);
        perf1 = hard_partition(x1,x2,y1,u1);
        [cent2,~,~,u2] = Multi_IT2_FCM_withchange(data,Multi_IT2_para);
        perf2 = hard_partition(x2,y1,u2);


        arr1 = [arr1;i+1 norm(cent1-center)];
        arr2 = [arr2;i+1 norm(cent2-center)];
        
        pause(3);
    end


    figure(2);
    plot(arr1(:,1),arr1(:,2),'-or','MarkerIndices',1:1:length(arr1),'LineWidth',1);
    hold on;
    plot(arr2(:,1),arr2(:,2),'-ob','MarkerIndices',1:1:length(arr2),'LineWidth',1);
    xlim([1 11]);
    ylim([0 0.02]);
    xlabel("Test cases #");
    ylabel("Error in center");
    title("IT2 FCM vs Multi IT2 FCM");
    lgd = legend('IT2 FCM','Multi IT2 FCM','Location','northwest');
    lgd.FontSize = 13;
    % disp(lgd);
    hold  off;
end

function perf =  hard_partition(X,u)
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

%     check3 = cent_test([101:150],1);
%     temp1 = check3(check3 == mode(check3));
%     temp2 = temp1(temp1 ~= mode(check2));
%     temp3 = temp2(temp2 ~= mode(check1));
%     count = count + length(temp3);
    
    perf = count/length(X)*100;

    fprintf('Performance of algorithm is %f',count/length(X)*100);

end


% functioncent_diff(cent1,cent2)
% 
% 
% 
% disp(mean(x1));
% disp(mean(y1));
% center = [mean(x1) mean(y1);mean(x2) mean(y1)];
% 
% scatter(x,y,15,'k');
% hold on;
% scatter(center(1,1),center(1,2),15,'r','filled');
% scatter(center(2,1),center(2,2),15,'r','filled');
% xlim([0 1]);
% ylim([0 1]);






