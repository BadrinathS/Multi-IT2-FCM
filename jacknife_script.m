% need to define function for each dataset
avg = 0;
for i=1:10
   avg = avg + run_cluster_Iris();
end
fprintf("Avg performance is %f",avg/10);

function perf = run_cluster_Iris()
    data = csvread("./DATA/Iris.csv");
    % data(:,[1,2]) = [];
    
    %Normalization done on data
    data = ( data-min(data) ) ./ ( max(data) - min(data) );
    
    %Define no. of clusters
    num_clust = 3;
    
    %Set parameters for IT2,GT2,Multi_IT2,Multi_GT2
    IT2_para = set_IT2_params("Iris");              %<--- 1
    %Multi_IT2_para = set_Multi_IT2_params("Iris");  %<--- 2
    %GT2_para = set_GT2_params("Iris");              %<--- 3
    %Multi_GT2_para = set_Multi_GT2_params("Iris");  %<--- 4
     
    %Jacknifing starts here
    pre_data = preprocess_data(data,"Iris");
    count = 0;
    for w=1:5
        data_org = pre_data;
        test1 = data_org([10*(w-1)+1:10*w],:);
        test2 = data_org([10*(w-1)+51:10*w+50],:);
        test3 = data_org([10*(w-1)+101:10*w+100],:);
        data_test = [test1;test2;test3;];

        data_org([10*(w-1)+1:10*w],:) = [];
        data_org([10*(w-1)+41:10*w+40],:) = [];
        data_org([10*(w-1)+81:10*w+80],:) = [];
        
        [V,~] = with_cent_passed(data_org,IT2_para,num_clust);
        [U,~] = cal_U(data_test,V,IT2_para.m,num_clust);
        count = count + hard_partition(data_test,U); 
    end
    perf = count/length(pre_data)*100;
    fprintf('Performance of algorithm is %f',count/length(pre_data)*100);
end

function pre_data = preprocess_data(data,str)
    if(str == "Iris")
        data1 = data([1:50],:);
        data2 = data([51:100],:);
        data3 = data([101:150],:);

        shuffledArray1 = data1(randperm(50),:);
        shuffledArray2 = data2(randperm(50),:);
        shuffledArray3 = data3(randperm(50),:);

        pre_data = [shuffledArray1;shuffledArray2;shuffledArray3];
    end
end
function IT2_para = set_IT2_params(str)
    if(str == "Iris")
        IT2_para.eps = 10^(-5);
        IT2_para.m1 = 2;
        IT2_para.m2 = 7;
        IT2_para.m = 3;
        IT2_para.cent_option = 3;
    end
end
function Multi_IT2_para = set_Multi_IT2_params(str)
    if(str == "Iris")
        
    end
end
function GT2_para = set_GT2_params(str)
    if(str == "Iris")
        
    end
end
function Multi_GT2_para = set_Multi_GT2_params(str)
    if(str == "Iris")
        
    end
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
function count = hard_partition(X,u)
    cent_test = zeros(length(X),1);
    num_clust = size(u,1);
    count = 0;
    
    for i=1:length(X)
        [~,cent_test(i)] = max(u(:,i));
    end
    
    check1 = cent_test([1:10],1);
    count = count + length(find(check1 == mode(check1)));

    check2 = cent_test([11:20],1);
    temp1 = check2(check2 == mode(check2));
    temp2 = temp1(temp1 ~= mode(check1));
    count = count + length(temp2); 

    check3 = cent_test([21:30],1);
    temp1 = check3(check3 == mode(check3));
    temp2 = temp1(temp1 ~= mode(check2));
    temp3 = temp2(temp2 ~= mode(check1));
    count = count + length(temp3);

%     fprintf('Performance of algorithm is %f',count/length(X)*100);

end