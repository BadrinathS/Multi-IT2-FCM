function [V,U] = FCM(X,m1,U0,eps,c)
% Remember U = [c * N](N = no. of data points , c = No. of clusters )
%          V = [c * D](D = no. of dimensions)

% step -------> 1
% Randomly initialize the cluster membership values, ?ij.(U0)
U = U0;
J_old = 0;
itr = 0;
while(1)
    itr = itr + 1;
    % step -------> 2
    % Calculate the cluster centers:
    V = (U .^ (m1)) * X ./ sum( (U .^ (m1)) , 2 ); 


    % step -------> 3
    % Update ?ij according to the following:
    [U,dist] = cal_U(X,V,m1,c);

    % step -------> 4
    % Calculate the objective function, Jm.
    J_new = obje_fn(dist,U,m1,c);
%     fprintf("Iteration no. %d and obj. fn is %f\n",itr,J_new);

    % step -------> 5
    if( (abs(J_new - J_old) < eps) || (itr > 30) )
        break;
    else
        J_old = J_new;
    end
  
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

function J = obje_fn(dist,U,m,c)
    temp = 0;
    for i=1:length(dist)
        for j=1:c
            temp = temp + (U(j,i)^ m) * (dist(i,j) ^ 2);
        end
    end
    J = temp;           
end