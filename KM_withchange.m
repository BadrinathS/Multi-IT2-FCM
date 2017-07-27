function [y,yl,yr,u_left,u_right] = KM_withchange(X,Wl,Wr,h_id,l_id)

%Type reduction by using KM algorithm
% computing yl and yr
% Step ------> 1
    %peprocessing on weights 
    Wr = Wr .^ h_id;
    Wl = Wl .^ l_id;
    
    u_left = zeros(1,length(X));
    u_right = zeros(1,length(X));
    [X,index] = sort(X);
    Wl = Wl(index);
    Wr = Wr(index);
    itr1 = 0;
    itr2 = 0;
    
% Step ------> 2    
    W = (Wl + Wr)/2 ; 
    y_old = X*(W') / sum(W) ;
    
    while(1)
        itr1=itr1+1;
    % Step ------> 3   
        for i=1:length(X)-1
            if( (X(i)<y_old) &&  (y_old <= X(i+1)) )
                l = i;
                break;
            end
        end

    % Step ------> 4   
        for i=1:length(X)
            if( i <= l )
                W(i) = Wr(i);
            else
                W(i) = Wl(i);
            end
        end

        y_new = X*(W') / sum(W) ;

    % Step ------> 5 
        if( y_new == y_old || itr1>20)
            break;
        else
            y_old = y_new;
        end
    end
    
    yl = y_old;
    u_left(index) = W; 

    
%   For calculating yr step 1 remains the same hence we stat fom step 2

% Step ------> 2    
    W = (Wl + Wr)/2 ; 
    y_old = X*(W') / sum(W) ;
    
    while(1)
        itr2=itr2+1;
    % Step ------> 3   
        for i=1:length(X)-1
            if( (X(i)<y_old) &&  (y_old <= X(i+1)) )
                r = i;
                break;
            end
        end
        
    % Step ------> 4   
        for i=1:length(X)
            if( i <= r )
                W(i) = Wl(i);
            else
                W(i) = Wr(i);
            end
        end
        
        y_new = X*(W') / sum(W) ;
        
    % Step ------> 5 
        if( y_new == y_old  || itr2>20)
            break;
        else
            y_old = y_new;
        end
    end
    
    yr = y_old;
    u_right(index) = W;
%     fprintf("final switch points at l = %d,r=%d\n",l,r);
  
    
%     Defuzzification step
    y = ( yl + yr ) / 2;
    
end