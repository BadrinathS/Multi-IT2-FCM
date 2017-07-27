function [V_final,V_l,V_r,U_final,U_l,U_r,centers_l,centers_r] = Multi_KM(X_const,Wl_const,Wr_const,m,dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remember y = [1 * d] (d = no. of dimensions)
%          X = [n * d] (n = no. of pattens)
%          Wl,Wr = [1 * n]
%          dv = [1 * d]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    no_dir = size(dir,1); 
    centers = zeros(no_dir,size(X_const,2));
    centers_l = zeros(no_dir,size(X_const,2));
    centers_r = zeros(no_dir,size(X_const,2));

    for t=1:no_dir
        
        %direction 
        dv = dir(t,:);
        
        %preprocessing on weights 
        Wr = Wr_const .^ m; 
        Wl = Wl_const .^ m;
        X = X_const;

        % sorting the points on the basis of distance
        D = abs (dv * X' ) ./ sqrt(dv * dv');
        X = [X D'];
        [X,index] = sortrows(X,size(X,2));
        X(:,size(X,2)) = [];

        Wl = Wl(index);
        Wr = Wr(index);
        
        u_left = Wl;
        u_right = Wl;


        % step--------> 2 and 3
        yl=X(end,:);
        l=0;
        a=Wl*X;
        b=sum(Wl);
        while ( l<length(X) && (yl*dv' > X(l+1,:)*dv') )
            l=l+1;
            a=a+(Wr(l)-Wl(l))*X(l,:);
            b=b+Wr(l)-Wl(l);
            yl=a./b;
            u_left(l) = Wr(l);
        end

        % Similarly calculating for yr
        % step--------> 2
        r=length(X); 
        yr=X(1,:);
        a=Wl*X; 
        b=sum(Wl);
        while (r>0 && (yr*dv' <  X(r,:)*dv') )
            a=a+(Wr(r)-Wl(r))*X(r,:);
            b=b+Wr(r)-Wl(r);
            yr=a./b;
            r=r-1;
            u_right(r) = Wr(r); 
        end

        y=(yl+yr)/2;
        centers(t,:) = y;
        centers_l(t,:) = yl;
        centers_r(t,:) = yr;
        
        u_above(index) = u_right;
        u_below(index) = u_left;
        u(t,:) = (u_above + u_below)/2;
        u_r(t,:) = u_above;
        u_l(t,:) = u_below;
        
%         fprintf("iteration no. is %d and center is\n",t);
%         disp(y);
%         fprintf("\n\n");
    %     plot(V_new(1,1),V_new(1,2),'xb','MarkerSize',15,'LineWidth',3);
    %     plot(V_new(2,1),V_new(2,2),'xr','MarkerSize',15,'LineWidth',3);
    %     plot(V_new(3,1),V_new(3,2),'xg','MarkerSize',15,'LineWidth',3);
    end
    
    V_final = sum(centers,1) ./ no_dir;
    V_l = sum(centers_l,1) ./ no_dir;
    V_r = sum(centers_r,1) ./ no_dir;
    U_l = sum(u_l,1) ./ no_dir;
    U_r = sum(u_r,1) ./ no_dir;
    U_final = sum(u,1) ./ no_dir;
    
end