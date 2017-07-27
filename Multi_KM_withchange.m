function [V_final,yl,yr] = Multi_KM_withchange(X_const,Wl_const,Wr_const,h_id,l_id)
% Remember y = [1 * d] (d = no. of dimensions)
%          X = [n * d] (n = no. of pattens)
%          Wl,Wr = [1 * n]
%          dv = [1 * d]
%     eps = 10^(-3);

    dv = zeros(1,size(X_const,2));
    var = 100;
    centers = zeros(100,size(X_const,2));
    incr = pi/180;
    itr = 0;
    for t=0:incr:pi/2
        itr = itr + 1;
%         for d=1:size(X_const,2)
%             dv(d) = rand();
%         end

        dv(1) = cos(t);
        dv(2) = sin(t);
        
        %preprocessing on weights 
        Wr = Wr_const .^ h_id;
        Wl = Wl_const .^ l_id;
        X = X_const;

        % sorting the points on the basis of distance
        D = abs (dv * X' ) ./ sqrt(dv * dv');
        X = [X D'];
        [X,index] = sortrows(X,size(X,2));
        X(:,size(X,2)) = [];

%         for i=1:size(X,2)
%             X(:,i) = X(index,i);
%         end

        Wl = Wl(index);
        Wr = Wr(index);


        % step--------> 2 and 3
        yl=X(end,:);
        l=0;
        a=Wl*X;
        b=sum(Wl);
        while ( l<length(X) && (yl*dv' > X(l+1,:)*dv') )
%             disp(yl*dv' - X(l+1,:)*dv');
            l=l+1;
            a=a+(Wr(l)-Wl(l))*X(l,:);
            b=b+Wr(l)-Wl(l);
            yl=a./b;
        end
%         disp(yl*dv' - X(l+1,:)*dv');
        % Similarly calculating for yr
        % step--------> 2
        r=length(X); 
        yr=X(1,:);
        a=Wl*X; 
        b=sum(Wl);
        while (r>0 && (yr*dv' <  X(r,:)*dv') )
%             disp(yr*dv' <  X(r,:)*dv');
            a=a+(Wr(r)-Wl(r))*X(r,:);
            b=b+Wr(r)-Wl(r);
            yr=a./b;
            r=r-1;
        end

        y=(yl+yr)/2;
        centers(itr,:) = y;
%         fprintf("iteration no. is %d and center is\n",t);
%         disp(y);
%         fprintf("\n\n");
%         plot(yl(1,1),yl(1,2),'ob','MarkerSize',5,'LineWidth',1);
%         hold on;
%         plot(yr(1,1),yr(1,2),'or','MarkerSize',5,'LineWidth',1);
%         plot(y(1,1),y(1,2),'og','MarkerSize',5,'LineWidth',1);

        %     plot(V_new(3,1),V_new(3,2),'xg','MarkerSize',15,'LineWidth',3);
    end
    V_final = sum(centers,1) ./ itr;
    
end