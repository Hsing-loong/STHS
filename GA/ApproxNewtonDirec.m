function newtondirection = ApproxNewtonDirec(f0,x0,f,x,gama,epsilon)
f=[f0;f];
x=[x0,x];
[f,idx,~]=unique(f,'rows','stable');
x=x(:,idx);
idx=BestRank(f,epsilon);
rankx=x(:,idx);
k=find(idx==1);
if k==1
    xworse=rankx(:,2:end);
    disworse=pdist2(x0',xworse');
    [disworse,index]=sort(disworse);
    xw=xworse(:,index(1));
    deltax=(gama+abs(xw-x0))*0.5;
    xb=x0-deltax;
elseif k==length(idx)
    xbetter=rankx(:,1:end-1);
    disbetter=pdist2(x0',xbetter');
    [disbetter,index]=sort(disbetter);
    xb=xbetter(:,index(1));
    deltax=(gama+abs(x0-xb))*0.5;
    xw=x0+deltax;
else
    xworse=rankx(:,k+1:end);
    xbetter=rankx(:,1:k-1);
    disworse=pdist2(x0',xworse');
    [disworse,index]=sort(disworse);
    xw=xworse(:,index(1));
    disbetter=pdist2(x0',xbetter');
    [disbetter,index]=sort(disbetter);
    xb=xbetter(:,index(1));
    deltax=(abs(xw-x0)+abs(x0-xb))*0.5;
end
newtondirection=0.5*deltax.*(xw-xb)./(xw-2*x0+xb+1e-12);
% distance=pdist2(x0',x');
% [distance,idx]=sort(distance);
% idx=idx(1:num);
% xneighbor=x(:,idx);
% fneighbor=f(idx,:);
% [rankf,rankx]=BestRank(fneighbor,xneighbor,epsilon);
% xbetter=rankx(:,1);
% xworse=rankx(:,end);
% if xbetter==xgbest
%     deltax=(gama+abs(xworse-x0))*0.5;
%     xb=x0-deltax;
%     xw=xworse;
% elseif xworse==xgworst
%     deltax=(gama+abs(x0-xbetter))*0.5;
%     xw=x0+deltax;
%     xb=xbetter;
% else
%     deltax=(abs(xworse-x0)+abs(x0-xbetter))*0.5;
%     xw=xworse;
%     xb=xbetter;
% end
% newtondirection=0.5*deltax.*(xw-xb)./(xw-2*x0+xb);


end

