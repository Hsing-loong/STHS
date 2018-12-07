function newtondirection = ApproxNewtonDirec(x0,x,gama) %#codegen
[D,num]=size(x);
dist=sqrt(sum((x0-x).^2));
k=find(dist<eps);
if k==1
    [~,index]=min(dist(k+1:end));
    xw=x(:,index+k);
    deltax=(gama+abs(xw-x0))*0.5;
    xb=x0+deltax;
elseif k==num
    [~,index]=min(dist(1:k-1));
    xb=x(:,index);
    deltax=(gama+abs(x0-xb))*0.5;
    xw=x0-deltax;
else
    [~,index]=min(dist(k+1:end));
    xw=x(:,index+k);
    [~,index]=min(dist(1:k-1));
    xb=x(:,index);
    deltax=(abs(xw-x0)+abs(x0-xb))*0.5;
end
newtondirection=0.5*deltax.*(xw-xb)./(xw-2*x0+xb+1e-12);

end

