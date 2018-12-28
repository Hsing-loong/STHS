function newx = GAMutation(x,p,xmin,xmax,iter,itermax)
[D,popsize]=size(x);
% X'(iDim)=X(iDim) + Δ( t,UB-X(iDim) ) ，随机变量是0
% X(iDim) – Δ( t,X(iDim)-LB ) ，随机变量是1
% 
% 其中 Δ(t,y) = y*( 1 – r^((1-t/T)^b) )
% 
% 上式中t表示当前代数，T是最大进化，r是(0,1)间均匀产生的随机小数，b是系统参数，决定算法的收敛压力。
r2=rand(D,popsize);
delta1=(xmax-x).*(1-rand(D,popsize).^(1-iter/itermax));
delta2=-(x-xmin).*(1-rand(D,popsize).^(1-iter/itermax));
delta= ((r2>=0.5).*delta1+(r2<0.5).*delta2);

newx=x+(rand(1,popsize)<p).*(randi(2,D,popsize)-1).*delta;

% Cauchy(0,1)


[~,IA,~]=uniquetol(newx',1e-1,'ByRows',true,'DataScale',1);
idex=setdiff((1:D)',IA);
num=numel(idex);
if num/popsize>0.01   
%     newx(:,idex)=newx(:,idex)+(randi(2,D,num)-1).*randn(D,num);
    newx(:,idex)=newx(:,idex)+(randi(2,D,num)-1).*(0+1*tan(pi*(rand(D,num)-0.5)));
end

end


