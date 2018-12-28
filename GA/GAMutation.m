function newx = GAMutation(x,p,xmin,xmax,iter,itermax)
[D,popsize]=size(x);
% X'(iDim)=X(iDim) + ��( t,UB-X(iDim) ) �����������0
% X(iDim) �C ��( t,X(iDim)-LB ) �����������1
% 
% ���� ��(t,y) = y*( 1 �C r^((1-t/T)^b) )
% 
% ��ʽ��t��ʾ��ǰ������T����������r��(0,1)����Ȳ��������С����b��ϵͳ�����������㷨������ѹ����
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


