clc,clear;
global params;
params=load('Data_4Hydro1Thermal.mat');
[T,Nh]=size(params.I);
params.DecVarTrans1=@(a) reshape(a,[],Nh);
params.DecVarTrans2=@(a) a(:);
params.ObFunc=@ObFunc_4Hydro1Thermal;

params.PDZ=false;
params.VP=false;



% w=sum(params.I)+params.Vini-params.Vend;
% w(3)=w(3)+w(1)-sum(params.Qmin(end-params.Td(1)+1:end,1))+...
%      w(2)-sum(params.Qmin(end-params.Td(2)+1:end,2));
% w(4)=w(4)+w(3)-sum(params.Qmin(end-params.Td(3)+1:end,3));
% 
% Rmax=params.Qmin;
% Imax=params.I;
% for i=1:Nh
%     if i<3
%         
%         
%     elseif i<4
%         Imax(:,i)=Imax(:,i)+...
%             [zeros(params.Td(1),1);Rmax(1:T-params.Td(1),1)]+...
%             [zeros(params.Td(2),1);Rmax(1:T-params.Td(2),2)];
%     else
%         Imax(:,i)=Imax(:,i)+...
%             [zeros(params.Td(3),1);Rmax(1:T-params.Td(3),3)];
%     end
%     Vmax=[params.Vini(i);zeros(T-1,1)];
%     Vmin=[zeros(T-1,1);params.Vend(i)];
%     for t=2:T
%         Vmax(t)=Vmax(t-1)+Imax(t,i)-params.Qmin(t,i);
%         Vmax(t)=min(Vmax(t),params.Vmax(t-1,i));
%     end
%     for t=T-1:-1:1
%         Vmin(t)=Vmin(t+1)-Imax(t+1,i)+params.Qmin(t+1,i);
%         Vmin(t)=max(Vmin(t),params.Vmin(t,i));
%     end  
%     Rmax(:,i)=Vmax+Imax(:,i)-Vmin;    
% end
% sumR=sum(params.I)+params.Vini-params.Vend;
% sumR(3)=sumR(3)+sumR(1)+sumR(2);
% sumR(4)=sumR(4)+sumR(3);
% summinQ=sum(params.Qmin);
% params.Rmax=zeros(T,Nh);
% for t=1:T
%     params.Rmax(t,:)=sumR-summinQ+params.Qmin(t,:);
% end




popsize=300;
itermax=2000;
for i=1:20
xmin=params.Qmin(:);
xmax=params.Qmax(:);
x0=xmin+rand(length(xmin),popsize).*(xmax-xmin);
[xgbest,fgbest] = IGA(x0,xmin,xmax,popsize,itermax);


R=reshape(xgbest,[],4);
[obvalue_viol,R,V,Q,Ph,SP,Ps]= ObFunc_4Hydro1Thermal(R);

save(['IGA--',num2str(fgbest(end,2)),'.mat']);
end

