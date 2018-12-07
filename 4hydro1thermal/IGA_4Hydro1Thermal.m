clc,clear;
global params;
params=load('Data_4Hydro1Thermal.mat');
params.PDZ=false;
params.VP=false;

[T,Nh]=size(params.I);

global funhandle;
funhandle.DecVarTrans1=@(a) reshape(a,[],Nh);
funhandle.DecVarTrans2=@(a) a(:);
funhandle.ObFunc=@ObFunc_4Hydro1Thermal;
% funhandle.Fitness=@Fitness;
funhandle.GradObFunc=@GradObFunc;


% sumR=sum(params.I)+params.Vini-params.Vend;
% sumR(3)=sumR(3)+sumR(1)+sumR(2);
% sumR(4)=sumR(4)+sumR(3);
% summinQ=sum(params.Qmin);
% params.Rmax=zeros(T,Nh);
% for t=1:T
%     params.Rmax(t,:)=sumR-summinQ+params.Qmin(t,:);
% end
% Imax=params.I;
% for i=1:Nh
%     if i<3
%
%
%     elseif i<4
%         Imax(:,i)=Imax(:,i)+...
%             [zeros(params.Td(1),1);params.Rmax(1:T-params.Td(1),1)]+...
%             [zeros(params.Td(2),1);params.Rmax(1:T-params.Td(2),2)];
%     else
%         Imax(:,i)=Imax(:,i)+...
%             [zeros(params.Td(3),1);params.Rmax(1:T-params.Td(3),3)];
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
%     params.Rmax(:,i)=min(params.Rmax(:,i),Vmax+Imax(:,i)-Vmin);
% end

Vmax=zeros(T-1,2);
Vmin=zeros(T-1,2);
for t=1:T-1
    if t<2
        Vmax(t,1:2)=params.Vini(1:2)+params.I(t,1:2)-params.Qmin(t,1:2);
    else
        Vmax(t,1:2)=Vmax(t-1,1:2)+params.I(t,1:2)-params.Qmin(t,1:2);
    end
    params.Vmax(t,1:2)=min(Vmax(t,1:2),params.Vmax(t,1:2));
end
for t=T-1:-1:1
    if t>T-2
        Vmin(t,1:2)=params.Vend(1:2)-params.I(t+1,1:2)+params.Qmin(t+1,1:2);
    else
        Vmin(t,1:2)=Vmin(t+1,1:2)-params.I(t+1,1:2)+params.Qmin(t+1,1:2);
    end
    params.Vmin(t,1:2)=max(Vmin(t,1:2),params.Vmin(t,1:2));
end




popsize=300;
itermax=2000;

xmin=params.Vmin(:);
xmax=params.Vmax(:);
% f=@Initial;
for i=1:20
%     x0=cell2mat(arrayfun(f,1:popsize,'UniformOutput',false));
    x0=xmin+rand(length(xmin),popsize).*(xmax-xmin);
    [xgbest,fgbest] = IGA(x0,xmin,xmax,popsize,itermax);
    
    
    V=reshape(xgbest,[],4);
    [obvalue_viol,V,R,Q,Ph,SP,Ps]= ObFunc_4Hydro1Thermal(V);
    
    save(['IGA--',num2str(fgbest(end,2)),'.mat']);
end

