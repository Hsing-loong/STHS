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

popsize=300;
itermax=2000;

xmin=params.Qmin(:);
xmax=params.Qmax(:);
f=@Initial;
% x=cell2mat(arrayfun(f,1:10000,'UniformOutput',false));
for i=1:20
    x0=cell2mat(arrayfun(f,1:popsize,'UniformOutput',false));
%     x0=xmin+rand(length(xmin),popsize).*(xmax-xmin);
%     x0=x(:,randperm(10000,popsize));
    [xgbest,fgbest] = IGA(x0,xmin,xmax,popsize,itermax);
    
    
    V=reshape(xgbest,[],4);
    [obvalue_viol,V,R,Q,Ph,SP,Ps]= ObFunc_4Hydro1Thermal(V);
    
    save(['IGA--',num2str(fgbest(end,2)),'.mat']);
end

