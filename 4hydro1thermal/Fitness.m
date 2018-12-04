function [f,x] = Fitness(x)
global params;
x=num2cell(x,1);
V=cellfun(params.DecVarTrans1,x,'UniformOutput',false);
[obvalue_viol,V,~,~,~,~,~]=cellfun(params.ObFunc,V,'UniformOutput',false);
x=cellfun(params.DecVarTrans2,V,'UniformOutput',false);
x=cell2mat(x);
f=cell2mat(obvalue_viol');
end

