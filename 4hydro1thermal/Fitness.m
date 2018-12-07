function [f,x] = Fitness(x) 
global funhandle;
x=num2cell(x,1);
V=cellfun(funhandle.DecVarTrans1,x,'UniformOutput',false);
[obvalue_viol,V,~,~,~,~,~]=cellfun(funhandle.ObFunc,V,'UniformOutput',false);
x=cellfun(funhandle.DecVarTrans2,V,'UniformOutput',false);
x=cell2mat(x);
f=cell2mat(obvalue_viol');
end

