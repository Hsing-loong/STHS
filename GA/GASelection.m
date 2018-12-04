function [self,selx]= GASelection( f,x, epsilon)
% GA select operator
% Tournament selection
[D,popsize]=size(x);
% toursize=fix(0.6*popsize);
toursize=2;
selsize=popsize/2;

index=zeros(toursize,selsize);
for i=1:selsize
index(:,i)=randperm(popsize,toursize);
end
index=num2cell(index,1);
fun1=@(a) f(a,:);
fun2=@(a) x(:,a);
tempf=cellfun(fun1,index,'UniformOutput',false);
tempx=cellfun(fun2,index,'UniformOutput',false);

fun3=@(a) BestRank(a,epsilon);
idx=cellfun(fun3,tempf,'UniformOutput',false);
fun4=@(a,b) a(b(1),:);
fun5=@(a,b) a(:,b(1));
self=cell2mat(cellfun(fun4,tempf,idx,'UniformOutput',false)');
selx=cell2mat(cellfun(fun5,tempx,idx,'UniformOutput',false));

end


