function [newx] = GACrossover(self,selx,f,x,p,gama,epsilon)
[Dx,popsize]=size(selx);
[~,Df]=size(self);
parentsidx=randperm(popsize);
parentsx=selx(:,parentsidx);
parentsf=self(parentsidx,:);

[uniparentsx,unia,unic]=unique(parentsx','rows','stable');
uniparentsx=uniparentsx';
uniparentsf=parentsf(unia,:);

fun=@(a,b) ApproxNewtonDirec(a,b,f,x,gama,epsilon);
parentsndir=cell2mat(cellfun(fun,num2cell(uniparentsf,2)',num2cell(uniparentsx,1),'UniformOutput',false));
parentsndir=parentsndir(:,unic);

parentsx1=parentsx(:,1:popsize/2);
parentsndir1=parentsndir(:,1:popsize/2);
parentsx2=parentsx(:,popsize/2+1:end);
parentsndir2=parentsndir(:,popsize/2+1:end);


iscross=binornd(1,p,1,popsize/2);
x1=parentsx1-iscross.*(rand(1,popsize/2).*parentsndir1+0.5*rand(1,popsize/2).*(parentsx2-parentsx1));
x2=parentsx2-iscross.*(rand(1,popsize/2).*parentsndir2+0.5*rand(1,popsize/2).*(parentsx1-parentsx2));



newx=[x1,x2];


% parentsx=mat2cell(parentsx,Dx,2*ones(1,popsize/2));
% parentsf=mat2cell(parentsf,2*ones(1,popsize/2),Df)';
% 
% newx=cellfun(@(a,b) Crossover(a,b,f,x,p,gama,epsilon),parentsf,parentsx,'UniformOutput',false);
% newx=cell2mat(newx);






end
function offspringsx=Crossover(parentsf,parentsx,f,x,p,gama,epsilon)
[D,~]=size(parentsx);
offspringsx=parentsx;
if binornd(1,p)
    %%
%     r=rand(D,2);
%     offspringsx=r.*parentsx(:,1)+(1-r).*parentsx(:,2);
    %%
    newtondirection1=ApproxNewtonDirec(parentsf(1,:),parentsx(:,1),f,x,gama,epsilon);
    newtondirection2=ApproxNewtonDirec(parentsf(2,:),parentsx(:,2),f,x,gama,epsilon);
    
    x1=parentsx(:,1)-rand*newtondirection1+rand*(parentsx(:,2)-parentsx(:,1));
    x2=parentsx(:,2)-rand*newtondirection2+rand*(parentsx(:,1)-parentsx(:,2));
    offspringsx=[x1,x2];
    %%
    
end
end








