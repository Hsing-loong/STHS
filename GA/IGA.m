function [ xgbest,fgbest ] = IGA(x0,xmin,xmax,popsize,itermax)
%The improved GA algorithm solves the minimum optimization problem
%xgbest fgbest Global optimum
%fitness fitness function
%xmax xmin  The lower and upper bounds of the solution space

%GA parameters
% popsize=300;
% itermax=2000;
D=length(xmin);%Number of dimensions
pc=0.8;
pm=0.05;
Tc=0.8*itermax;
theta=0.9;
cpmin=3;
cpmax=10;
Tlambda=0.95*Tc;
elambda=1e-5;
repnum=0;


x=x0;% individual initialization
[f,x]=Fitness(x);
viol=sortrows(f);
epsilon0=viol(floor(theta*popsize),1);
fgbest=zeros(itermax,2);
for iter=1:itermax



    
    viol=sortrows(f);
    
    lambda=1-sum(viol(:,1)>0)/popsize;
    cp=cpmin+lambda*(cpmax-cpmin);
    if iter>=Tc
        epsilon=0;
    else
        epsilon=epsilon0*exp(-cp*(iter/Tc));
    end
    
    [self,selx]= GASelection(f,x,epsilon);
    
    gama=0.8*(xmax-xmin)-(iter-1)/(itermax-1)*0.5*(xmax-xmin);  
    crox= GACrossover(self,selx,f,x,pc,gama,epsilon);
    
    mutx = GAMutation(crox,pm,xmin,xmax,iter,itermax);
    [newf,newx]=Fitness(mutx);
    mergf=[f;newf];
    mergx=[x,newx];

%     [unix,idx,~]=unique(mergx','rows','stable');
%     unix=unix';
%     unif=mergf(idx,:);
    
    idx=BestRank(mergf,epsilon);
    fgbest(iter,:)=mergf(idx(1),:);
    xgbest=mergx(:,idx(1));
    
    renum=sum(ismember(mergx',xgbest','rows'));
%     if renum>0
%         rex=xgbest+(randi(2,D,popsize)-1).*trnd(1,D,popsize);
%         [ref,rex]=Fitness(rex);
%         mergx(:,end+1:end+renum)=rex;
%         mergf(end+1:end+renum,:)=ref;
%     end 
    x=mergx(:,idx(1:popsize));
    f=mergf(idx(1:popsize),:);
    
    
    
    
    disp([num2str(iter),':',num2str([fgbest(iter,:),epsilon])]);
    
    
    
    
    
end

end


