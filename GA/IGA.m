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
repnum=0;


x=x0;% individual initialization
[f,x]=Fitness(x);
mergf=f;
mergx=x;
fgbest=zeros(itermax,size(f,2));

for iter=1:itermax      
    epsilon = epsilonControl(f(:,1),iter,Tc);    
    gama=0.8*(xmax-xmin)-(iter-1)/(itermax-1)*0.5*(xmax-xmin);
    crox= GACrossover(f,x,mergf,mergx,pc,gama,epsilon);
    mutx = GAMutation(crox,pm,xmin,xmax,iter,itermax);
    if repnum>popsize*0.01
        mutx=mutx+(rand(1,popsize)<pm).*(randi(2,D,popsize)-1).*trnd(1,D,popsize);
    end
    mutx=(mutx>xmax).*xmax+(mutx<xmin).*xmax+(mutx<=xmax&mutx>=xmin).*mutx;
    [newf,newx]=Fitness(mutx);
    mergf=[f;newf];
    mergx=[x,newx];
    idx=BestRank(mergf,epsilon);
    fgbest(iter,:)=mergf(idx(1),:);
    xgbest=mergx(:,idx(1));
    repnum=ismember(mergx',xgbest','rows');
    
    [f,x]= GASelection(mergf,mergx, epsilon); 


    
        
disp([num2str(iter),':',num2str(fgbest(iter,:))]);
    
    
    
    
    
end

end


