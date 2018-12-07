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
Tc=0.6*itermax;




x=x0;% individual initialization
[f,x]=Fitness(x);
newf=[];
newx=[];
fgbest=zeros(itermax,2);
% epsilon0=0.8*(sum(f(:,1))/popsize+min(f(:,1)));
epsilon0=0.6*max(f(:,1));
cp=log(epsilon0/1e-4);
% cpmin=cp-3;
% cpmax=cp+3;
for iter=1:itermax
%     lambda=1-sum(f(:,1)>0)/popsize;
%     cp=cpmin+lambda*(cpmax-cpmin);
    if iter>=Tc
        epsilon=0;
    else
        epsilon=epsilon0*exp(-cp*(iter/Tc));
    end
    
    mergf=[f;newf];
    mergx=[x,newx];
    
    idx=BestRank(mergf,epsilon);
    fgbest(iter,:)=mergf(idx(1),:);
    xgbest=mergx(:,idx(1));   
    
    x=mergx(:,idx(1:popsize));
    f=mergf(idx(1:popsize),:);

    [self,selx]= GASelection(f,x,epsilon);
    
    gama=0.6*(xmax-xmin)-(iter-1)/(itermax-1)*0.4*(xmax-xmin);  
    crox= GACrossover(self,selx,mergf,mergx,pc,gama,epsilon);
    
    mutx = GAMutation(crox,pm,xmin,xmax,iter,itermax);

    [newf,newx]=Fitness(mutx);
    
    dist=sqrt(sum((fgbest(iter,:)-f).^2,2));
    logic=dist<1;
    repnum=sum(logic);
    logic(1)=0;
    if repnum>popsize*0.9
        x(:,logic)=xgbest+(randi(2,D,repnum-1)-1).*randn(D,repnum-1);
        [f(logic,:),x(:,logic)]=Fitness(x(:,logic));
    end
    
    
    disp([num2str(iter),':',num2str([fgbest(iter,:),epsilon])]);
    
    
      
end

end


