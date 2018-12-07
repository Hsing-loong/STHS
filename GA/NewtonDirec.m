function newtondirection = NewtonDirec(f0,x0) %#codegen
error=0.0005;
num=numel(x0);
x0=repmat(x0,1,num);
x=[x0+error*eye(num),x0-error*eye(num)];
[f,~]=GradFitness(x);
f1=f(1:num,2);
f2=f(num+1:end,2);

firstd=0.5*(f1-f2)/error;
secondd=(f1-2*f0(2)+f2)/error^2;
newtondirection=firstd./(secondd+eps);
end





