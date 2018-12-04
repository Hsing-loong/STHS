function epsilon = epsilonControl(voil,iter,Tc)
cpmin=3;
cpmax=10;
theta=0.9;
popsize=length(voil);
lambda=1-sum(voil>0)/popsize;
cp=cpmin+lambda*(cpmax-cpmin);
epsilon0=voil(floor(theta*popsize));
if iter>Tc
    epsilon=0;
else
    epsilon=epsilon0*(1-iter/Tc)^cp;
end


end

