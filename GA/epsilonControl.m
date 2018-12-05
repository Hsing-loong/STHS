function epsilon = epsilonControl(epsilon0,iter,Tc)
cpmin=3;
Tlambda=0.95*Tc;
elambda=1e-5;
cp=max(cpmin,(log(elambda)-log(epsilon0))/(1-Tlambda/Tc));
if iter>Tlambda
    cp=0.3*cp+0.7*cpmin;
    epsilon=epsilon0*(1-iter/Tc)^cp;
elseif iter>Tc
    epsilon=0;
else
    epsilon=epsilon0*(1-iter/Tc)^cp;
end


end

