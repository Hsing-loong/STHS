function R = Initial(~)
global params;
[T,Nh]=size(params.I);
inflow=params.I;
R=zeros(T,Nh);
for i=1:Nh
    if i<3
        
    elseif i<4
        inflow(:,i)=inflow(:,i)+...
            [zeros(params.Td(1),1);R(1:T-params.Td(1),1)]+...
            [zeros(params.Td(2),1);R(1:T-params.Td(2),2)];
    else
        inflow(:,i)=inflow(:,i)+...
            [zeros(params.Td(3),1);R(1:T-params.Td(3),3)];
    end
    s=params.Vini(i)+sum(inflow(:,i))-params.Vend(i);
    x = sort(rand(T-1,1));
    x = vertcat(x(1,:),diff(x,1,1),1-x(end,:)); % the colums of x sum to 1
    R(:,i) = params.Qmin(:,i)+(s-sum(params.Qmin(:,i)))*x; % shift-and-scale x to have column-sum of s
end
% R0=R0(:);
V=zeros(T+1,Nh);
V(1,:)=params.Vini;
V(end,:)=params.Vend;
for t=1:T
    V(t+1,:)=V(t,:)+inflow(t,:)-R(t,:);
end

end

