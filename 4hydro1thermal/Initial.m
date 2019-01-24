function x0 = Initial(~)
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
    s=params.Vini(i)+sum(inflow(:,i)-params.Qmin(:,i))-params.Vend(i);
        
    x = [0;sort(rand(T-1,1));1];
    x = diff(x);
    R(:,i) = params.Qmin(:,i)+s*x; 
end
x0=R(:);
% V=zeros(T+1,Nh);
% V(1,:)=params.Vini;
% V(end,:)=params.Vend;
% for t=1:T
%     V(t+1,:)=V(t,:)+inflow(t,:)-R(t,:);
% end
% V=V(2:end-1,:);
% logic=V>params.Vmax|V<params.Vmin;
% V=logic.*(params.Vmin+rand(T-1,Nh).*(params.Vmax-params.Vmin))+(~logic).*V;
% x0=V(:);

% V=params.Vmin+rand(T-1,Nh).*(params.Vmax-params.Vmin);
% V0=zeros(T,Nh);
% V0(1,:)=params.Vini;
% V0(2:end,:)=V;
% V1=zeros(T,Nh);
% V1(1:end-1,:)=V;
% V1(end,:)=params.Vend;
% R=zeros(T,Nh);
% for i=1:Nh
%     if i<3
%         
%     elseif i<4
%         inflow(:,i)=inflow(:,i)+...
%             [zeros(params.Td(1),1);R(1:T-params.Td(1),1)]+...
%             [zeros(params.Td(2),1);R(1:T-params.Td(2),2)];
%     else
%         inflow(:,i)=inflow(:,i)+...
%             [zeros(params.Td(3),1);R(1:T-params.Td(3),3)];
%     end
%     VV=zeros(T+1,1);
%     VV(1)=params.Vini(i);
%     VV(end)=params.Vend(i);
%     
%     for t=1+randperm(T-1)
%         backwardV=VV;
%         backwardV(1:t)=0;
%         [tl,~,Vl]=find(backwardV);
%         tl=tl(1);
%         Vl=Vl(1);
%         forwardV=VV;
%         forwardV(t:end)=0;
%         [tu,~,Vu]=find(forwardV);
%         tu=tu(end);
%         Vu=Vu(end);
%         
%         Vl=min(max(params.Vmin(t-1,i),Vl-sum(inflow(t:tl-1,i)-params.Qmin(t:tl-1,i))),params.Vmax(t-1,i));
%         Vu=max(min(params.Vmax(t-1,i),Vu+sum(inflow(tu:t-1,i)-params.Qmin(tu:t-1,i))),params.Vmin(t-1,i));
% %         if V(t-1,i)<Vl
% %             V(t-1,i)=Vl;
% %         elseif V(t-1,i)>Vu
% %             V(t-1,i)=Vu;
% %         else
% %         end
%         
%         if V(t-1,i)<Vl||V(t-1,i)>Vu
%             V(t-1,i)=Vl+rand*(Vu-Vl);
%         end
%         VV(t)=V(t-1,i);
%     end
%     V0(:,i)=VV(1:end-1);
%     V1(:,i)=VV(2:end);
% 
%     R(:,i)=V0(:,i)+inflow(:,i)-V1(:,i);
% end
% x0=V(:);
end

