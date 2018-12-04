% function [rankf,rankx]=BestRank(f,x,epsilon)
% tempfcompare=f;
% tempfcompare(tempfcompare(:,1)<epsilon,1)=0;
% [~,idx]=sortrows(tempfcompare);
% rankf=f(idx,:);
% rankx=x(:,idx);
% end
function idx=BestRank(f,epsilon)
f(f(:,1)<epsilon,1)=0;
[~,idx]=sortrows(f);
end
