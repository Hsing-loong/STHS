function  [obvalue_viol,R,V,Q,Ph,SP,Ps]= ObFunc_4Hydro1Thermal(R)
global params;
inflow=params.I;
[T,Nh]=size(inflow);
V=zeros(T+1,Nh);
V(1,:)=params.Vini;
V(end,:)=params.Vend;
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
    
    array=1:T;
    m=1;
    while m<100
        VT=params.Vini(i)+sum(inflow(:,i)-R(:,i));
        deltaV=VT-params.Vend(i);
        if abs(deltaV)<1e-5
            break;
        end
        
        l=length(array);
        k=randi(l);
        index=array(randperm(l,k));
        % Using the adjacent differences of sorted uniform RVs
        if k==1
            x=rand;
        else
            x = sort(rand(k-1,1));
            x = vertcat(x(1,:),diff(x,1,1),1-x(end,:)); % the colums of x sum to 1
        end
        deltaR = deltaV*x; % shift-and-scale x to have column-sum of s
        
%         index=array;
%         deltaR=deltaV/length(array);
        if deltaV<0
            R(index,i)=max(R(index,i)+deltaR,params.Qmin(index,i));
            array=find(R(:,i)>params.Qmin(:,i));
        else
            R(index,i)=R(index,i)+deltaR;
            array=1:T;
        end
        m=m+1;
    end
end
for t=1:T
    V(t+1,:)=V(t,:)+inflow(t,:)-R(t,:);
end
V=V(2:end,:);
a=params.C2;
b=params.C3.*V+params.C5;
c=params.C1.*(V.^2)+params.C4.*V+params.C6;

%note
%1.Ph can not reach params.Phmax with any V1
%2.There is a certain Q that makes Ph reach params.Phmin
%3.params.Qmin is less than extrpoint with any V1
discriminant=b.^2-4*a.*(c-params.Phmin);
extrpoint=-0.5*b./a;
Q1=0.5*(-b+sqrt(discriminant))./a;
Q2=0.5*(-b-sqrt(discriminant))./a;
Q=Q1.*(R<Q1)+R.*(R>=Q1&R<=Q2)+Q2.*(R>Q2);
Q=params.Qmin.*(Q<params.Qmin)+Q.*(Q>=params.Qmin&Q<=params.Qmax)+params.Qmax.*(Q>params.Qmax);
Q=Q.*(Q<extrpoint)+extrpoint.*(Q>=extrpoint);
if params.PDZ
    logic=Q<params.Qpu&Q>params.Qpl;
    Q=Q.*(~logic)+params.Qpl.*(logic);
end
SP=R-Q;
Ph=params.C1.*(V.^2)+params.C2.*(Q.^2)+params.C3.*(V.*Q)+...
    params.C4.*V+params.C5.*Q+params.C6;
Ps=params.load-sum(Ph,2);
if params.VP
    cost=sum(params.a.*Ps.^2+params.b.*Ps+params.c+abs(params.e.*sin(params.f.*(params.Psmin-Ps))));
else
    cost=sum(params.a.*Ps.^2+params.b.*Ps+params.c);
end

VminViol=max(0,params.Vmin-V(1:end-1,:));
VmaxViol=max(0,V(1:end-1,:)-params.Vmax);
VendViol=max(0,abs(V(end,:)-params.Vend));
PsminViol=max(0,params.Psmin-Ps);
PsmaxViol=max(0,Ps-params.Psmax);
viol=[VminViol(:);VmaxViol(:);VendViol(:);PsminViol(:);PsmaxViol(:)];
viol=sum(viol);
obvalue_viol=[viol,cost];
end
