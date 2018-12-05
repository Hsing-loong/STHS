function  [obvalue_viol,V,R,Q,Ph,SP,Ps]= ObFunc_4Hydro1Thermal(V)  %#codegen
global params;
inflow=params.I;
[T,Nh]=size(inflow);
V(V>params.Vmax)=params.Vmax(V>params.Vmax);
V(V<params.Vmin)=params.Vmin(V<params.Vmin);
% logic=V>params.Vmax|V<params.Vmin;
% V=logic.*(params.Vmin+rand(T-1,Nh).*(params.Vmax-params.Vmin))+(~logic).*V;

V0=zeros(T,Nh);
V0(1,:)=params.Vini;
V0(2:end,:)=V;
V1=zeros(T,Nh);
V1(1:end-1,:)=V;
V1(end,:)=params.Vend;
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
    VV=zeros(T+1,1);
    VV(1)=params.Vini(i);
    VV(end)=params.Vend(i);
    
    for t=1+randperm(T-1)
        backwardV=VV;
        backwardV(1:t)=0;
        [tl,~,Vl]=find(backwardV,1);
        forwardV=VV;
        forwardV(t:end)=0;
        [tu,~,Vu]=find(forwardV,1,'last');
        
        
        Vl=max(params.Vmin(t-1,i),Vl-sum(inflow(t:tl(1)-1,i)-params.Qmin(t:tl(1)-1,i)));
        Vu=min(params.Vmax(t-1,i),Vu+sum(inflow(tu(1):t-1,i)-params.Qmin(tu(1):t-1,i)));
        if V(t-1,i)<Vl(1)||V(t-1,i)>Vu(1)
            V(t-1,i)=Vl(1)+rand*(Vu(1)-Vl(1));
            
            if V(t-1,i)>params.Vmax(t-1,i)
                V(t-1,i)=params.Vmax(t-1,i);
            elseif V(t-1,i)>params.Vmin(t-1,i)
            else
                V(t-1,i)=params.Vmin(t-1,i);
            end
        end
        VV(t)=V(t-1,i);
    end
    V0(:,i)=VV(1:end-1);
    V1(:,i)=VV(2:end);

    R(:,i)=V0(:,i)+inflow(:,i)-V1(:,i);
end
a=params.C2;
b=params.C3.*V1+params.C5;
c=params.C1.*(V1.^2)+params.C4.*V1+params.C6;

%note
%1.Ph can not reach params.Phmax with any V1
%2.There is a certain Q that makes Ph reach params.Phmin
%3.params.Qmin is less than extrpoint with any V1
discriminant=b.^2-4*a.*(c-params.Phmin);
extrpoint=-0.5*b./a;
Q1=0.5*(-b+sqrt(discriminant))./a;
Q2=0.5*(-b-sqrt(discriminant))./a;

Qmin=max(Q1,params.Qmin);
Qmax=min(min(Q2,params.Qmax),extrpoint);
Q=params.Qmin.*(R<Qmin)+R.*(R>=Qmin&R<=Qmax)+Qmax.*(R>Qmax);
% Q=R.*(R<=Qmax)+Qmax.*(R>Qmax);
if params.PDZ
    logic=Q<params.Qpu&Q>params.Qpl;
    Q=Q.*(~logic)+params.Qpl.*(logic);
end
SP=R-Q;
Ph=params.C1.*(V1.^2)+params.C2.*(Q.^2)+params.C3.*(V1.*Q)+...
    params.C4.*V1+params.C5.*Q+params.C6;
Ps=params.load-sum(Ph,2);
if params.VP
    cost=sum(params.a.*Ps.^2+params.b.*Ps+params.c+abs(params.e.*sin(params.f.*(params.Psmin-Ps))));
else
    cost=sum(params.a.*Ps.^2+params.b.*Ps+params.c);
end

VminViol=max(0,params.Vmin-V);
VmaxViol=max(0,V-params.Vmax);
% VendViol=max(0,abs(V(end,:)-params.Vend)-1e-5);
RminViol=max(0,Qmin-R);
% SPViol=max(0,0-SP);
PsminViol=max(0,params.Psmin-Ps);
PsmaxViol=max(0,Ps-params.Psmax);

viol=[VminViol(:);VmaxViol(:);RminViol(:);PsminViol(:);PsmaxViol(:)];
viol=sum(viol);
obvalue_viol=[viol,cost];
end

