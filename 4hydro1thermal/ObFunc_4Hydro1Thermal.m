function  [obvalue_viol,V1,R,Q,Ph,SP,Ps]= ObFunc_4Hydro1Thermal(Q)  %#codegen
global params;
inflow=params.I;
[T,Nh]=size(inflow);

V0=zeros(T,Nh);
V0(1,:)=params.Vini;

V1=zeros(T,Nh);
V1(end,:)=params.Vend;
R=Q;
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
%     m=1;
%     t=randperm(T);
%     while m<T+1
%         VT=params.Vini(i)+sum(inflow(:,i)-R(:,i));
%         deltaV=VT-params.Vend(i);
%         if deltaV<-eps
%             R(t(m),i)=max(R(t(m),i)+deltaV,params.Qmin(t(m),i));
%         elseif deltaV>eps            
%             R(t(m),i)=R(t(m),i)+deltaV;    
%         else
%             break;
%         end
%         m=m+1;
%     end
%     for t=1:T-1
%         V1(t,i)=V0(t,i)+inflow(t,i)-R(t,i);
%         V0(t+1,i)=V1(t,i);
%     end


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
%
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

VminViol=max(0,params.Vmin-V1(1:end-1,:)-1e-5);
VmaxViol=max(0,V1(1:end-1,:)-params.Vmax-1e-5);
VendViol=max(0,abs(V1(end,:)-params.Vend)-1e-5);
RminViol=max(0,params.Qmin-R-1e-5);
% SPViol=max(0,0-SP);
PhminViol=max(0,params.Phmin-Ph-1e-5);
PhmaxViol=max(0,Ph-params.Phmax-1e-5);
PsminViol=max(0,params.Psmin-Ps-1e-5);
PsmaxViol=max(0,Ps-params.Psmax-1e-5);

viol=[VminViol(:);VmaxViol(:);VendViol(:);RminViol(:);PhminViol(:);PhmaxViol(:);PsminViol(:);PsmaxViol(:)];
viol=sum(viol);
obvalue_viol=[viol,cost];
end

