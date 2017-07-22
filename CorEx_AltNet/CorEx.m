function Sol=CorEx(GEM,C,Z,Vmax,epsilon,time_limit)

%***************************The CorEx method*******************************
%**************************************************************************

%This function implements the CorEx method to obtain a context-specific
%model from a GEM and a core set of reactions.

%Required arguments are:
  %GEM: COBRA-like structure with the GEM
  
  %C: an array with the indexes in the GEM of the reactions in the core set
  
%Optional values are:
  %Z: the number of added non-core reactions, if added it must be equal or
  %greater than the optimal value found by a previous unconstrained CorEx
  %optimization.
  
  %Vmax: an upper bound for the maximum flux value of any reaction (in
  %reversibles the bound would be -Vmax <= v <= Vmax), default 1e3
  
  %epsilon: the numerical threshold to consider a reaction active (default
  %1e-3)
  
  %time-limit: a time limit (in seconds) for an early termination of the
  %MILP (default 60s)
  
%**************************************************************************
%           Semidán (robaina@mpimp-golm.mpg.de), May, 2016
%**************************************************************************

if nargin<3 || isempty(Z),
    Z=[];
end
if nargin<4 || isempty(Vmax),
    Vmax=1e3;
end
if nargin<5 || isempty(epsilon),
    epsilon=0.001;
end
if nargin<6 || isempty(time_limit),
    time_limit=60;
end

RevRxns=find(GEM.rev==1);
IrrRxns=find(GEM.rev==0);
NRev=length(RevRxns);
NIrr=length(IrrRxns);
P=setdiff(1:length(GEM.rxns),C);

%Stoichiometric Matrix reorganization
Irr_C=intersect(C,IrrRxns);
Irr_P=intersect(P,IrrRxns);
Rev_C=intersect(C,RevRxns);
Rev_P=intersect(P,RevRxns);
NIrr_C=length(Irr_C);
NIrr_P=length(Irr_P);
NRev_C=length(Rev_C);
NRev_P=length(Rev_P);
S=GEM.S;
Sam=[S(:,Irr_C),S(:,Irr_P),S(:,Rev_C),S(:,Rev_P),-S(:,Rev_C),-S(:,Rev_P)];  
Rxns=size(Sam,2);
Mets=size(Sam,1);

%Construction of Amat matrix (virr_C,virr_P,vfor_C,vfor_P,vrev_C,vrev_P,xirr(P),xRev(P),y)
A1=[Sam,zeros(Mets,NIrr_P+NRev_P+NRev)]; %SV=0
A2=[eye(NIrr_C),zeros(NIrr_C,NIrr_P+2*NRev+NIrr_P+NRev_P+NRev)]; %virr_C>=epsilon 
A3=[zeros(NRev_C,NIrr),eye(NRev_C),zeros(NRev_C,NRev_P),eye(NRev_C),zeros(NRev_C,NRev_P+NIrr_P+NRev_P+NRev)]; %vfor_C+vrev_C>=epsilon

A4a=[zeros(NIrr_P,NIrr_C),eye(NIrr_P),zeros(NIrr_P,2*NRev),-Vmax*eye(NIrr_P),zeros(NIrr_P,NRev_P+NRev)]; %virr_p-xirr*Vmax<=0
A4b=[zeros(NIrr_P,NIrr_C),eye(NIrr_P),zeros(NIrr_P,2*NRev),-epsilon*eye(NIrr_P),zeros(NIrr_P,NRev_P+NRev)]; %virr_p-xirr*epsilon>=0

A5a=[zeros(NRev_P,NIrr+NRev_C),eye(NRev_P),zeros(NRev_P,NRev_C),eye(NRev_P),zeros(NRev_P,NIrr_P),-Vmax*eye(NRev_P),zeros(NRev_P,NRev)]; %(vfor_p + vback_p)-xRev*Vmax<=0
A5b=[zeros(NRev_P,NIrr+NRev_C),eye(NRev_P),zeros(NRev_P,NRev_C),eye(NRev_P),zeros(NRev_P,NIrr_P),-epsilon*eye(NRev_P),zeros(NRev_P,NRev)]; %(vfor_p + vback_p)-xRev*epsilon>=0

A6=[zeros(NRev,NIrr),eye(NRev),zeros(NRev,NRev+NIrr_P+NRev_P),Vmax*eye(NRev)]; %vfor+y*Vmax<=Vmax
A7=[zeros(NRev,NIrr+NRev),eye(NRev),zeros(NRev,NIrr_P+NRev_P),-Vmax*eye(NRev)]; %vrev-y*Vmax<=0
Amat=[A1;A2;A3;A4a;A4b;A5a;A5b;A6;A7];

%Construction of c, lb, ub, sense and b vectors
b=[zeros(Mets,1);epsilon*ones(NIrr_C+NRev_C,1);zeros(2*NIrr_P+2*NRev_P,1);Vmax*ones(NRev,1);zeros(NRev,1)];
vsense=[repmat('=',Mets,1);repmat('>',NIrr_C+NRev_C,1);repmat('<',NIrr_P,1);repmat('>',NIrr_P,1);repmat('<',NRev_P,1);repmat('>',NRev_P,1);repmat('<',2*NRev,1)];
if ~isempty(Z),
   A9a=[zeros(1,Rxns),ones(1,NIrr_P+NRev_P),zeros(1,NRev)]; %||x||_1<Z+epsilon
   A9b=[zeros(1,Rxns),ones(1,NIrr_P+NRev_P),zeros(1,NRev)]; %||x||_1>Z-epsilon
   Amat=[Amat;A9a;A9b];
   b=[b;Z+0.5;Z-0.5];
   vsense=[vsense;'<';'>'];
end
c=[zeros(Rxns,1);ones(NIrr_P+NRev_P,1);zeros(NRev,1)];
vtype=[repmat('C',Rxns,1);repmat('B',NIrr_P+NRev_P+NRev,1)];
lb=zeros(size(Amat,2),1);
ub=[Vmax*ones(Rxns,1);ones(NIrr_P+NRev_P+NRev,1)];

%Solve MILP
model.A=sparse(Amat);
model.modelsense='min';
model.obj=c;
model.sense=vsense;
model.rhs=b;
model.lb=lb;
model.ub=ub;
model.vtype=vtype;
params.OutputFlag=0;
params.FeasibilityTol=1e-9;
params.IntFeasTol=1e-9;
params.TimeLimit=time_limit;
gur=gurobi(model,params);

V=zeros(length(GEM.rxns),1);
V(Irr_C)=gur.x(1:NIrr_C);
V(Irr_P)=gur.x(NIrr_C+1:NIrr);
V(Rev_C)=gur.x(NIrr+1:NIrr+NRev_C)-gur.x(NIrr+NRev+1:NIrr+NRev+NRev_C);
V(Rev_P)=gur.x(NIrr+NRev_C+1:NIrr+NRev_C+NRev_P)-gur.x(NIrr+NRev+NRev_C+1:NIrr+NRev+NRev_C+NRev_P);
 
% X(Irr_P)=gur.x(Rxns+1:Rxns+NIrr_P);
% X(Rev_P)=gur.x(Rxns+NIrr_P+1:Rxns+NIrr_P+NRev_P);
% X=round(X);
V(abs(V)<(epsilon-0.01*epsilon))=0;
A=zeros(length(V),1);A(abs(V)>0)=1;

QualityCheck(1,1)=length(find(abs(V(C))>0));
QualityCheck(1,2)=length(find(abs(V(P))>0));
QualityCheck(1,3)=length(find(abs(V)>0));
Sol.QualityCheck=[{'C Active Rxns','P Active Rxns','Total Active Rxns'};num2cell(QualityCheck)];
Sol.V=V;
Sol.Vcore=V(C);
% Sol.Vadded=V((X==1));
% Sol.Vnoadded=V(intersect(P,find(X==0)));
Sol.AddedRxns=length(find(abs(V)>0))-length(find(abs(V(C))>0));
Sol.Xopt=gur.x(Rxns+1:Rxns+NIrr_P+NRev_P);
Sol.Abinary=A;

end