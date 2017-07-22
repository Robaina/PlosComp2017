function Sol=CorExCORDA(GEM,HC,MC,NC,OT,Z_HC,Z_MC,Z_NC,Z_OT,Vmax,epsilon,time_limit)

%*******************The CorEx method(CORDA adaptation)*********************
%**************************************************************************

%This function implements the CorEx method to obtain a context-specific
%model from a GEM and a core set of reactions. This is the adpatation of
%CorEx to obtain a model that is also a CORDA optimal model. To this end,
%the core set of reactions is divided into the four categories used in
%CORDA (HC, MC, NC and OT).

%Required arguments are:
  %GEM: COBRA-like structure with the GEM
  
  %HC,MC,NC,OT: arrays with the indexes in the GEM of the reactions in the
  %HC, MC, NC and OT sets.
  
%Optional values are:
  %Z_HC, Z_MC, Z_NC, Z_OT: the number of added reaction in each one of the
  %four categories. These values should equal the ones obtained by CORDA at
  %the optimum.

  %Vmax: an upper bound for the maximum flux value of any reaction (in
  %reversibles the bound would be -Vmax <= v <= Vmax), default 1e3
  
  %epsilon: the numerical threshold to consider a reaction active (default
  %1e-3)
  
  %time-limit: a time limit (in seconds) for an early termination of the
  %MILP (default 60s)
  
%**************************************************************************
%             Semidán (robaina@mpimp-golm.mpg.de), May, 2016
%**************************************************************************


if nargin<6 || isempty(Z_HC),
    Z_HC=[];
end
if nargin<7 || isempty(Z_MC),
    Z_MC=[];
end
if nargin<8 || isempty(Z_NC),
    Z_NC=[];
end
if nargin<9 || isempty(Z_OT),
    Z_OT=[];
end
if nargin<10 || isempty(Vmax),
    Vmax=1e3;
end
if nargin<11 || isempty(epsilon),
    epsilon=0.001;
end
if nargin<12 || isempty(time_limit),
    time_limit=60;
end

RevRxns=find(GEM.rev==1);
IrrRxns=find(GEM.rev==0);
NRev=length(RevRxns);
NIrr=length(IrrRxns);
P=setdiff(1:length(GEM.rxns),HC);

%Classify reactions
Irr_HC=intersect(HC,IrrRxns);
Irr_MC=intersect(MC,IrrRxns);
Irr_NC=intersect(NC,IrrRxns);
Irr_OT=intersect(OT,IrrRxns);
Irr_P=intersect(P,IrrRxns);
Rev_HC=intersect(HC,RevRxns);
Rev_MC=intersect(MC,RevRxns);
Rev_NC=intersect(NC,RevRxns);
Rev_OT=intersect(OT,RevRxns);
Rev_P=intersect(P,RevRxns);
NIrr_HC=length(Irr_HC);
NIrr_MC=length(Irr_MC);
NIrr_NC=length(Irr_NC);
NIrr_OT=length(Irr_OT);
NIrr_P=length(Irr_P);
NRev_HC=length(Rev_HC);
NRev_MC=length(Rev_MC);
NRev_NC=length(Rev_NC);
NRev_OT=length(Rev_OT);
NRev_P=length(Rev_P);

%Stoichiometric Matrix reorganization
S=GEM.S;
Sam=[S(:,Irr_HC),S(:,[Irr_MC;Irr_NC;Irr_OT]),S(:,Rev_HC),S(:,[Rev_MC;Rev_NC;Rev_OT]),-S(:,Rev_HC),-S(:,[Rev_MC;Rev_NC;Rev_OT])];  
Rxns=size(Sam,2);
Mets=size(Sam,1);

%Construction of Amat matrix (virr_C,virr_P,vfor_C,vfor_P,vrev_C,vrev_P,xirr(P),xRev(P),y)
A1=[Sam,zeros(Mets,NIrr_P+NRev_P+NRev)]; %SV=0
A2=[eye(NIrr_HC),zeros(NIrr_HC,NIrr_P+2*NRev+NIrr_P+NRev_P+NRev)]; %virr_C>=epsilon 
A3=[zeros(NRev_HC,NIrr),eye(NRev_HC),zeros(NRev_HC,NRev_P),eye(NRev_HC),zeros(NRev_HC,NRev_P+NIrr_P+NRev_P+NRev)]; %vfor_C+vrev_C>=epsilon

A4a=[zeros(NIrr_P,NIrr_HC),eye(NIrr_P),zeros(NIrr_P,2*NRev),-Vmax*eye(NIrr_P),zeros(NIrr_P,NRev_P+NRev)]; %virr_p-xirr*Vmax<=0
A4b=[zeros(NIrr_P,NIrr_HC),eye(NIrr_P),zeros(NIrr_P,2*NRev),-epsilon*eye(NIrr_P),zeros(NIrr_P,NRev_P+NRev)]; %virr_p-xirr*epsilon>=0

A5a=[zeros(NRev_P,NIrr+NRev_HC),eye(NRev_P),zeros(NRev_P,NRev_HC),eye(NRev_P),zeros(NRev_P,NIrr_P),-Vmax*eye(NRev_P),zeros(NRev_P,NRev)]; %(vfor_p+vback_p)-xRev*Vmax<=0
A5b=[zeros(NRev_P,NIrr+NRev_HC),eye(NRev_P),zeros(NRev_P,NRev_HC),eye(NRev_P),zeros(NRev_P,NIrr_P),-epsilon*eye(NRev_P),zeros(NRev_P,NRev)]; %(vfor_p+vback_p)-xRev*epsilon>=0

A7=[zeros(NRev,NIrr),eye(NRev),zeros(NRev,NRev+NIrr_P+NRev_P),Vmax*eye(NRev)]; %vfor+y*Vmax<=Vmax
A8=[zeros(NRev,NIrr+NRev),eye(NRev),zeros(NRev,NIrr_P+NRev_P),-Vmax*eye(NRev)]; %vrev-y*Vmax<=0
Amat=[A1;A2;A3;A4a;A4b;A5a;A5b;A7;A8];

%Construction of c, lb, ub, sense and b vectors
b=[zeros(Mets,1);epsilon*ones(NIrr_HC+NRev_HC,1);zeros(2*NIrr_P+2*NRev_P,1);Vmax*ones(NRev,1);zeros(NRev,1)];
vsense=[repmat('=',Mets,1);repmat('>',NIrr_HC+NRev_HC,1);repmat('<',NIrr_P,1);repmat('>',NIrr_P,1);repmat('<',NRev_P,1);repmat('>',NRev_P,1);repmat('<',2*NRev,1)];
if ~isempty(Z_HC),
    %These constraints are modified to account for the CORDA reaction
    %categories
    
   A9a=[zeros(1,Rxns),ones(1,NIrr_MC),zeros(1,NIrr_NC+NIrr_OT),ones(1,NRev_MC),zeros(1,NRev_NC+NRev_OT+NRev)]; %||x_mc||_1<Z_MC+epsilon
   A9b=[zeros(1,Rxns),ones(1,NIrr_MC),zeros(1,NIrr_NC+NIrr_OT),ones(1,NRev_MC),zeros(1,NRev_NC+NRev_OT+NRev)]; %||x_mc||_1>Z_MC-epsilon

   A9c=[zeros(1,Rxns+NIrr_MC),ones(1,NIrr_NC),zeros(1,NIrr_OT+NRev_MC),ones(1,NRev_NC),zeros(1,NRev_OT+NRev)]; %||x_nc||_1<Z_NC+epsilon
   A9d=[zeros(1,Rxns+NIrr_MC),ones(1,NIrr_NC),zeros(1,NIrr_OT+NRev_MC),ones(1,NRev_NC),zeros(1,NRev_OT+NRev)]; %||x_nc||_1>Z_NC-epsilon

   A9e=[zeros(1,Rxns+NIrr_MC+NIrr_NC),ones(1,NIrr_OT),zeros(1,NRev_MC+NRev_NC),ones(1,NRev_OT),zeros(1,NRev)]; %||x_ot||_1<Z_OT+epsilon
   A9f=[zeros(1,Rxns+NIrr_MC+NIrr_NC),ones(1,NIrr_OT),zeros(1,NRev_MC+NRev_NC),ones(1,NRev_OT),zeros(1,NRev)]; %||x_ot||_1>Z_OT-epsilon
   
   Amat=[Amat;A9a;A9b;A9c;A9d;A9e;A9f];
   b=[b;Z_MC+0.5;Z_MC-0.5;Z_NC+0.5;Z_NC-0.5;Z_OT+0.5;Z_OT-0.5];
   vsense=[vsense;'<';'>';'<';'>';'<';'>'];
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
params.OutputFlag=1;
params.FeasibilityTol=1e-9;
params.IntFeasTol=1e-9;
params.TimeLimit=time_limit;
gur=gurobi(model,params);

V=zeros(length(GEM.rxns),1);
V(Irr_HC)=gur.x(1:NIrr_HC);
V(Irr_MC)=gur.x((NIrr_HC+1):(NIrr_HC+NIrr_MC));
V(Irr_NC)=gur.x((NIrr_HC+NIrr_MC+1):(NIrr_HC+NIrr_MC+NIrr_NC));
V(Irr_OT)=gur.x((NIrr_HC+NIrr_MC+NIrr_NC+1):(NIrr_HC+NIrr_MC+NIrr_NC+NIrr_OT));
V(Rev_HC)=gur.x(NIrr+1:NIrr+NRev_HC)-gur.x(NIrr+NRev+1:NIrr+NRev+NRev_HC);
V(Rev_MC)=gur.x(NIrr+NRev_HC+1:NIrr+NRev_HC+NRev_MC)-gur.x(NIrr+NRev+NRev_HC+1:NIrr+NRev+NRev_HC+NRev_MC);
V(Rev_NC)=gur.x(NIrr+NRev_HC+NRev_MC+1:NIrr+NRev_HC+NRev_MC+NRev_NC)-gur.x(NIrr+NRev+NRev_HC+NRev_MC+1:NIrr+NRev+NRev_HC+NRev_MC+NRev_NC);
V(Rev_OT)=gur.x(NIrr+NRev_HC+NRev_MC+NRev_NC+1:NIrr+NRev_HC+NRev_MC+NRev_NC+NRev_OT)-gur.x(NIrr+NRev+NRev_HC+NRev_MC+NRev_NC+1:NIrr+NRev+NRev_HC+NRev_MC+NRev_NC+NRev_OT);
V(abs(V)<(epsilon-0.01*epsilon))=0;
A=zeros(length(V),1);A(abs(V)>0)=1;

QualityCheck(1,1)=length(find(abs(V(HC))>0));
QualityCheck(1,2)=length(find(abs(V(MC))>0));
QualityCheck(1,3)=length(find(abs(V(NC))>0));
QualityCheck(1,4)=length(find(abs(V(OT))>0));
QualityCheck(1,5)=length(find(abs(V(:))>0));

Sol.V=V;
Sol.QualityCheck=[{'HC Active Rxns','MC Active Rxns','NC Active Rxns','OT Active Rxns','Total Active Rxns'};num2cell(QualityCheck)];
Sol.addedMC=sum(gur.x([zeros(1,Rxns),ones(1,NIrr_MC),zeros(1,NIrr_NC+NIrr_OT),ones(1,NRev_MC),zeros(1,NRev_NC+NRev_OT+NRev)]==1));
Sol.addedNC=sum(gur.x([zeros(1,Rxns+NIrr_MC),ones(1,NIrr_NC),zeros(1,NIrr_OT+NRev_MC),ones(1,NRev_NC),zeros(1,NRev_OT+NRev)]==1));
Sol.addedOT=sum(gur.x([zeros(1,Rxns+NIrr_MC+NIrr_NC),ones(1,NIrr_OT),zeros(1,NRev_MC+NRev_NC),ones(1,NRev_OT),zeros(1,NRev)]==1));
Sol.A=A;
% Sol.AddedRxns=find(X==1);
Sol.Xopt=gur.x(Rxns+1:Rxns+NIrr_P+NRev_P);
Sol.Abinary=A;

end