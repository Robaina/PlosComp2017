function Sol=AltNet(GEM,C,SolCorEx,Vmax,epsilon,Nmodels,time_limit,IterLimit)

%***************************The AltNet procedure***************************
%**************************************************************************

%This function implements the AltNet procedure to generate alternative
%optimal networks to CorEx and FastCORE.

%Required arguments are:
  %GEM: COBRA-like structure with the GEM
  
  %C: an array with the indexes in the GEM of the reactions in the core set
  
  %Z: the number of added non-core reactions, it must be equal or
  %greater than the optimal value found by a previous unconstrained CorEx
  %optimization.
  
  %Xopt: a binary array produced by the CorEx function. It contains
  %information about the non-core reactions that added at the optimum
  
  %A: an array generated by the CorEx function, it contains the indexes in
  %the GEM of the reactions in the optimal context-specific model.
  
%Optional values are:
  %Vmax: an upper bound for the maximum flux value of any reaction (in
  %reversibles the bound would be -Vmax <= v <= Vmax), default 1e3
  
  %epsilon: the numerical threshold to consider a reaction active (default
  %1e-3)
  
  %Nmodels: number of generated alternative optimal models
  
  %time-limit: a time limit (in seconds) for an early termination of the
  %MILP (default 60s)
  
  %IterLimit: an upper bound on the number of iterations in the main loop.
  %The loop continues until Nmodels is satisfied or the number of
  %iterations reaches IterLimit
  
%**************************************************************************
%           Semidán (robaina@mpimp-golm.mpg.de), May, 2016
%**************************************************************************

Z = SolCorEx.AddedRxns;
Xopt = SolCorEx.Xopt;
A = SolCorEx.Abinary;
RevRxns = find(GEM.rev==1);
IrrRxns = find(GEM.rev==0);
NRev = length(RevRxns);
NIrr = length(IrrRxns);
P = setdiff(1:length(GEM.rxns),C);
S = GEM.S;

Irr_C=intersect(C,IrrRxns);
Irr_P=intersect(P,IrrRxns);
Rev_C=intersect(C,RevRxns);
Rev_P=intersect(P,RevRxns);
NIrr_C=length(Irr_C);
NIrr_P=length(Irr_P);
NRev_C=length(Rev_C);
NRev_P=length(Rev_P);

%Stoichiometric Matrix reorganization
Sam=[S(:,Irr_C),S(:,Irr_P),S(:,Rev_C),S(:,Rev_P),-S(:,Rev_C),-S(:,Rev_P)];  
Rxns=size(Sam,2);
Mets=size(Sam,1);

%Construction of Amat matrix (virr_C,virr_P,vfor_C,vfor_P,vrev_C,vrev_P,xirr(P),xRev(P),y,delta+,delta-)
A1=[Sam,zeros(Mets,NIrr_P+NRev_P+NRev+2*(NIrr_P+NRev_P))]; %SV=0
A2=[eye(NIrr_C),zeros(NIrr_C,NIrr_P+2*NRev+NIrr_P+NRev_P+NRev+2*(NIrr_P+NRev_P))]; %virr_C>=epsilon 
A3=[zeros(NRev_C,NIrr),eye(NRev_C),zeros(NRev_C,NRev_P),eye(NRev_C),zeros(NRev_C,NRev_P+NIrr_P+NRev_P+NRev+2*(NIrr_P+NRev_P))]; %vfor_C+vrev_C>=epsilon

A4a=[zeros(NIrr_P,NIrr_C),eye(NIrr_P),zeros(NIrr_P,2*NRev),-Vmax*eye(NIrr_P),zeros(NIrr_P,NRev_P+NRev+2*(NIrr_P+NRev_P))]; %virr_p-xirr*Vmax<=0
A4b=[zeros(NIrr_P,NIrr_C),eye(NIrr_P),zeros(NIrr_P,2*NRev),-epsilon*eye(NIrr_P),zeros(NIrr_P,NRev_P+NRev+2*(NIrr_P+NRev_P))]; %virr_p-xirr*epsilon>=0

A5a=[zeros(NRev_P,NIrr+NRev_C),eye(NRev_P),zeros(NRev_P,NRev_C),eye(NRev_P),zeros(NRev_P,NIrr_P),-Vmax*eye(NRev_P),zeros(NRev_P,NRev+2*(NIrr_P+NRev_P))]; %(vfor_p + vback_p)-xRev*Vmax<=0
A5b=[zeros(NRev_P,NIrr+NRev_C),eye(NRev_P),zeros(NRev_P,NRev_C),eye(NRev_P),zeros(NRev_P,NIrr_P),-epsilon*eye(NRev_P),zeros(NRev_P,NRev+2*(NIrr_P+NRev_P))]; %(vfor_p + vback_p)-xRev*epsilon>=0

A7=[zeros(NRev,NIrr),eye(NRev),zeros(NRev,NRev+NIrr_P+NRev_P),Vmax*eye(NRev),zeros(NRev,+2*(NIrr_P+NRev_P))]; %vfor+y*Vmax<=Vmax
A8=[zeros(NRev,NIrr+NRev),eye(NRev),zeros(NRev,NIrr_P+NRev_P),-Vmax*eye(NRev),zeros(NRev,+2*(NIrr_P+NRev_P))]; %vrev-y*Vmax<=0
A9a=[zeros(1,Rxns),ones(1,NIrr_P+NRev_P),zeros(1,NRev+2*(NIrr_P+NRev_P))]; %||x||_1<Z+epsilon
A9b=[zeros(1,Rxns),ones(1,NIrr_P+NRev_P),zeros(1,NRev+2*(NIrr_P+NRev_P))]; %||x||_1>Z-epsilon

A10=[zeros(NIrr_P+NRev_P,Rxns),eye(NIrr_P+NRev_P),zeros(NIrr_P+NRev_P,NRev),eye(NIrr_P+NRev_P),-eye(NIrr_P+NRev_P)]; %delta+ - delta- +x =Xopt
A11=[zeros(NIrr_P+NRev_P,Rxns+NIrr_P+NRev_P+NRev),eye(NIrr_P+NRev_P),eye(NIrr_P+NRev_P)]; %delta+ + delta- <= 1

Amat=[A1;A2;A3;A4a;A4b;A5a;A5b;A7;A8;A9a;A9b;A10;A11];

%Construction of c, lb, ub, sense and b vectors
c=[zeros(Rxns+NIrr_P+NRev_P+NRev,1);ones(2*(NIrr_P+NRev_P),1)];
vsense=[repmat('=',Mets,1);repmat('>',NIrr_C+NRev_C,1);repmat('<',NIrr_P,1);repmat('>',NIrr_P,1);repmat('<',NRev_P,1);repmat('>',NRev_P,1);repmat('<',2*NRev,1);'<';'>';repmat('=',NIrr_P+NRev_P,1);repmat('<',NIrr_P+NRev_P,1)];
vtype=[repmat('C',Rxns,1);repmat('B',NIrr_P+NRev_P+NRev,1);repmat('B',2*(NIrr_P+NRev_P),1)];
lb=zeros(size(Amat,2),1);
ub=[Vmax*ones(Rxns,1);ones(NIrr_P+NRev_P+NRev+2*(NIrr_P+NRev_P),1)];

%Prepare gurobi structure
model.A=sparse(Amat);
model.modelsense='max';
model.obj=c;
model.sense=vsense;
model.lb=lb;
model.ub=ub;
model.vtype=vtype;
params.OutputFlag=1;
params.FeasibilityTol=1e-9;
params.IntFeasTol=1e-9;
params.TimeLimit=time_limit;
% params.presolve=2;

%Reconstruct vector x of binary variables
Abinary=zeros(length(GEM.rxns),1);
Abinary(A)=1;
% Xopt2=Abinary([intersect(P,IrrRxns);intersect(P,RevRxns)]);

%Start Iteration
Xopt2=Xopt;
Vmatrix=zeros(length(GEM.rxns),Nmodels);
Modmatrix=zeros(length(GEM.rxns),Nmodels+1);
QualityCheck=zeros(Nmodels,2);
MaxDifferences=zeros(Nmodels,1);
Xoptmatrix=zeros(length(Xopt2),Nmodels+1);
wbar=waitbar(0,'Getting Alternative Networks...');
netcounter=1;itercounter=1;
Modmatrix(:,1)=Abinary;
Xoptmatrix(:,1)=Xopt2;
tic
while netcounter<Nmodels && itercounter<IterLimit,
    waitbar(netcounter/Nmodels)
    %Solve MILP  
    model.rhs=[zeros(Mets,1);epsilon*ones(NIrr_C+NRev_C,1);zeros(2*NIrr_P+2*NRev_P,1);Vmax*ones(NRev,1);zeros(NRev,1);Z+0.5;Z-0.5;Xopt2;ones(NIrr_P+NRev_P,1)];
    model.start=[NaN*ones(Rxns,1);Xopt2;NaN*ones(NRev+2*(NIrr_P+NRev_P),1)];
    try
      gur=gurobi(model,params);
      Xoptmatrix(:,netcounter+1)=Xopt2;
      Xopt2=gur.x(Rxns+1:Rxns+NIrr_P+NRev_P); 
      Vmatrix(Irr_C,netcounter)=gur.x(1:NIrr_C);
      Vmatrix(Irr_P,netcounter)=gur.x(NIrr_C+1:NIrr);
      Vmatrix(Rev_C,netcounter)=gur.x(NIrr+1:NIrr+NRev_C)-gur.x(NIrr+NRev+1:NIrr+NRev+NRev_C);
      Vmatrix(Rev_P,netcounter)=gur.x(NIrr+NRev_C+1:NIrr+NRev_C+NRev_P)-gur.x(NIrr+NRev+NRev_C+1:NIrr+NRev+NRev_C+NRev_P);
      Vmatrix(abs(Vmatrix)<(epsilon-0.01*epsilon))=0;
      Modmatrix((abs(Vmatrix(:,netcounter))>0),netcounter+1)=1;
    catch
        Modmatrix(:,netcounter+1)=Modmatrix(:,netcounter);
    end
        MaxDifferences(netcounter)=sum(abs(Modmatrix(:,netcounter+1)-Modmatrix(:,netcounter)));
        netcounter=netcounter+1;

end

%Assesst quality of the networks: all networks must include the core set
%and be consistent (have non-zero flux only in all selected reactions)
for i=1:Nmodels,
    QualityCheck(i,1)=length(find(abs(Vmatrix(C,i))>0));
    QualityCheck(i,2)=length(find(abs(Vmatrix(:,i))>0));
end
Modmatrix(:,sum(Modmatrix)==0)=[];
QualityCheck=[{'Core Active Rxns','Total Active Rxns'};num2cell(QualityCheck)];
close(wbar)
Sol.Vmatrix=Vmatrix;
Sol.Modmatrix=unique(Modmatrix','rows')';
Sol.maxDiff=MaxDifferences;
Sol.QualityCheck=QualityCheck;
Sol.Xoptmatrix=Xoptmatrix;
Sol.TotalTime=toc;
Sol.Abinary=Abinary;

end
