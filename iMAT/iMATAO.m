function Sol = iMATAO(model,Zopt,rxn_exp,coreset,eps,samplesize,FVArange,deltamax)

%*********************iMAT Alternative Optima Sampling*********************
%**************************************************************************

%This function implements the iMAT alternative optima sampling method. The
%original iMAT implementation has been copied from the available code
%inside the COBRA toolbox:
%createTissueSpecificModel.m (http://opencobra.sf.net/).
%This code has been altered to adapt the optimization program to the new
%one developed in this study, as well as adapt its structure to the gurobi
%solver. Also the possibility of using a predefined core set as the RH
%group of reactions instead of expression data has been addded.

%Required arguments are:
  %model: COBRA-like structure with the GEM
  
  %Zopt: the optimum value of a previous iMAT solution (objective value)
  
  %rxn_exp or coreset: either an array containing the gene expression
  %values mapped to every reaction in the GEM or an array containing the
  %indexes in the GEM of the reactions that should be included in the RH
  %set
 
%Optional values are:
  %eps: the numerical threshold to consider a reaction active (default
  %1e-6)
  
  %samplesize: the number of sampled alternative optimal solutions (default
  %1e3)
  
  %FVArange: the min and max allowable flux values per each reaction, if
  %empty they are calculated through the function "FVA"
  
  %deltamax: maximum value of the error term in the objective function
  %(default 10)
  
%**************************************************************************
%         Semidán (robaina@mpimp-golm.mpg.de), May, 2016
%**************************************************************************
  
if nargin<3 || isempty(rxn_exp),
    rxn_exp=[];
end
if nargin<4 || isempty(coreset),
    coreset=[];
end
if nargin<5 || isempty(eps),
    eps = 1e-6;
end
if nargin<6 || isempty(samplesize),
    samplesize = 1e3;
end
if nargin<7 || isempty(FVArange),
    [FVArange(:,1),FVArange(:,2)]=FVA(GEM,1,'gurobi','F');
end
if nargin<8 || isempty(deltamax),
    deltamax=10;
end
if isempty(coreset) && ~isempty(rxn_exp),
   RHindex=find(rxn_exp >= quantile(rxn_exp(rxn_exp>0),0.75));
   RLindex=setdiff(1:length(model.rxns),ExpressedRxns);
elseif ~isempty(coreset) && isempty(rxn_exp),
    RHindex=coreset;
    RLindex=setdiff(1:length(model.rxns),coreset);
end

% below copied from createTissueSpecificModel.m, part of the Cobra toolbox
% http://opencobra.sf.net/

S = model.S;
lb = model.lb;
ub = model.ub;
v_sol = zeros(size(S,2),samplesize);

% Creating A matrix
A = sparse(size(S,1)+2*length(RHindex)+2*length(RLindex),size(S,2)+2*length(RHindex)+length(RLindex));
[m,n,s] = find(S);
for i = 1:length(m)
    A(m(i),n(i)) = s(i); %#ok<SPRIX>
end

for i = 1:length(RHindex)
    A(i+size(S,1),RHindex(i)) = 1; %#ok<SPRIX>
    A(i+size(S,1),i+size(S,2)) = lb(RHindex(i)) - eps; %#ok<SPRIX>
    A(i+size(S,1)+length(RHindex),RHindex(i)) = 1; %#ok<SPRIX>
    A(i+size(S,1)+length(RHindex),i+size(S,2)+length(RHindex)+length(RLindex)) = ub(RHindex(i)) + eps; %#ok<SPRIX>
end

for i = 1:length(RLindex)
    A(i+size(S,1)+2*length(RHindex),RLindex(i)) = 1; %#ok<SPRIX>
    A(i+size(S,1)+2*length(RHindex),i+size(S,2)+length(RHindex)) = lb(RLindex(i)); %#ok<SPRIX>
    A(i+size(S,1)+2*length(RHindex)+length(RLindex),RLindex(i)) = 1; %#ok<SPRIX>
    A(i+size(S,1)+2*length(RHindex)+length(RLindex),i+size(S,2)+length(RHindex)) = ub(RLindex(i)); %#ok<SPRIX>
end
% Add set of new constraints to force v=vrand+delta and v \in AO of iMAT
B=[eye(size(S,2)),zeros(size(S,2),2*length(RHindex)+length(RLindex)),-eye(size(S,2))]; %v-delta = vrand
C1=[zeros(1,size(S,2)),ones(1,2*length(RHindex)+length(RLindex)),zeros(1,size(S,2))]; %sum(yplus + yminus + yplus(RL)) < Zopt + epsilon
C2=[zeros(1,size(S,2)),ones(1,2*length(RHindex)+length(RLindex)),zeros(1,size(S,2))]; %sum(yplus + yminus + yplus(RL)) > Zopt -epsilon
A=[[A,zeros(size(A,1),size(S,2))];B;C1;C2];A=sparse(A);

% Creating csense
sense1(1:size(S,1)) = '=';
sense2(1:length(RHindex)) = '>';
sense3(1:length(RHindex)) = '<';
sense4(1:length(RLindex)) = '>';
sense5(1:length(RLindex)) = '<';
sense6(1:(size(S,2))) = '=';
sense = [sense1 sense2 sense3 sense4 sense5 sense6 '<' '>'];

% Creating lb and ub
lb_y = zeros(2*length(RHindex)+length(RLindex),1);
ub_y = ones(2*length(RHindex)+length(RLindex),1);
lb = [lb;lb_y];
ub = [ub;ub_y];

% Creating c
c = zeros(size(A,2),1);

% Creating Q matrix
Q = sparse([zeros(2*size(S,2)+2*length(RHindex)+length(RLindex),size(S,2)+2*length(RHindex)+length(RLindex)),[zeros(size(S,2)+2*length(RHindex)+length(RLindex),size(S,2));eye(size(S,2))]]);

% Creating b
b_s = zeros(size(S,1),1);
lb_rh = lb(RHindex);
ub_rh = ub(RHindex);
lb_rl = lb(RLindex);
ub_rl = ub(RLindex);

% Creating vartype
vartype1(1:size(S,2),1) = 'C';
vartype2(1:2*length(RHindex)+length(RLindex),1) = 'B';
vartype3(1:size(S,2),1) = 'C';
vartype = [vartype1;vartype2;vartype3];

MIQPproblem.A = A;
MIQPproblem.Q = Q;
MIQPproblem.lb = [lb;-deltamax*ones(size(S,2),1)];
MIQPproblem.ub = [ub;deltamax*ones(size(S,2),1)];
MIQPproblem.vartype = vartype;
MIQPproblem.obj=c;
MIQPproblem.sense=sense;
MIQPproblem.modelsense='min';
params.OutputFlag=0;
params.Presolve=2; 
params.TimeLimit=60;

%Generate the random sample of AO flux distributions
i=1;
while i<=samplesize,
    %create random flux distribution
    vrand = (FVArange(:,2)-FVArange(:,1)).*rand(size(S,2),1)+FVArange(:,1);
    MIQPproblem.rhs=[b_s;lb_rh;ub_rh;lb_rl;ub_rl;vrand;Zopt+0.1;Zopt-0.1];
    try
       solution = gurobi(MIQPproblem,params);
       v_sol(:,i) = solution.x(1:size(S,2),1);
       Z(i,1) = round(sum(solution.x((size(S,2)+1):(size(S,2)+2*length(RHindex)+length(RLindex)))));
       i=i+1;
    catch
       if params.OutputFlag==1,
            sprintf('Failed attempt to solve QP in n=%d',i);
       end
    end
end
v_sol(v_sol<eps)=0;
Modmatrix=zeros(size(v_sol));
Modmatrix(abs(v_sol)>0)=1;
Sol.Modmatrix=unique(Modmatrix','rows')';
Sol.Vsample=v_sol;
Sol.Time=toc;
Sol.Z=Z;


end