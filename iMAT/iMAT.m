function Sol = iMAT(model,rxn_exp,coreset,eps,context)

%Modfied by Semidan, November 2014
%some lines are included to implement iMAT in gurobi, to allow defining the
%RH sert directly (to use the same core set used in CorEx) and to generate
%Cobra compatible structures

% cutoffs at lower and upper quantiles
%This part has been modified to directly find the reaction indeces through
%the quantile value of the data vector (instead of looking at reaction
%names). In addition, only data values over 0 are taken into account to
%calculate the quantile (This is necessary to accomodate the result of the
%function mapgene2rxn, which outputs 0 when no data is avaialable).

if isempty(coreset) && ~isempty(rxn_exp),
   RHindex=find(rxn_exp >= quantile(rxn_exp(rxn_exp>0),0.75));
   RLindex=setdiff(1:length(model.rxns),ExpressedRxns);
elseif ~isempty(coreset) && isempty(rxn_exp),
    RHindex=coreset;
    RLindex=setdiff(1:length(model.rxns),coreset);
end

clear toc
tic

% below copied from createTissueSpecificModel.m, part of the Cobra toolbox
% http://opencobra.sf.net/

S = model.S;
lb = model.lb;
ub = model.ub;

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

% Creating csense
sense1(1:size(S,1)) = '=';
sense2(1:length(RHindex)) = '>';
sense3(1:length(RHindex)) = '<';
sense4(1:length(RLindex)) = '>';
sense5(1:length(RLindex)) = '<';
sense = [sense1 sense2 sense3 sense4 sense5];

% Creating lb and ub
lb_y = zeros(2*length(RHindex)+length(RLindex),1);
ub_y = ones(2*length(RHindex)+length(RLindex),1);
lb = [lb;lb_y];
ub = [ub;ub_y];

% Creating c
c_v = zeros(size(S,2),1);
c_y = ones(2*length(RHindex)+length(RLindex),1);
c = [c_v;c_y];

% Creating b
b_s = zeros(size(S,1),1);
lb_rh = lb(RHindex);
ub_rh = ub(RHindex);
lb_rl = lb(RLindex);
ub_rl = ub(RLindex);
b = [b_s;lb_rh;ub_rh;lb_rl;ub_rl];

% Creating vartype
vartype1(1:size(S,2),1) = 'C';
vartype2(1:2*length(RHindex)+length(RLindex),1) = 'B';
vartype = [vartype1;vartype2];

MILPproblem.A = A;
MILPproblem.b = b;
MILPproblem.lb = lb;
MILPproblem.ub = ub;
MILPproblem.vartype = vartype;
MILPproblem.obj=c;
MILPproblem.sense=sense;
MILPproblem.rhs=b;
MILPproblem.modelsense='max';
params.OutputFlag=0;
params.Presolve=2; 
params.TimeLimit=60;
solution = gurobi(MILPproblem,params);
x = solution.x(1:size(S,2),1);

for i = 1:length(x)
    if abs(x(i)) < 1e-6
        x(i,1) = 0;
    end
end

%Added by Semidan, November 2014
v_sol = x;
v_info=zeros(6,1);
v_info(1)=length(find(abs(v_sol)>=eps));
v_info(2)=length(RHindex);
v_info(3)=length(RLindex);
v_info(4)=length(intersect(find(abs(v_sol)>=eps),find(rxn_exp>0)));
Corr=corrcoef(v_sol(find(rxn_exp>0)),rxn_exp(find(rxn_exp>0)));
v_info(5)=Corr(1,2);
v_info(6)=toc;

Sol.Flux=v_sol;
Sol.Info=v_info;
Sol.Zopt=round(sum(solution.x((size(S,2)+1):length(solution.x))));
Sol.Mod=contextmodel2COBRA(v_sol,model,context);

end