function [ObjVal,FluxValInac,FluxValAc] = iMATFVA(model,rxn_exp,coreset,epsilon)

%*********************iMAT Flux Variability Analysis***********************
%**************************************************************************

%This function implements the iMAT FVA procedure to evaluate the
%alternative optima space of the iMAT method. The main iMAT function has
%been adapted from part of the code in createTissueSpecificModel.m, in the
%Cobra toolbox(http://opencobra.sf.net/).This code has been altered to
%adapt the optimization program to the new one developed in this study, as
%well as adapt its structure to the gurobi solver. Also the possibility of
%using a predefined core set as the RH group of reactions instead of
%expression data has been addded.

%Required arguments are:
  %model: COBRA-like structure with the GEM
  
  %Zopt: the optimum value of a previous iMAT solution (objective value)
  
  %rxn_exp or coreset: either an array containing the gene expression
  %values mapped to every reaction in the GEM or an array containing the
  %indexes in the GEM of the reactions that should be included in the RH
  %set
 
%Optional values are:
  %epsilon: the numerical threshold to consider a reaction active (default
  %1e-6)

%**************************************************************************
%         Semidán (robaina@mpimp-golm.mpg.de), May, 2016
%**************************************************************************

clear toc
tic
if isempty(coreset) && ~isempty(rxn_exp),
   RHindex=find(rxn_exp >= quantile(rxn_exp(rxn_exp>0),0.75));
   RLindex=setdiff(1:length(model.rxns),ExpressedRxns);
elseif ~isempty(coreset) && isempty(rxn_exp),
    RHindex=coreset;
    RLindex=setdiff(1:length(model.rxns),coreset);
end

S = model.S;
FluxValInac = zeros(length(model.rxns));
FluxValAc = FluxValInac;

%**************************************************************************
%Evaluate the iMAT objective function upon inactivation of each reaction
%**************************************************************************
for nk=1:length(model.rxns),
    lb = model.lb;
    ub = model.ub;
    lb(nk)=0;
    ub(nk)=0;
    % Creating A matrix
    A = sparse(size(S,1)+2*length(RHindex)+2*length(RLindex),size(S,2)+2*length(RHindex)+length(RLindex));
    [m,n,s] = find(S);
    for i = 1:length(m)
        A(m(i),n(i)) = s(i); %#ok<SPRIX>
    end

    for i = 1:length(RHindex)
        A(i+size(S,1),RHindex(i)) = 1; %#ok<SPRIX>
        A(i+size(S,1),i+size(S,2)) = lb(RHindex(i)) - epsilon; %#ok<SPRIX>
        A(i+size(S,1)+length(RHindex),RHindex(i)) = 1; %#ok<SPRIX>
        A(i+size(S,1)+length(RHindex),i+size(S,2)+length(RHindex)+length(RLindex)) = ub(RHindex(i)) + epsilon; %#ok<SPRIX>
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
    try
       ObjVal(nk,1) = solution.objval;
       FluxValInac(:,nk) = solution.x(1:size(S,2));
    catch
        ObjVal(nk,1) = nan;
    end
end
%**************************************************************************
%Evaluate the iMAT objective function upon activation of each reaction
%**************************************************************************
for nk=1:length(model.rxns),
    lb = model.lb;
    ub = model.ub;

    % Creating A matrix
    A = sparse(size(S,1)+2*length(RHindex)+2*length(RLindex),size(S,2)+2*length(RHindex)+length(RLindex));
    [m,n,s] = find(S);
    for i = 1:length(m)
        A(m(i),n(i)) = s(i); %#ok<SPRIX>
    end
    
    %Add constraint on y variables to force activation of reaction "nk"
    %(Added by Semidan,May 2016)
    B = zeros(1,size(A,2));
    if ismember(nk,RHindex),
        B([size(S,2)+find(RHindex==nk),size(S,2)+length(RHindex)+find(RHindex==nk)])=1;
%         csense6 = 'G';
        sense6 = '>';
    elseif ismember(nk,RLindex),
        B(size(S,2)+2*length(RHindex)+find(RLindex==nk))=1;
        lb(nk)=epsilon;
%         csense6 = 'E';
        sense6 = '=';
    end
    
    for i = 1:length(RHindex)
        A(i+size(S,1),RHindex(i)) = 1; %#ok<SPRIX>
        A(i+size(S,1),i+size(S,2)) = lb(RHindex(i)) - epsilon; %#ok<SPRIX>
        A(i+size(S,1)+length(RHindex),RHindex(i)) = 1; %#ok<SPRIX>
        A(i+size(S,1)+length(RHindex),i+size(S,2)+length(RHindex)+length(RLindex)) = ub(RHindex(i)) + epsilon; %#ok<SPRIX>
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
    sense = [sense1 sense2 sense3 sense4 sense5 sense6];

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
    b = [b_s;lb_rh;ub_rh;lb_rl;ub_rl;0];

    % Creating vartype
    vartype1(1:size(S,2),1) = 'C';
    vartype2(1:2*length(RHindex)+length(RLindex),1) = 'B';
    vartype = [vartype1;vartype2];

    MILPproblem.A = [A;B];
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
    try
       ObjVal(nk,2) = solution.objval;
       FluxValAc(:,nk) = solution.x(1:size(S,2));
    catch
        ObjVal(nk,2) = nan;
    end  
end
ObjVal = round(ObjVal);
end