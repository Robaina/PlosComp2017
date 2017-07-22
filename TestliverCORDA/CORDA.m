%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [tissue_model, rescue, HCtoMC, HCtoNC, MCtoNC] = CORDA(model,metTests,...
%   ES,PR,NP,PRtoNP,constraint,constrainby,om,ntimes,nl)
% CORDA calculates tissue specific reconstruction by determining reaction
% dependency using a cost association. If you wish to use default
% parameters (when possible) set them as empty
% INPUTS:
%   model - general model from which the tissue specific model will be
%       calculated
%   metTests - metabolic tests to be included in the reconstruction. This
%       argument should be a cell array of strings of size nx2, where n is 
%       the number of metabolic tests to be performed. Column 1 should be 
%       the name of the reaction to be included and the corresponding column 
%       2 should be the reaction to be included. For example, to test for
%       the production of pep and pyruvate, metTests should be equal to
%       {'DM_pep[c]' 'pep[c] -> ';'DM_pyr[c]' 'pyr[c] -> '}. May be left
%       empty.
%   ES - High confidence reactions. Reactions to be included in the model. Cell
%       array of strings.
%   PR - Medium confidence reactions to be included in the model if they do 
%       not depend on too many NP reactions to carry a flux. Cell array of 
%       strings.
%   NP - Negatice confidence reactions not to be included in the model. 
%       These reactions will be included in the tissue model only if they 
%       are necessary for the flux of ES reactions or for the flux of PRtoNP 
%       or more PR reactions. Cell array of strings.
% OPTIONAL INPUTS
%   PRtoNP - Define the threshold to include NP reactions based on PR
%       dependency. NP reactions will be included in the tissue model only
%       if they are deemed relevant for the flux of PRtoNP or more PR
%       reactions (or necessary for ES reactions). Default 2.
%   constraint - constraint value. Numerical. Default 1. Negative
%       percentage values will be used as positive.
%   constrainby - type of constraint used when defining reaction dependency. 
%       String. Either 'perc' or 'val'. 'perc' constrains reactions based 
%       on percentage from their optimal value. The reactions is first 
%       optimized and then held at the percentage from optima defined by 
%       the constraint value. 'val' constrains the reactions based on 
%       numerical values. Forward fluxes will be held numerically at the 
%       constraint value and negative fluxes at -constraint. Default 'val'
%   om - cost assigned to reactions while determining reaction dependency.
%       In step 1, PR reactions will have a cost of sqrt(om). In steps 1
%       and 2 NP reactions, as well as OT reactions in step 3, will be
%       assigned a cost of om. Default 1e+04.
%   ntimes - While determining reaction dependecy, simulations are
%       performed ntimes with the associated costs plus a small level of
%       noise, in order to determine pathways with the same associated
%       cost. Default value of 5.
%   nl - Noise level. Numeric. The noise level added to reaction costs will
%       be sampled uniformly between zero and nl.
% OUTPUTS:
%   tissue_model - tissue specific model produced (will be a subset of the
%       input model).
%   rescue - nx2 cell array. First colum is composed of strings of the
%       names of PR reactions not included in the model. The second column
%       is composed of cell arrays of the NP reactions implicated as being
%       necesary for the flux of the deleted PR reaction. This output will
%       help in manually curating the model if desired.
%   HCtoMC - dataset. Describes reaction dependencies between HC and MC
%       reactions. That is, if HCtoMC (i,j) = 1, MC reaction j was 
%       implicated as associated with HC reaction i. If HCtoMC(i,j) = 0
%       then the reactions are not associated.
%   HCtoNP - dataset. Similat to HCtoMC but defined the association between
%       HC and NP reactions defined in step one.
%   MCtoNP - dataset. Similat to HCtoMC but defined the association between
%       MC and NP reactions defined in step two.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tissue_model, rescue, HCtoMC, HCtoNC, MCtoNC] = CORDA(model,metTests,...
  ES,PR,NP,PRtoNP,constraint,constrainby,om,ntimes,nl)
h = waitbar(0,'initializing waitbar');
OT = model.rxns(~ismember(model.rxns,cat(1,ES,PR,NP)));

if (nargin < 6) || isempty(PRtoNP)
    PRtoNP = 2;
end
if (nargin < 7) || isempty(constraint)
    constraint = 1;
end
if (nargin < 8) || isempty(constrainby)
    constrainby = 'val';
end
if (nargin < 9) || isempty(om)
    om = 1e+04;
end
if (nargin < 10) || isempty(ntimes)
    ntimes = 5;
end
if (nargin < 11) || isempty(nl)
    nl = 1e-02;
end

%% Step 1
costbase = zeros(length(model.rxns),1);
PRids = findRxnIDs(model,PR);
NPids = findRxnIDs(model,NP);
PRpres = false(length(PR),1);
NPpres = false(length(NP),1);
costbase(PRids) = sqrt(om);
costbase(NPids) = om;
EStodelf = false(length(ES),1);
EStodelb = false(length(ES),1);
HCtoMC = zeros(length(ES),length(PR));
HCtoNC = zeros(length(ES),length(NP));

for j = 1:length(ES)
    waitbar(j/length(ES),h,'Step 1 - finding support reactions')
 
    %If reaction is sink then add it
    rem = false;
    if findRxnIDs(model,ES{j}) == 0
        id = find(strcmp(metTests(:,1),ES{j}));
        model = addReaction(model,metTests{id,1},metTests{id,2});
        costbase = [costbase; 0];
        rem = true;
    end
    
    %Test see which reactions are needed forward
    model = changeObjective(model,ES{j},1);
    flux = corsoFBA2(model,'max',constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
    if abs(flux.f) > 1e-6
        PRpres(abs(flux.x(PRids)) > 1e-6) = true;
        NPpres(abs(flux.x(NPids)) > 1e-6) = true;
        HCtoMC(j,abs(flux.x(PRids)) > 1e-6) = 1;
        HCtoNC(j,abs(flux.x(NPids)) > 1e-6) = 1;
        for k = 1:(ntimes-1)
            flux = corsoFBA2(model,'max',constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
            PRpres(abs(flux.x(PRids)) > 1e-6) = true;
            NPpres(abs(flux.x(NPids)) > 1e-6) = true;
            HCtoMC(j,abs(flux.x(PRids)) > 1e-6) = 1;
            HCtoNC(j,abs(flux.x(NPids)) > 1e-6) = 1;
        end
    else
        EStodelf(j) = true;
    end
    
    %Test see which reactions are needed backwards
    if model.lb(findRxnIDs(model,ES{j})) < 0
        flux = corsoFBA2(model,'min',-constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
        if abs(flux.f) > 1e-6
            PRpres(abs(flux.x(PRids)) > 1e-6) = true;
            NPpres(abs(flux.x(NPids)) > 1e-6) = true;
            HCtoMC(j,abs(flux.x(PRids)) > 1e-6) = 1;
            HCtoNC(j,abs(flux.x(NPids)) > 1e-6) = 1;
            for k = 1:(ntimes-1)
                flux = corsoFBA2(model,'min',-constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
                PRpres(abs(flux.x(PRids)) > 1e-6) = true;
                NPpres(abs(flux.x(NPids)) > 1e-6) = true;
                HCtoMC(j,abs(flux.x(PRids)) > 1e-6) = 1;
                HCtoNC(j,abs(flux.x(NPids)) > 1e-6) = 1;
            end
        else
            EStodelb(j) = true;
        end
    else
        EStodelb(j) = true;
    end
    
    if EStodelf(j) && EStodelb(j)
        fprintf(['Reaction ' ES{j} ' was deleted\n'])
    end
    
    %Remove sink if sink was added
    if rem
        model = removeRxns(model,ES{j});
        costbase(end) = [];
    end
end

%Tailor model
EStodel = EStodelf & EStodelb;
fprintf([num2str(length(find(EStodel))) ' blocked reactions removed from ES\n'])
ES(EStodel) = [];
HCtoMC(EStodel,:) = [];
HCtoNC(EStodel,:) = [];
HCtoMC = mat2dataset(HCtoMC,'varnames',PR,'obsnames',ES);
HCtoNC = mat2dataset(HCtoNC,'varnames',NP,'obsnames',ES);

fprintf([num2str(length(find(PRpres))) ' PR reactions added to ES\n'])
ES = cat(1,ES,PR(PRpres));
PR(PRpres) = [];
fprintf([num2str(length(find(NPpres))) ' NP reactions added to ES\n'])
ES = cat(1,ES,NP(NPpres));
NP(NPpres) = [];

%% Step 2
%% Step 2.1
%Assign costs
costbase = zeros(length(model.rxns),1);
costbase(findRxnIDs(model,NP)) = om;

%Initialize variables
NPid = findRxnIDs(model,NP);
PRxNP = zeros(length(PR),length(NP));
PRtodelf = false(length(PR),1);
PRtodelb = false(length(PR),1);
for j = 1:length(PR)
    waitbar(j/length(PR),h,'Step 2.1 - Checking PR and NP co-occurence')
    model = changeObjective(model,PR{j},1);
    
    %Check forward flux
    flux = corsoFBA2(model,'max',constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
    if isempty(flux.f)
        PRtodelf(j) = true;
    else
        for k = 1:(ntimes-1)
            PRxNP(j,abs(flux.x(NPid)) > 1e-6) = 1;
            flux = corsoFBA2(model,'max',constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
        end
        PRxNP(j,abs(flux.x(NPid)) > 1e-6) = 1;
    end
    
    %Check backwards flux
    if model.lb(findRxnIDs(model,PR{j})) < 0
        flux = corsoFBA2(model,'min',-constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
        if isempty(flux.f)
            PRtodelb(j) = true;
        else
            for k = 1:(ntimes-1)
                PRxNP(j,abs(flux.x(NPid)) > 1e-6) = 1;
                flux = corsoFBA2(model,'min',-constraint,constrainby,costbase+(nl*floor(10000*rand(length(costbase),1))/10000));
            end
            PRxNP(j,abs(flux.x(NPid)) > 1e-6) = 1;
        end
    else
        PRtodelb(j) = true;
    end

    if PRtodelf(j) && PRtodelb(j)
        fprintf(['Reaction ' PR{j} ' was deleted.\n'])
    end
end

PRxNP(PRtodelf & PRtodelb,:) = [];
PR(PRtodelf & PRtodelb) = [];
MCtoNC = mat2dataset(PRxNP,'varnames',NP,'obsnames',PR);

%% Step 2.2
%Add high occuring NPs to PR
t = sum(PRxNP);
ind = NP(t >= PRtoNP);
fprintf([num2str(length(ind)) ' reactions from NP are added to ES\n'])
%Fix occurence matrix
PR = cat(1,PR,ind);
PRxNP = [PRxNP;zeros(length(ind),length(NP))];
PRxNP(:,ismember(NP,ind)) = [];
NP(ismember(NP,ind)) = [];

%See which reactions from PR are no longer feasible
model = changeRxnBounds(model,NP,0,'b');
PRtodelf = false(length(PR),1);
PRtodelb = false(length(PR),1);
res1 = {};
res2 = {};

for j = 1:length(PR)
    waitbar(j/length(PR),h,'Step 2.2 - Checking PR feasibility')
    model = changeObjective(model,PR{j},1);
    
    %Check forward flux
    flux = optimizeCbModel(model,'max');
    if abs(flux.f) < 1e-6
        PRtodelf(j) = true;
    end
    
    %Check backwards flux
    if model.lb(findRxnIDs(model,PR{j})) < 0
        flux = optimizeCbModel(model,'min');
        if abs(flux.f) < 1e-6
            PRtodelb(j) = true;
        end
    else
        PRtodelb(j) = true;
    end
    
    if PRtodelf(j) && PRtodelb(j)
        fprintf([PR{j} ' was deleted. Dependent on: '])
        res1 = cat(1,res1,PR{j});
        tmp = find(PRxNP(j,:));
        if isempty(tmp)
            % NP reaction added can be dependent on other NP reactions that
            % appear less than PRtoNP times, and have thus been blocked.
            fprintf('Undefined')
            res2 = cat(1,res2,' ');
        else
            for k = 1:length(tmp)
                fprintf(NP{tmp(k)})
                if k ~= length(tmp)
                    fprintf(', ')
                end
                if k == 1
                    res2 = cat(1,res2,NP{tmp(k)});
                else
                    res2{length(res2)} = strcat(res2{length(res2)},',',NP{tmp(k)});
                end
            end
        end
        fprintf('\n')
    end
end
fprintf([num2str(length(find(PRtodelf & PRtodelb))) ' reactions deleted from PR\n'])
PR(PRtodelf & PRtodelb) = [];
ES = cat(1,ES,PR);
rescue = cat(2,res1,res2);

%% Step 3
%Block reactions not in ES or OT
model = changeRxnBounds(model,...
    model.rxns(~ismember(model.rxns,cat(1,ES,OT))),0,'b');
%Define cost
costbase = zeros(length(model.rxns),1);
costbase(findRxnIDs(model,OT)) = om;
%Parse through ES
OTid = findRxnIDs(model,OT);
ESxOT = zeros(length(ES),length(OT));
for j = 1:length(ES)
    waitbar(j/length(ES),h,'Step 3 - Define remaining reactions')
    %Add sink if sink is needed
    rem = false;
    if findRxnIDs(model,ES{j}) == 0
        id = find(strcmp(metTests(:,1),ES{j}));
        model = addReaction(model,metTests{id,1},metTests{id,2});
        costbase = [costbase; 0];
        rem = true;
    end
    model = changeObjective(model,ES{j},1);
    
    %optimize ES forward
    for k = 1:ntimes
        flux = corsoFBA2(model,'max',constraint,constrainby,costbase+((nl*floor(10000*rand(length(costbase),1))/10000)));
        ESxOT(j,abs(flux.x(OTid)) > 1e-6) = 1;
    end
    
    %optimize ES backwards
    if model.lb(findRxnIDs(model,ES{j})) < 0
        for k = 1:ntimes
            flux = corsoFBA2(model,'min',-constraint,constrainby,costbase+((nl*floor(10000*rand(length(costbase),1))/10000)));
            ESxOT(j,abs(flux.x(OTid)) > 1e-6) = 1;
        end
    end
    
    %Remove sink if sink was added
    if rem
        model = removeRxns(model,ES{j});
        costbase(end) = [];
    end
end

fprintf([num2str(length(find(sum(ESxOT)))) ' reactions added to ES for final model\n'])
ES = cat(1,ES,OT(sum(ESxOT) ~= 0));
tissue_model = removeRxns(model,model.rxns(~ismember(model.rxns,ES)));
close(h)
end

function flux = corsoFBA2(model,onstr,constraint,constrainby,costas)

if strcmp(constrainby,'perc')
    constraint = abs(constraint);
end
%determine objective function
flux1 = optimizeCbModel(model,onstr);
if abs(flux1.f) < 1e-6
%     warning('FBA problem infeasible')
    flux.f = [];
    flux.x = zeros(length(model.rxns),1);
    return
end

%relax results to avoid computational error
if strcmp(constrainby,'perc')
    flux1.f = flux1.x(model.c ~= 0)*(constraint/100);   %Bound
elseif strcmp(constrainby,'val')
    if (flux1.f < constraint) && strcmp(onstr,'max')
        error('Objective Flux not attainable')
    elseif (flux1.f > constraint) && strcmp(onstr,'min')
        error('Objective Flux not attainable')
    else
        flux1.f = constraint;   %Bound
    end
else
    error('Invalid Constraint option');
end
%save original model
model1 = model;

%See if cost is of right length
if length(costas) == 1
    costas = ones(length(model.rxns),1);
end
if ~iscolumn(costas)
    costas = costas';
end
if length(costas)==length(model.rxns)
    costas = [costas; costas];
elseif length(costas) ~= 2*length(model.rxns)
    fprintf('Invalid length of costs\n');
    flux = [];
    return
end

%find internal reactions that are actively reversible
orlen = length(model.rxns);
leng = find(model.lb<0 & model.ub>=0); 

%Tailor model
model.S = [model.S -model.S(:,leng); sparse(zeros(1,orlen+length(leng)))];
model.mets{length(model.mets)+1} = 'pseudomet';
model.b = zeros(length(model.mets),1);
model.c = zeros(orlen+length(leng),1);
model.ub = [model.ub; -model.lb(leng)];
model.lb = zeros(orlen+length(leng),1);
model.S(end,:) = [costas(1:orlen); costas(orlen+leng)];
model.rxns = cat(1,model.rxns,strcat(model.rxns(leng),'added'));

%add reaction for pseudomet consumption
model.rxns = cat(1,model.rxns,'EX_pseudomet');
model.ub = [model.ub; 1e20];
model.lb = [model.lb; 0];
temp = zeros(length(model.mets),1);
temp(end) = -1;
model.S = [model.S temp];
model.c = [model.c; 1];

%change bounds on original optimized reaction
t = find(model1.c);
for k = 1:length(t)
    model = changeRxnBounds(model,model1.rxns(t(k)),flux1.f(k),'b');
    if findRxnIDs(model,[model1.rxns{t(k)} 'added']) ~= 0
        model = changeRxnBounds(model,[model1.rxns{t(k)} 'added'],...
            0,'b');
    end
end

%perform FBA
flux2 = optimizeCbModel(model,'min');

flux.x = flux2.x(1:orlen);
flux.x(leng) = flux.x(leng) - flux2.x((orlen+1):(end-1));
flux.x(abs(flux.x) < 1e-8) = 0;

if isfield(flux1,'y')
    flux.y = flux1.y;
end
if isfield(flux1,'f')
    flux.f = flux1.f;
end
if isfield(flux2,'f')
    flux.fm = flux2.f;
end
end