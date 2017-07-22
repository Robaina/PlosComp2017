function contextCOBRA = contextmodel2COBRA(V,GEM,contextname,d)

% Function to transform a flux distribution into a COBRA model structure
% Arguments
% V: vector of flux values
% GEM: genome-scale metabolic model COBRA structure
% contextname: character array displaying the name that the context-specific 
% model should have (default set to empty name)
% d: binary, if 1, V must be the index set of active reactions
%
% Semidan, October, 2014
%**************************************************************************

if nargin<4,
    d=0;
end
if isempty(d),
    d=0;
end
eps=0;
if d==0,
   contextrxns=find(abs(V)>eps);
elseif d==1,
    contextrxns=V;
end

fieldname=fields(GEM);
% exGEMgenes=zeros(length(GEM.genes),1);
exGEMmets=zeros(length(GEM.mets),1);
if nargin<3,
    contextname=[];
end
if isempty(contextname),
    contextname=[];
end

% for j=1:length(exGEMgenes),
%     exGEMgenes(j)=sum(GEM.rxnGeneMat(contextrxns,j));
% end
% contextgenes=find(exGEMgenes~=0);
contextgenes=1:length(GEM.genes);

for j=1:length(exGEMmets),
    exGEMmets(j)=sum(abs(GEM.S(j,contextrxns)));
end
contextmets=find(exGEMmets~=0);
    

for i=1:length(fieldname),
    fieldnametype=size(GEM.(fieldname{i}));
    if fieldnametype(2)==1 && length(GEM.(fieldname{i}))==length(GEM.rxns),
        contextCOBRA.(fieldname{i})=GEM.(fieldname{i})(contextrxns);
    elseif fieldnametype(2)==1 && length(GEM.(fieldname{i}))==length(GEM.mets),
        contextCOBRA.(fieldname{i})=GEM.(fieldname{i})(contextmets);
    elseif fieldnametype(2)==1 && length(GEM.(fieldname{i}))==length(GEM.genes),
        contextCOBRA.(fieldname{i})=GEM.(fieldname{i})(contextgenes);
    elseif fieldnametype(2)>1 && size(GEM.(fieldname{i}),2)==length(GEM.rxns),
        contextCOBRA.(fieldname{i})=GEM.(fieldname{i})(contextmets,contextrxns);
    elseif fieldnametype(2)>1 && size(GEM.(fieldname{i}),2)==length(GEM.genes),
        contextCOBRA.(fieldname{i})=GEM.(fieldname{i})(contextrxns,contextgenes);
    elseif size(GEM.(fieldname{i}),1)>length(contextrxns) || size(GEM.(fieldname{i}),1)>length(contextmets),
        contextCOBRA.(fieldname{i})=[];
    else
        contextCOBRA.(fieldname{i})=contextname;
    end
end
end
    

