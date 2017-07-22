function [RxnMechanism]=createRxnMechanism(GEM,Flux,addcomp)

%Function to create reaction mechanisms for a GEM, also displays the
%cellular compartment where the reaction occurs
%
%Arguments:
%GEM: model COBRA structure
%Flux: predicted flux distribution (optional), if Flux is provided then the
%direction of the reversible reactions is fixed according to the predicted
%flux value
%
%Semidán November, 2014
%
%**************************************************************************
if nargin<2,
    Flux=0;
end
if isempty(Flux)==1,
    Flux=0;
end

if nargin<3,
    addcomp='T';
end
if isempty(addcomp)==1,
    addcomp='T';
end

%Add compartment to metabolite name

if addcomp=='T',
    for i=1:length(GEM.metNames),
    GEM.metNamesC{i}=strcat(GEM.metNames{i},GEM.mets{i}(end-2:end));
    end
elseif addcomp=='F',
    GEM.metNamesC=GEM.metNames;
end

for i=1:length(GEM.rxns),
    
    Substrates=GEM.metNamesC(GEM.S(:,i)<0);
    Products=GEM.metNamesC(GEM.S(:,i)>0);
    
    if isempty(Substrates)==0,
    S=Substrates(1);
    if length(Substrates)>1,
    for k=2:length(Substrates),
        S=strcat(S,{' '},'+',{' '},Substrates(k));
    end
    end
    elseif isempty(Substrates)==1,
        S='Extracellular';
    end
    
    if isempty(Products)==0,
    P=Products(1);
    if length(Products)>1,
    for j=2:length(Products),
        P=strcat(P,{' '},'+',{' '},Products(j));
    end
    end
    elseif isempty(Products)==1,
        P='Extracellular';
    end
    
    %Assigns reaction direction when required
    
    if length(Flux)==1,
        if GEM.rev(i)==0,
        RxnMechanism{i}=strcat(S,{' '},'==>',{' '},P);
        elseif GEM.rev(i)==1,
        RxnMechanism{i}=strcat(S,{' '},'<==>',{' '},P);
        end
    end
    
    if length(Flux)>1,
        if Flux(i)>0,
            RxnMechanism{i}=strcat(S,{' '},'==>',{' '},P);
        elseif Flux(i)<0,
            RxnMechanism{i}=strcat(S,{' '},'<==',{' '},P);
        elseif Flux(i)==0,
            if GEM.rev(i)==0,
               RxnMechanism{i}=strcat(S,{' '},'==>',{' '},P,{' '},', INACTIVE REACTION');
            elseif GEM.rev(i)==1,
               RxnMechanism{i}=strcat(S,{' '},'<==>',{' '},P,{' '},', INACTIVE REACTION');
            end
        end
    end     
            
end

RxnMechanism=RxnMechanism';

end

