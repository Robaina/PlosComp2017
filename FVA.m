function [Vmin,Vmax,Val]=FVA(GEM,L,solver,split,Blocked,eps)

%Calculates the Vmin, Vmax and flux range per each reaction of the model
%Inputs:
%GEM: a model in cobra structure format (at least S and rev fields)
%L: the maximum absolute value used as boudary condition, default is 1
%solver: string displaying either 'glpk' or 'gurobi'
%split: if TRUE, FVA returns the ranges of the extended matrix S', 
%constructed by splitting reversible reactions, i.e. S'=[S,-S(:,RevRxns)]
%Blocked: an index vector indicating reactions that should be blocked, i.e. zero
%flux
%Tasks: an index vector indicating reactions that should be always
%operative, i.e. non-zero flux

%Semidán, December 2014

%**************************************************************************

    if nargin<2,
        L=1;
    end
    if isempty(L),
        L=1;
    end
    if nargin<3,
        solver='gurobi';
    end
    if isempty(solver),
        solver='gurobi';
    end

    if nargin<4,
        split='F';
    end
    if isempty(split),
        split='F';
    end
    
    if nargin<5,
        Blocked=[];
    end
    if isempty(Blocked),
        Blocked=[];
    end
    
    if nargin<6,
        eps=1e-6*L;
    end
    if isempty(eps),
        eps=1e-6*L;
    end
    RevRxns=find(GEM.rev==1);
   
    if split=='F',
        S=full(GEM.S);
        Rxns=size(S,2);
        Mets=size(S,1);
        lb=zeros(Rxns,1);lb(RevRxns)=-L;
        ub=L*ones(Rxns,1);
%         if ~isempty(Blocked),
%             lb(Blocked)=0;
%             ub(Blocked)=0;
%         end
    elseif split=='T',
        S=[full(GEM.S),-full(GEM.S(:,RevRxns))];
        Rxns=size(S,2);
        Mets=size(S,1);
        
    end

    if strcmpi('gurobi',solver),
        n=1;m=1;
    %main loop: gurobi
        for i=1:Rxns,
            c=zeros(Rxns,1);c(i)=1;
            if split=='T', 
                lb=zeros(Rxns,1);
                ub=L*ones(Rxns,1);
                if ismember(i,RevRxns)
                    ub(size(GEM.S,2)+n)=0;
                    n=n+1;
                elseif i>size(GEM.S,2),
                    ub(RevRxns(m))=0;
                    m=m+1;
                end
            end
                
            model.A=sparse(S);
            model.modelsense='min';
            model.obj=c;
            model.sense=repmat('=',Mets,1);
            model.rhs=zeros(Mets,1);
            model.lb=lb;
            model.ub=ub;
            params.OutputFlag=0;

            gur=gurobi(model,params);
            Vmin(i)=gur.objval;

            model.A=sparse(S);
            model.modelsense='max';
            model.obj=c;
            model.sense=repmat('=',Mets,1);
            model.rhs=zeros(Mets,1);
            model.lb=lb;
            model.ub=ub;
            params.OutputFlag=0;

            gur=gurobi(model,params);
            Vmax(i)=gur.objval;
        end
    end

    if strcmpi('glpk',solver),
        n=1;
    %main loop: glpk
        for i=1:Rxns,
            c=zeros(Rxns,1);c(i)=1;
            if split=='T', 
                 lb=zeros(Rxns,1);
                 ub=L*ones(Rxns,1);
                if ismember(i,RevRxns),
                    ub(size(GEM.S,2)+n)=0;
                    n=n+1;
                elseif i>size(GEM.S,2),
                    ub(RevRxns(i-size(GEM.S,2)))=0;
                end
            end

            a=sparse(S);
            sense=1;
            ctype=repmat('S',Mets,1);
            vartype=repmat('C',Rxns,1);
            b=zeros(Mets,1);
            [~,fopt,~,~]=glpk(c,a,b,lb,ub,ctype,vartype,sense);
            Vmin(i)=fopt;

            a=sparse(S);
            sense=-1;
            ctype=repmat('S',Mets,1);
            vartype=repmat('C',Rxns,1);
            b=zeros(Mets,1);
            [~,fopt,~,~]=glpk(c,a,b,lb,ub,ctype,vartype,sense);
            Vmax(i)=fopt;
        end
    end

       Val.blocked=find(abs(Vmin)<=eps & abs(Vmax)<=eps);
       Val.VminMax(:,1)=Vmin';
       Val.VminMax(:,2)=Vmax';
       Val.RevRxns=RevRxns;
       Val.Rxns=Rxns;
       if split=='T', 
           VminMaxSplit=[Vmin(1:length(GEM.rxns))',Vmax(1:length(GEM.rxns))'];
           VminMaxSplit(GEM.rev==1,3:4)=[Vmin(length(GEM.rxns)+1:end)',Vmax(length(GEM.rxns)+1:end)'];
           VminMaxSplit=round(1e4*VminMaxSplit)/1e4;
           Val.VminMaxSplit=VminMaxSplit;
       end

end