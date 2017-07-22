function Sol=RegrExLAD(GEM,D,Options)
      
%************************RegrExLAD Implementation*****************************

%Depends on:
%
% Gurobi, which can be found in www.gurobi.com (free academic licenses 
% available)
% contextmodel2COBRA function (provided in supplementary material)
% FVA function (provided in supplementary material)
%
%Arguments:
%
%     Required
%
% GEM: Metabolic model in COBRA structure. Blocked reactions should be 
% removed first, for instance using the reduceModel function of the COBRA 
% toolbox
% D: Experimental data vector (values already mapped to reactions) must 
% be of same length as the number of reactions. It can be provided by the 
% genetorxn function
% 
%     Options:
%
% Options must be a structure with the following fields (not all required)
%
%     Options.Dstd: a vector containing the standard deviation of data, it is 
%     intended to be used as weighting factor for reactions with associated  
%     data, in the form (1/CV(D))*||V-D||_1, where CV=std(D)/mean(D) 
%     Options.essmet: A vector indicating the indexes of the metabolites that should be 
%     present in the final model
%     Options.tasks: A vector indicating the indexes of the reactions that should be
%     active in the final model (i.e. flux value above epsilon)
%     Options.blocked: A vector indicating the indexes of the reactions that should not
%     be active in the final model (i.e. zero flux value)
%     Options.L_min: Minimum lambda value in the lambda sequence (default is set to 0)
%     Options.L_max: Maximum lambda value in the lambda sequence (default is set to .2)
%     Options.L_int: Step size used to generate the lambda sequence (default is set to
%     0.1)
%     Options.rxns_lb: a vector containing the lower bounds for each reaction (default
%     is set to 0, min value is 0, reversible reactions are splitted)
%     Options.rxns_ub: a vector containing the upper bounds for each reaction (default 
%     is set to 1)
%     Options.time_limit: scalar indicating the time limit in seconds for the MIQP
%     (default is set to 60)
%     Options.OutFlag: selects if the gurobi solution's flag should be
%     included in the console, must be 0 or 1 (default 0)
%     Options.epsilon: scalar indicating the threshold value to consider a reaction active,
%     i.e. a reaction i with flux value vi is active if vi>=epsilon
%     Options.E_max: vector of length equal to the number of reactions indicating the 
%     maximum allowed error value between data and flux prediction, i.e. E=d-v,
%     (default is set to 1000)
%     Options.contextname: character string indicating the name of the context to
%     display in the COBRA structure (default value is empty)
% 
%     Value:
%
%    Sol: A structure with the fields,
%    Sol.Flux: flux distribution
%    Sol.Cor: correlation value between fluxes and data
%    Sol.ZE: total error value (residual)
%    Sol.ActRxn: set of active (i.e. |Vi|>epsilon) reactions
%    Sol.ActRxnData: set of active reactions with associated data
%    Sol.Card: Number of active reactions
%    Sol.CardData: Number of active reactions with associated data
%    Sol.Lambda: used lambda-sequence
%
%    ContextCOBRA: A structure in COBRA format (see contextModel2COBRA
%    function for details
%
%**************************************************************************
%         Semidán (robaina@mpimp-golm.mpg.de), May, 2016
%**************************************************************************


%Argument evaluation

if ~exist('GEM','var'),
    error('The Genome-Scale model is missing...')
end

if ~exist('D','var'),
    error('The data vector is missing...')
end

if ~exist('Options','var'),
    Options=struct;
end

if ~isfield(Options,'Dstd'),
    Dstd=1;
else
    Dstd=Options.Dstd;
end

if ~isfield(Options,'essmet'), 
    essmet=[];
else
    essmet=Options.essmet;
end

if ~isfield(Options,'tasks'), 
    tasks=[];
else
    tasks=Options.tasks;
end

if ~isfield(Options,'blocked'), 
    blocked=[];
else
    blocked=Options.blocked;
end

if ~isfield(Options,'L_min'), 
    L_min=0;
else
    L_min=Options.L_min;
end

if ~isfield(Options,'L_max'), 
    L_max=0.1;
else
    L_max=Options.L_max;
end

if ~isfield(Options,'L_int'), 
    L_int=0.01;
else
    L_int=Options.L_int;
end

if ~isfield(Options,'MaxCapacity'), 
    MaxCapacity=1;
else
    MaxCapacity=Options.MaxCapacity;
end

if ~isfield(Options,'rxns_lb'), 
    rxns_lb=0;
else
    rxns_lb=Options.rxns_lb;
end

if ~isfield(Options,'rxns_ub'), 
    rxns_ub=MaxCapacity;
else
    rxns_ub=Options.rxns_ub;
end

if ~isfield(Options,'MIPGap'), 
    MIPGap=0.005;
else
    MIPGap=Options.MIPGap;
end

if ~isfield(Options,'OutFlag'), 
    OutFlag=0;
else
    OutFlag=Options.OutFlag;
end

if ~isfield(Options,'epsilon'), 
    epsilon=MaxCapacity*1e-6;
else
    epsilon=Options.epsilon;
end

if ~isfield(Options,'E_max'), 
    E_max=MaxCapacity*1e3;
else
    E_max=Options.E_max;
end

if ~isfield(Options,'contextname'), 
    contextname=[];
else
    contextname=Options.contextname;
end

    %Parse Data
    
    S=GEM.S;
    Rev=find(GEM.rev==1);
    Irr=(setdiff(1:size(S,2),Rev))';
    
    for i=1:length(D),
        if isnan(D(i)) || isinf(D(i)),
            D(i)=0;
            sprintf('Error: Data point D(%d) is NaN or Inf and has been converted to 0',i);
        end
        if length(Dstd)>1,
            if isnan(Dstd(i)) || isinf(Dstd(i)),
               Dstd(i)=1;
               sprintf('Error: Data std point Dstd(%d) is NaN or Inf and has been converted to 1',i);
            end
        end
    end
    if max(D)>MaxCapacity,
       D=D/max(D);
       Dstd=Dstd/max(D);
    end
    D=MaxCapacity*D;
    Dor=D;
    
    %Reaction partition: Irr_D, Irr_nD, Bfordat, Bfornodat, Rev_D, Rev_nD
    
    Irr_D=Irr((D(Irr)~=0));
    Irr_nD=Irr((D(Irr)==0));
    Rev_D=Rev((D(Rev)~=0));
    Rev_nD=Rev((D(Rev)==0));
    NIrr_D=length(Irr_D);
    NRev_D=length(Rev_D);
    NRev_nD=length(Rev_nD);
    Rxns_Or=[Irr_D;Irr_nD;Rev_D;Rev_nD;Rev_D;Rev_nD];
    
    %Stoichiometric Matrix reorganization
    
    Sam=[S(:,Irr_D),S(:,Irr_nD),S(:,Rev_D),S(:,Rev_nD),-S(:,Rev_D),-S(:,Rev_nD)];  
    Rxns=size(Sam,2);Mets=size(Sam,1);
    NIrr=length(Irr);NRev=length(Rev);
    DIrr=D(Irr_D);
    DRev=D(Rev_D);
    DFor=DRev;
    
    if length(Dstd)>1,
        DstdIrr=Dstd(Irr_D);
        DstdFor=Dstd(Rev_D);
        DstdRev=DstdFor;
    end
        
    %Construction of L1 sequence
    
    if L_min==L_max,
        L=L_min;
    else
    L=L_min:L_int:L_max;
    end
 
    V=zeros(size(S,2),length(L));
    Dcor=zeros(length(L),1);
    Ddis=Dcor;
    Card=Dcor;
    CardData=Dcor;
    IterTime=Dcor;
    
    %Construction of bounding vectors
    
    ubrxns=MaxCapacity*ones(Rxns,1);
    lbrxns=zeros(Rxns,1);
    
    if length(rxns_lb)~=1,
      for h=1:size(S,2),
          lbrxns(ismember(Rxns_Or,h))=rxns_lb(h);
      end
    end
    
    if length(rxns_ub)~=1,
        for h=1:size(S,2),
            ubrxns(ismember(Rxns_Or,h))=rxns_ub(h);
        end
    end
    
    if ~isempty(tasks),
       lbrxns((ismember(Rxns_Or,tasks)))=epsilon;
    end
    
    if ~isempty(blocked),
       ubrxns((ismember(Rxns_Or,blocked)))=epsilon;
    end
    
    Revlb=lbrxns(Rev);
    Revub=ubrxns(Rev);
    lb=[zeros(2*NIrr_D+4*NRev_D,1);lbrxns;zeros(NRev,1)];
    ub=[repmat(E_max,(2*NIrr_D+4*NRev_D),1);ubrxns;ones(NRev,1)];
    
    %Include constraints on essential metabolites
    
    if ~isempty(essmet),
        essmetmat=zeros(length(essmet),Rxns);
        for k=1:length(essmet),
            essmetmat(k,Sam(essmet(k),:)>0)=1;
        end
        vecsense=[repmat('=',Mets,1);repmat('=',(NIrr_D+NRev_D),1);repmat('<',NRev,1);repmat('=',NRev_D,1);repmat('<',NRev,1);repmat('>',NRev,1);repmat('>',NRev,1);repmat('>',length(essmet),1)];
        b=[zeros(Mets,1);DIrr;DFor;Revub;zeros(NRev_D,1);zeros(NRev,1);Revlb;zeros(NRev,1);epsilon*ones(length(essmet),1)];
    else
        vecsense=[repmat('=',Mets,1);repmat('=',(NIrr_D+NRev_D),1);repmat('<',NRev,1);repmat('=',NRev_D,1);repmat('<',NRev,1);repmat('>',NRev,1);repmat('>',NRev,1)];
        b=[zeros(Mets,1);DIrr;DFor;Revub;zeros(NRev_D,1);zeros(NRev,1);Revlb;zeros(NRev,1)];
    end
    
    %Construction of A matrix: E+irrdat, E-irrdat, E+for, E-for, E+rev, E-rev,Irr_D, Irr_nD, For_D,
    %For_nD, Rev_D, Rev_nD, Y
    
    A0=[zeros(Mets,2*NIrr_D+4*NRev_D),Sam,zeros(Mets,NRev)]; %SV=0
    A1=[diag(ones(NIrr_D,1)),-diag(ones(NIrr_D,1)),zeros(NIrr_D,4*NRev_D),diag(ones(NIrr_D,1)),zeros(NIrr_D,(NIrr-NIrr_D)),zeros(NIrr_D,3*NRev)]; %V+E=D, IrrRxns
    A2=[zeros(NRev_D,2*NIrr_D),diag(ones(NRev_D,1)),-diag(ones(NRev_D,1)),zeros(NRev_D,2*NRev_D),zeros(NRev_D,NIrr),diag(ones(NRev_D,1)),zeros(NRev_D,NRev_nD),zeros(NRev_D,NRev),diag(DFor),zeros(NRev_D,NRev_nD)]; %V+E+y*D=D, VforRxns
    A3=[zeros(NRev,2*NIrr_D+4*NRev_D+NIrr),diag(ones(NRev,1)),zeros(NRev,NRev),diag(Revub)];%Vfor+y*Vformax<=Vformax
    A4=[zeros(NRev_D,2*NIrr_D),zeros(NRev_D,2*NRev_D),diag(ones(NRev_D,1)),-diag(ones(NRev_D,1)),zeros(NRev_D,NIrr),zeros(NRev_D,NRev),diag(ones(NRev_D,1)),zeros(NRev_D,NRev_nD),diag(-DRev),zeros(NRev_D,NRev_nD)]; %V+E-y*D=0, VrevRxns
    A5=[zeros(NRev,2*NIrr_D+4*NRev_D+NIrr+NRev),diag(ones(NRev,1)),-diag(Revub)]; %Vrev-y*Vrevmax<=0
    A6=[zeros(NRev,2*NIrr_D+4*NRev_D+NIrr),diag(ones(NRev,1)),zeros(NRev,NRev),diag(Revlb)];%Vfor+y*Vformin>=Vformin
    A7=[zeros(NRev,2*NIrr_D+4*NRev_D+NIrr+NRev),diag(ones(NRev,1)),-diag(Revlb)]; %Vrev-y*Vrevmin>=0

    if ~isempty(essmet),
        A14=zeros(length(essmet),size(A0,2));
        for k=1:length(essmet),
            A14(k,:)=[zeros(1,2*NIrr_D+4*NRev_D),essmetmat(k,:),zeros(1,NRev)];
        end
    Amat=[A0;A1;A2;A3;A4;A5;A6;A7];
    elseif isempty(essmet), 
       Amat=[A0;A1;A2;A3;A4;A5;A6;A7];
    end
        
    %Integrates weighting based on standard deviation of each data point
    %(W=mean/std)
   
    if length(Dstd)>1,
%         W=[DIrr;DIrr;DFor;DFor;DRev;DRev]./[DstdIrr;DstdIrr;DstdFor;DstdFor;DstdRev;DstdRev];
        W=1./[DstdIrr;DstdIrr;DstdFor;DstdFor;DstdRev;DstdRev];
        W(W==Inf)=max(W);
    else
        W=Dstd;
    end

    %Iteration to determine optimum L
    
    for i=1:length(L),
        tic
        Lvec=repmat(L(i),Rxns,1);
        cvec=[W.*ones(2*NIrr_D+4*NRev_D,1);Lvec;zeros(NRev,1)];
        
        m.obj=cvec;
        m.A=sparse(Amat);
        m.rhs=b;
        m.sense=vecsense;
        m.vtype=[repmat('C',2*NIrr_D+4*NRev_D+Rxns,1);repmat('B',NRev,1)];
        m.lb=lb;
        m.ub=ub;
        params.OutputFlag=OutFlag;
        params.Presolve=2;
        params.MIPGap=MIPGap;
        gur=gurobi(m,params);
        
        X=gur.x;
        Virr_D=X((2*NIrr_D+4*NRev_D+1):(2*NIrr_D+4*NRev_D+NIrr_D));
        Virr_nD=X((2*NIrr_D+4*NRev_D+NIrr_D+1):(2*NIrr_D+4*NRev_D+NIrr));
        Vfor_D=X((2*NIrr_D+4*NRev_D+NIrr+1):(2*NIrr_D+4*NRev_D+NIrr+NRev_D));
        Vfor_nD=X((2*NIrr_D+4*NRev_D+NIrr+NRev_D+1):(2*NIrr_D+4*NRev_D+NIrr+NRev));
        Vrev_D=X((2*NIrr_D+4*NRev_D+NIrr+NRev+1):(2*NIrr_D+4*NRev_D+NIrr+NRev+NRev_D));
        Vrev_nD=X((2*NIrr_D+4*NRev_D+NIrr+NRev+NRev_D+1):(2*NIrr_D+4*NRev_D+NIrr+2*NRev));
        V(:,i)=zeros(size(S,2),1);V(Irr_D,i)=Virr_D;V(Irr_nD,i)=Virr_nD;
        V(Rev_D,i)=Vfor_D-Vrev_D;V(Rev_nD,i)=Vfor_nD-Vrev_nD;
        V((abs(V(:,i))<epsilon),i)=0;
        
        if max(abs(V(:,i)))>0,
            Cor=corrcoef(abs(V((Dor~=0),i)),Dor((Dor~=0)));
            Dcor(i)=Cor(1,2);
        end
        Ddis(i)=sum(sqrt((abs(V((Dor~=0),i))-Dor((Dor~=0))).^2))/(NIrr_D+NRev_D);
        Card(i)=length(find(abs(V(:,i))>0));
        CardData(i)=length(intersect(find(abs(V(:,i))>0),find(Dor>0)));
        IterTime(i)=toc; 
    end
    
    L_opt=find(Dcor==max(Dcor));
    if length(Dstd)>1,
        Weight=zeros(size(V,1),1);
        Weight(Irr_D)=W(1:NIrr_D);Weight(Rev_D)=W(2*NIrr_D+1:2*NIrr_D+NRev_D);
        Sol.W=Weight;
        ZE=sum(W.*gur.x(1:2*NIrr_D+4*NRev_D));
    elseif length(Dstd)==1,
        ZE=sum(gur.x(1:2*NIrr_D+4*NRev_D));
    end
    Sol.Flux=V(:,L_opt);
    Sol.Cor=Dcor(L_opt);
    Sol.Res=Ddis(L_opt);
    Sol.ActRxn=find(abs(V(:,L_opt))>0);
    Sol.ActRxnData=intersect(find(abs(V(:,L_opt))>0),find(Dor>0));
    Sol.Card=Card(L_opt);
    Sol.CardData=CardData(L_opt);
    Sol.Lambda=L(L_opt);
    Sol.ZE=ZE;
    Sol.ZV=sum(abs(Sol.Flux));
    Sol.Time=IterTime;
    if ~isempty(contextname),
       Sol.ContextCOBRA=contextmodel2COBRA(V(:,L_opt),GEM,contextname,0,0);
    end
end
