function Val=coneSampling(GEM,split,Vrange,Vop,n,k,solver)
  %Performs a sampling of the feasible space determined by NV=0  
  %Argument evaluation
    h=0;
    if nargin<2 || isempty(split),
        split='F';
    end
    if nargin<3 || isempty(Vrange),
        [~,~,Val]=FVA(GEM,[],solver,split);
        Vrange=Val.VminMax;
    end
  
    if nargin<4 || isempty(Vop),
        Vop=Vrange(:,1)+((Vrange(:,2)-Vrange(:,1))/2);
        h=1;
    end

    if nargin<5,
        n=1e4;
    end
    if isempty(n),
        n=1e4;
    end
    
    if nargin<6,
        k=1;
    end
    if isempty(k),
        k=1;
    end
    
    if split=='T',
        GEM.S=[GEM.S,-GEM.S(:,find(GEM.rev==1))];
        if h==0,
            Vop=[Vop;Vop(find(GEM.rev==1))];
            if min(Vop(1:length(GEM.rxns)))<0,
                Vop(find(Vop(1:length(GEM.rxns))<0))=0;
            end
            if max(Vop((length(GEM.rxns)+1):length(Vop)))>0,
                Vop(find(Vop((length(GEM.rxns)+1):length(Vop)))>0)=0;
            end
        end
    end
    Rxns=size(GEM.S,2);
    Mets=size(GEM.S,1);
    S=full(GEM.S);
    Vss=zeros(Rxns,n);
    Idx=zeros(Rxns,2);
    Idx(:,1)=1:Rxns;
    Idx(GEM.rev==1,2)=Rxns+find(GEM.rev==1);

    for i=1:n,
    
        %Generate random perturbation in a region scaled by k
        
            em=k*(Vop-Vrange(:,1)).*rand(Rxns,1);
            ep=k*(Vrange(:,2)-Vop).*rand(Rxns,1);
            y=round(rand(Rxns,1));
            Vp=Vop+y.*ep-(1-y).*em;
        
        %Solve QP
        
        if strcmpi('gurobi',solver),
            model.Q=sparse(eye(Rxns));
            model.A=sparse(S);
            model.modelsense='min';
            model.obj=zeros(Rxns,1);
            model.sense=repmat('=',Mets,1);
            model.rhs=-S*Vp;
            if split=='F',
                model.lb=Vrange(:,1)-Vp;
                model.ub=Vrange(:,2)-Vp;
            elseif split=='T',
                model.lb=Vrange(:,1)-Vp;
                model.ub=Vrange(:,2)-Vp;
            end
            params.OutputFlag=0;
            gur=gurobi(model,params);
            X=gur.x;
        elseif strcmpi('cplex',solver),
            H=sparse(eye(Rxns));
            f=zeros(Rxns,1);
            Aeq=sparse(S);
            Aineq=[];
            beq=-S*Vp;
            bineq=[];
            if split=='F',
                lb=Vrange(:,1)-Vp;
                ub=Vrange(:,2)-Vp;
            elseif split=='T',
                lb=Vrange(:,1)-Vp;
                ub=Vrange(:,2)-Vp;
            end
            cplexoptions=cplexoptimset;
            cplexoptions.Display='off';
            X=cplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub,[],cplexoptions);
        end
            
        %Steady-state V in the vicinity
        Vss(:,i)=Vp+X;
        V(:,i)=Vp;
    
    end
Val.Vsample=Vss;
Val.Vrand=V;
Val.Vrange=Vrange;
Val.Idx=Idx;
end