

%map expression data to AraCOREred
MapGeneValues=zeros(length(AraCOREred.rxns),3);
GeneValues=GeneValues/max(max(GeneValues));
for i=1:3,
  MapGeneValues(:,i)=mapgene2rxn(GeneNames,GeneValues(:,i),AraCOREred.genes,AraCOREred.grRules,0);
end
% LeafData=zeros(length(AraCOREred.rxns),2);

for i=1:length(AraCOREred.rxns),
    LeafData(i,1)=mean(MapGeneValues(i,:));
    LeafData(i,2)=std(MapGeneValues(i,:));
end
LeafData=LeafData/max(LeafData(:,1));

%**************************************************************************
%****************************CorEx AO Analysis*****************************
%**************************************************************************

%LEAF (AraCOREred)

CLeaf=find(LeafData(:,1)>=quantile(LeafData(LeafData(:,1)>0,1),0.7));
%get CorEx model
LeafCorEx=CorEx(AraCOREred,CLeaf,[],1000,0.001)

%get FastCORE model
LeafFastCORE=fastcore(CLeaf,AraCOREred,1e-5);
PleafFast=length(setdiff(LeafFastCORE,CLeaf));

%get alternative optimal networks
LeafCorExAO=AltNet(AraCOREred,CLeaf,LeafCorEx.AddedRxns,LeafCorEx.Abinary,1000,0.001,100,5*60,100)
LeafFastCOREAO=AltNet(AraCOREred,CLeaf,PleafFast,LeafFastCORE,1000,0.001,100,5*60,100)

%get the fixed and variable fractions of the total number of non-core
%reactions among alternative networks
LeafCorExtotnoC=LeafCorExAO.Modmatrix;LeafCorExtotnoC(CLeaf,:)=-1;LeafCorExtotnoC=sum(LeafCorExtotnoC')';FleafCorExnoC=find(LeafCorExtotnoC==size(LeafCorExAO.Modmatrix,2));VleafCorExnoC=find(LeafCorExtotnoC>0 & LeafCorExtotnoC<size(LeafCorExAO.Modmatrix,2));
LeafFastCOREtotnoC=LeafFastCOREAO.Modmatrix;LeafFastCOREtotnoC(CLeaf,:)=-1;LeafFastCOREtotnoC=sum(LeafFastCOREtotnoC')';FleafCorExFastnoC=find(LeafFastCOREtotnoC==size(LeafFastCOREAO.Modmatrix,2));VleafCorExFastnoC=find(LeafFastCOREtotnoC>0 & LeafFastCOREtotnoC<size(LeafFastCOREAO.Modmatrix,2));

%rank reactions according to its frequency of ocurrence in the AO space of
%models
[LeafFreqCorEx,idx1]=sort(LeafCorExtotnoC/max(LeafCorExtotnoC),'descend');
[LeafFreqFastCORE,idx2]=sort(LeafFastCOREtotnoC/max(LeafFastCOREtotnoC),'descend');
idx1(LeafFreqCorEx<0)=[];LeafFreqCorEx(LeafFreqCorEx<0)=[];
idx2(LeafFreqFastCORE<0)=[];LeafFreqFastCORE(LeafFreqFastCORE<0)=[];
LeafCorExRxnRanking=[{'Frequency in AO','Rxn Name','Rxn SubSystem'};[num2cell(LeafFreqCorEx),AraCOREred.rxnNames(idx1),AraCOREred.subSystems(idx1)];{'Fixed Active','Fixed Inactive','Variable'};num2cell([length(find(LeafFreqCorEx==1)),length(find(LeafFreqCorEx==0)),length(find(LeafFreqCorEx~=1 & LeafFreqCorEx~=0))])];
LeafFastCORERxnRanking=[{'Frequency in AO','Rxn Name','Rxn SubSystem'};[num2cell(LeafFreqFastCORE),AraCOREred.rxnNames(idx2),AraCOREred.subSystems(idx2)];{'Fixed Active','Fixed Inactive','Variable'};num2cell([length(find(LeafFreqFastCORE==1)),length(find(LeafFreqFastCORE==0)),length(find(LeafFreqFastCORE~=1 & LeafFreqFastCORE~=0))])];

%get the distribution of reactions in subsystems in the fixed and the
%variable non-core set
B=[];C=[];D=[];E=[];
NCSubSys=unique(LeafCorExRxnRanking(2:(end-2),3));
for i=1:length(NCSubSys),
    E(i,1)=sum(strcmp(LeafCorExRxnRanking(2:(end-2),3),NCSubSys{i}));
end

ActiveNCSubSys=LeafCorExRxnRanking(find(cell2mat(LeafCorExRxnRanking(2:(end-2),1))==1)+1,3);
VarNCSubSys=LeafCorExRxnRanking(find(cell2mat(LeafCorExRxnRanking(2:(end-2),1))<1 & cell2mat(LeafCorExRxnRanking(2:(end-2),1))>0)+1,3);
InactiveNCSubSys=LeafCorExRxnRanking(find(cell2mat(LeafCorExRxnRanking(2:(end-2),1))==0)+1,3);

UniqActiveNCSubSys=unique(ActiveNCSubSys);UniqVarNCSubSys=unique(VarNCSubSys);UniqInactiveNCSubSys=unique(InactiveNCSubSys);
for i=1:length(UniqInactiveNCSubSys),
    B(i,1)=sum(strcmp(InactiveNCSubSys,UniqInactiveNCSubSys{i}))/E(strcmp(NCSubSys,UniqInactiveNCSubSys{i}));
end
for i=1:length(UniqActiveNCSubSys),
    C(i,1)=sum(strcmp(ActiveNCSubSys,UniqActiveNCSubSys{i}))/E(strcmp(NCSubSys,UniqActiveNCSubSys{i}));
end
for i=1:length(UniqVarNCSubSys),
    D(i,1)=sum(strcmp(VarNCSubSys,UniqVarNCSubSys{i}))/E(strcmp(NCSubSys,UniqVarNCSubSys{i}));
end
[B,idx2]=sort(B,'descend');[C,idx3]=sort(C,'descend');[D,idx4]=sort(D,'descend');
LeafCorExInactSubSys=[UniqInactiveNCSubSys(idx2),num2cell(B)];
LeafCorExActSubSys=[UniqActiveNCSubSys(idx3),num2cell(C)];
LeafCorExVarSubSys=[UniqVarNCSubSys(idx4),num2cell(D)];

B=[];C=[];D=[];E=[];
NCSubSys=unique(LeafFastCORERxnRanking(2:(end-2),3));
for i=1:length(NCSubSys),
    E(i,1)=sum(strcmp(LeafFastCORERxnRanking(2:(end-2),3),NCSubSys{i}));
end

ActiveNCSubSys=LeafFastCORERxnRanking(find(cell2mat(LeafFastCORERxnRanking(2:(end-2),1))==1)+1,3);
VarNCSubSys=LeafFastCORERxnRanking(find(cell2mat(LeafFastCORERxnRanking(2:(end-2),1))<1 & cell2mat(LeafFastCORERxnRanking(2:(end-2),1))>0)+1,3);
InactiveNCSubSys=LeafFastCORERxnRanking(find(cell2mat(LeafFastCORERxnRanking(2:(end-2),1))==0)+1,3);

UniqActiveNCSubSys=unique(ActiveNCSubSys);UniqVarNCSubSys=unique(VarNCSubSys);UniqInactiveNCSubSys=unique(InactiveNCSubSys);
for i=1:length(UniqInactiveNCSubSys),
    B(i,1)=sum(strcmp(InactiveNCSubSys,UniqInactiveNCSubSys{i}))/E(strcmp(NCSubSys,UniqInactiveNCSubSys{i}));
end
for i=1:length(UniqActiveNCSubSys),
    C(i,1)=sum(strcmp(ActiveNCSubSys,UniqActiveNCSubSys{i}))/E(strcmp(NCSubSys,UniqActiveNCSubSys{i}));
end
for i=1:length(UniqVarNCSubSys),
    D(i,1)=sum(strcmp(VarNCSubSys,UniqVarNCSubSys{i}))/E(strcmp(NCSubSys,UniqVarNCSubSys{i}));
end
[B,idx2]=sort(B,'descend');[C,idx3]=sort(C,'descend');[D,idx4]=sort(D,'descend');
LeafFastCOREInactSubSys=[UniqInactiveNCSubSys(idx2),num2cell(B)];
LeafFastCOREActSubSys=[UniqActiveNCSubSys(idx3),num2cell(C)];
LeafFastCOREVarSubSys=[UniqVarNCSubSys(idx4),num2cell(D)];

%calculate hamming distance
HmatLeafCorEx=HammingMat(LeafCorExAO.Modmatrix);
HmatLeafFastCORE=HammingMat(LeafFastCOREAO.Modmatrix);

%get distribution of hamming distances 
HDisLeafCorEx=HmatLeafCorEx(triu(logical(ones(size(HmatLeafCorEx,1))),1));
HDisLeafFastCORE=HmatLeafFastCORE(triu(logical(ones(size(HmatLeafFastCORE,1))),1));

%LIVER (Recon1red)

%Cliver <- liver CORE in FastCORE publication
%get CorEx model
LiverCorEx=CorEx(Recon1red,CLiver,[],1000,0.001)

%get FastCORE model
LiverFastCORE=fastcore(CLiver,Recon1red,1e-5);
PliverFast=length(setdiff(LiverFastCORE,CLiver));

%get alternative optimal networks
LiverCorExAO=AltNet(Recon1red,CLiver,LiverCorEx.AddedRxns,LiverCorEx.Abinary,1000,0.001,100,5*60,100)
LiverFastCOREAO=AltNet(Recon1red,CLiver,PliverFast,LiverFastCORE,1000,0.001,100,5*60,100)

%get the fixed and variable fractions of the total number of non-core
%reactions among alternative networks
LiverCorExtotnoC=LiverCorExAO.Modmatrix;LiverCorExtotnoC(CLiver,:)=-1;LiverCorExtotnoC=sum(LiverCorExtotnoC')';FliverCorExnoC=find(LiverCorExtotnoC==size(LiverCorExAO.Modmatrix,2));VliverCorExnoC=find(LiverCorExtotnoC>0 & LiverCorExtotnoC<size(LiverCorExAO.Modmatrix,2));
LiverFastCOREtotnoC=LiverFastCOREAO.Modmatrix;LiverFastCOREtotnoC(CLiver,:)=-1;LiverFastCOREtotnoC=sum(LiverFastCOREtotnoC')';FliverCorExFastnoC=find(LiverFastCOREtotnoC==size(LiverFastCOREAO.Modmatrix,2));VliverCorExFastnoC=find(LiverFastCOREtotnoC>0 & LiverFastCOREtotnoC<size(LiverFastCOREAO.Modmatrix,2));

%rank reactions according to its frequency of ocurrence in the AO space of
%models
[LiverFreqCorEx,idx1]=sort(LiverCorExtotnoC/max(LiverCorExtotnoC),'descend');
[LiverFreqFastCORE,idx2]=sort(LiverFastCOREtotnoC/max(LiverFastCOREtotnoC),'descend');
idx1(LiverFreqCorEx<0)=[];LiverFreqCorEx(LiverFreqCorEx<0)=[];
idx2(LiverFreqFastCORE<0)=[];LiverFreqFastCORE(LiverFreqFastCORE<0)=[];
LiverCorExRxnRanking=[{'Frequency in AO','Rxn Name','Rxn SubSystem'};[num2cell(LiverFreqCorEx),Recon1red.rxnNames(idx1),Recon1red.subSystems(idx1)];{'Fixed Active','Fixed Inactive','Variable'};num2cell([length(find(LiverFreqCorEx==1)),length(find(LiverFreqCorEx==0)),length(find(LiverFreqCorEx~=1 & LiverFreqCorEx~=0))])];
LiverFastCORERxnRanking=[{'Frequency in AO','Rxn Name','Rxn SubSystem'};[num2cell(LiverFreqFastCORE),Recon1red.rxnNames(idx2),Recon1red.subSystems(idx2)];{'Fixed Active','Fixed Inactive','Variable'};num2cell([length(find(LiverFreqFastCORE==1)),length(find(LiverFreqFastCORE==0)),length(find(LiverFreqFastCORE~=1 & LiverFreqFastCORE~=0))])];

%get the distribution of reactions in subsystems in the fixed and the
%variable non-core set
B=[];C=[];D=[];E=[];
NCSubSys=unique(LiverCorExRxnRanking(2:(end-2),3));
for i=1:length(NCSubSys),
    E(i,1)=sum(strcmp(LiverCorExRxnRanking(2:(end-2),3),NCSubSys{i}));
end

ActiveNCSubSys=LiverCorExRxnRanking(find(cell2mat(LiverCorExRxnRanking(2:(end-2),1))==1)+1,3);
VarNCSubSys=LiverCorExRxnRanking(find(cell2mat(LiverCorExRxnRanking(2:(end-2),1))<1 & cell2mat(LiverCorExRxnRanking(2:(end-2),1))>0)+1,3);
InactiveNCSubSys=LiverCorExRxnRanking(find(cell2mat(LiverCorExRxnRanking(2:(end-2),1))==0)+1,3);

UniqActiveNCSubSys=unique(ActiveNCSubSys);UniqVarNCSubSys=unique(VarNCSubSys);UniqInactiveNCSubSys=unique(InactiveNCSubSys);
for i=1:length(UniqInactiveNCSubSys),
    B(i,1)=sum(strcmp(InactiveNCSubSys,UniqInactiveNCSubSys{i}))/E(strcmp(NCSubSys,UniqInactiveNCSubSys{i}));
end
for i=1:length(UniqActiveNCSubSys),
    C(i,1)=sum(strcmp(ActiveNCSubSys,UniqActiveNCSubSys{i}))/E(strcmp(NCSubSys,UniqActiveNCSubSys{i}));
end
for i=1:length(UniqVarNCSubSys),
    D(i,1)=sum(strcmp(VarNCSubSys,UniqVarNCSubSys{i}))/E(strcmp(NCSubSys,UniqVarNCSubSys{i}));
end
[B,idx2]=sort(B,'descend');[C,idx3]=sort(C,'descend');[D,idx4]=sort(D,'descend');
LiverCorExInactSubSys=[UniqInactiveNCSubSys(idx2),num2cell(B)];
LiverCorExActSubSys=[UniqActiveNCSubSys(idx3),num2cell(C)];
LiverCorExVarSubSys=[UniqVarNCSubSys(idx4),num2cell(D)];

B=[];C=[];D=[];E=[];
NCSubSys=unique(LiverFastCORERxnRanking(2:(end-2),3));
for i=1:length(NCSubSys),
    E(i,1)=sum(strcmp(LiverFastCORERxnRanking(2:(end-2),3),NCSubSys{i}));
end

ActiveNCSubSys=LiverFastCORERxnRanking(find(cell2mat(LiverFastCORERxnRanking(2:(end-2),1))==1)+1,3);
VarNCSubSys=LiverFastCORERxnRanking(find(cell2mat(LiverFastCORERxnRanking(2:(end-2),1))<1 & cell2mat(LiverFastCORERxnRanking(2:(end-2),1))>0)+1,3);
InactiveNCSubSys=LiverFastCORERxnRanking(find(cell2mat(LiverFastCORERxnRanking(2:(end-2),1))==0)+1,3);

UniqActiveNCSubSys=unique(ActiveNCSubSys);UniqVarNCSubSys=unique(VarNCSubSys);UniqInactiveNCSubSys=unique(InactiveNCSubSys);
for i=1:length(UniqInactiveNCSubSys),
    B(i,1)=sum(strcmp(InactiveNCSubSys,UniqInactiveNCSubSys{i}))/E(strcmp(NCSubSys,UniqInactiveNCSubSys{i}));
end
for i=1:length(UniqActiveNCSubSys),
    C(i,1)=sum(strcmp(ActiveNCSubSys,UniqActiveNCSubSys{i}))/E(strcmp(NCSubSys,UniqActiveNCSubSys{i}));
end
for i=1:length(UniqVarNCSubSys),
    D(i,1)=sum(strcmp(VarNCSubSys,UniqVarNCSubSys{i}))/E(strcmp(NCSubSys,UniqVarNCSubSys{i}));
end
[B,idx2]=sort(B,'descend');[C,idx3]=sort(C,'descend');[D,idx4]=sort(D,'descend');
LiverFastCOREInactSubSys=[UniqInactiveNCSubSys(idx2),num2cell(B)];
LiverFastCOREActSubSys=[UniqActiveNCSubSys(idx3),num2cell(C)];
LiverFastCOREVarSubSys=[UniqVarNCSubSys(idx4),num2cell(D)];

%calculate hamming distance
HmatLiverCorEx=HammingMat(LiverCorExAO.Modmatrix);
HmatLiverFastCORE=HammingMat(LiverFastCOREAO.Modmatrix);

%get distribution of hamming distances 
HDisLiverCorEx = HmatLiverCorEx(triu(logical(ones(size(HmatLiverCorEx,1))),1));
HDisLiverFastCORE = HmatLiverFastCORE(triu(logical(ones(size(HmatLiverFastCORE,1))),1));

%compare Hamming distance distributions with ranksum test
[HpLeaf,~]=ranksum(HDisLeafCorEx,HDisLeafFastCORE,'tail','left');
[HpLiver,~]=ranksum(HDisLiverCorEx,HDisLiverFastCORE,'tail','left');
% [JpLeaf,~]=ranksum(JDisLeafCorEx,JDisLeafFastCORE,'tail','left');
% [JpLiver,~]=ranksum(JDisLiverCorEx,JDisLiverFastCORE,'tail','left');

%get CorEx Pie charts
Pielabels = {'Active Non-Core','Inactive Non-Core','Variable Non-Core'};
LeafCore=length(CLeaf);
LeafCorExFixed=length(find(LeafCorExtotnoC==round(max(LeafCorExtotnoC))));
LeafCorExVar=length(find(LeafCorExtotnoC>0 & LeafCorExtotnoC<round(max(LeafCorExtotnoC))));
LeafFastCOREFixed=length(find(LeafFastCOREtotnoC==round(max(LeafFastCOREtotnoC))));
LeafFastCOREVar=length(find(LeafFastCOREtotnoC>0 & LeafFastCOREtotnoC<round(max(LeafFastCOREtotnoC))));
LiverCore=length(CLiver);
LiverCorExFixed=length(find(LiverCorExtotnoC==round(max(LiverCorExtotnoC))));
LiverCorExVar=length(find(LiverCorExtotnoC>0 & LiverCorExtotnoC<round(max(LiverCorExtotnoC))));
LiverFastCOREFixed=length(find(LiverFastCOREtotnoC==round(max(LiverFastCOREtotnoC))));
LiverFastCOREVar=length(find(LiverFastCOREtotnoC>0 & LiverFastCOREtotnoC<round(max(LiverFastCOREtotnoC))));

figure(1);pie([LeafCorExFixed,size(AraCOREred.S,2)-(LeafCorExFixed+LeafCorExVar+LeafCore),LeafCorExVar]);title('LeafCorEx');colormap([51/256,153/256,255/256;0.1,0.85,0.15;1,1,0])
figure(2);pie([LeafFastCOREFixed,size(AraCOREred.S,2)-(LeafFastCOREFixed+LeafFastCOREVar+LeafCore),LeafFastCOREVar]);title('LeafFastCORE');colormap([51/256,153/256,255/256;0.1,0.85,0.15;1,1,0])
figure(3);pie([LiverCorExFixed,size(Recon1red.S,2)-(LiverCorExFixed+LiverCorExVar+LiverCore),LiverCorExVar]);title('LiverCorEx');colormap([51/256,153/256,255/256;0.1,0.85,0.15;1,1,0])
figure(4);pie([LiverFastCOREFixed,size(Recon1red.S,2)-(LiverFastCOREFixed+LiverFastCOREVar+LiverCore),LiverFastCOREVar]);title('LiverFastCORE');legend(Pielabels);colormap([51/256,153/256,255/256;0.1,0.85,0.15;1,1,0])

%**************************************************************************
%***************************CORDA AO Analysis*****************************
%**************************************************************************

%Obtain reaction indexes in Recon1red of the HC, MC, NC and OT groups of
%reactions
HCrecon1red=find(ismember(Recon1red.rxns,HC)==1);
MCrecon1red=find(ismember(Recon1red.rxns,MC)==1);
NCrecon1red=find(ismember(Recon1red.rxns,NC)==1);
OTrecon1red=find(ismember(Recon1red.rxns,OT)==1);OTrecon1red=[OTrecon1red;setdiff(1:length(Recon1red.rxns),[HCrecon1red;MCrecon1red;NCrecon1red;OTrecon1red])'];

%Obtain liver CorExCORDA models and evaluate metabolic tasks
LiverCORDA=CorExCORDA(Recon1red,HCrecon1red,MCrecon1red,NCrecon1red,OTrecon1red,279,369,11,1147,1000,1e-3,5*60)

%Use original liverCORDAnew model
load('RESULTSCorExCORDA.mat','liverCORDAnew')
LiverCORDAIdx=find(ismember(Recon1red.rxns,liverCORDAnew.rxns)==1);
LiverCORDAAO=AltNetCORDA(Recon1red,HCrecon1red,MCrecon1red,NCrecon1red,OTrecon1red,369,11,1147,LiverCORDAIdx,1000,1e-3,200,5*60)

%calculate hamming distance
HmatLiverCORDA=HammingMat(LiverCORDAAO.Modmatrix);
HmatLiverCORDATest=HammingMat(LiverCORDAAO.ModmatrixTest);

%get distribution of hamming distances 
HDisLiverCORDA=HmatLiverCORDA(triu(logical(ones(size(HmatLiverCORDA,1))),1));
HDisLiverCORDAtest=HmatLiverCORDATest(triu(logical(ones(size(HmatLiverCORDATest,1))),1));

%get the fixed and variable fractions of the total number of non-core
%reactions among alternative networks
LiverCORDAtotnoC=LiverCORDAAO.Modmatrix;LiverCORDAtotnoC(HCrecon1red,:)=-1;LiverCORDAtotnoC=sum(LiverCORDAtotnoC')';FliverCORDAnoC=find(LiverCORDAtotnoC==size(LiverCORDAAO.Modmatrix,2));VliverCORDAnoC=find(LiverCORDAtotnoC>0 & LiverCORDAtotnoC<size(LiverCORDAAO.Modmatrix,2));
LiverCORDATesttotnoC=LiverCORDAAO.ModmatrixTest;LiverCORDATesttotnoC(HCrecon1red,:)=-1;LiverCORDATesttotnoC=sum(LiverCORDATesttotnoC')';FliverCORDATestnoC=find(LiverCORDATesttotnoC==size(LiverCORDAAO.ModmatrixTest,2));VliverCORDATestnoC=find(LiverCORDATesttotnoC>0 & LiverCORDATesttotnoC<size(LiverCORDAAO.ModmatrixTest,2));

%rank reactions according to its frequency of ocurrence in the AO space of
%models
[LiverFreqCORDA,idx1]=sort(LiverCORDAtotnoC/max(LiverCORDAtotnoC),'descend');
[LiverFreqCORDATest,idx2]=sort(LiverCORDATesttotnoC/max(LiverCORDATesttotnoC),'descend');
idx1(LiverFreqCORDA<0)=[];LiverFreqCORDA(LiverFreqCORDA<0)=[];
idx2(LiverFreqCORDATest<0)=[];LiverFreqCORDATest(LiverFreqCORDATest<0)=[];
LiverCORDARxnRanking=[{'Frequency in AO','Rxn Name','Rxn SubSystem'};[num2cell(LiverFreqCORDA),Recon1red.rxnNames(idx1),Recon1red.subSystems(idx1)];{'Fixed Active','Fixed Inactive','Variable'};num2cell([length(find(LiverFreqCORDA==1)),length(find(LiverFreqCORDA==0)),length(find(LiverFreqCORDA~=1 & LiverFreqCORDA~=0))])];
LiverCORDATestRxnRanking=[{'Frequency in AO','Rxn Name','Rxn SubSystem'};[num2cell(LiverFreqCORDATest),Recon1red.rxnNames(idx2),Recon1red.subSystems(idx2)];{'Fixed Active','Fixed Inactive','Variable'};num2cell([length(find(LiverFreqCORDATest==1)),length(find(LiverFreqCORDATest==0)),length(find(LiverFreqCORDATest~=1 & LiverFreqCORDATest~=0))])];

%get the distribution of reactions in subsystems in the fixed and the
%variable non-core set
B=[];C=[];D=[];E=[];
NCSubSys=unique(LiverCORDARxnRanking(2:(end-2),3));
for i=1:length(NCSubSys),
    E(i,1)=sum(strcmp(LiverCORDARxnRanking(2:(end-2),3),NCSubSys{i}));
end

ActiveNCSubSys=LiverCORDARxnRanking(find(cell2mat(LiverCORDARxnRanking(2:(end-2),1))==1)+1,3);
VarNCSubSys=LiverCORDARxnRanking(find(cell2mat(LiverCORDARxnRanking(2:(end-2),1))<1 & cell2mat(LiverCORDARxnRanking(2:(end-2),1))>0)+1,3);
InactiveNCSubSys=LiverCORDARxnRanking(find(cell2mat(LiverCORDARxnRanking(2:(end-2),1))==0)+1,3);

UniqActiveNCSubSys=unique(ActiveNCSubSys);UniqVarNCSubSys=unique(VarNCSubSys);UniqInactiveNCSubSys=unique(InactiveNCSubSys);
for i=1:length(UniqInactiveNCSubSys),
    B(i,1)=sum(strcmp(InactiveNCSubSys,UniqInactiveNCSubSys{i}))/E(strcmp(NCSubSys,UniqInactiveNCSubSys{i}));
end
for i=1:length(UniqActiveNCSubSys),
    C(i,1)=sum(strcmp(ActiveNCSubSys,UniqActiveNCSubSys{i}))/E(strcmp(NCSubSys,UniqActiveNCSubSys{i}));
end
for i=1:length(UniqVarNCSubSys),
    D(i,1)=sum(strcmp(VarNCSubSys,UniqVarNCSubSys{i}))/E(strcmp(NCSubSys,UniqVarNCSubSys{i}));
end
[B,idx2]=sort(B,'descend');[C,idx3]=sort(C,'descend');[D,idx4]=sort(D,'descend');
LiverCORDAInactSubSys=[UniqInactiveNCSubSys(idx2),num2cell(B)];
LiverCORDAActSubSys=[UniqActiveNCSubSys(idx3),num2cell(C)];
LiverCORDAVarSubSys=[UniqVarNCSubSys(idx4),num2cell(D)];

%Compare Hamming distances in all models and the subset passing metabolic
%test
[HpLiverCorda,~]=ranksum(HDisLiverCORDAtest,HDisLiverCORDA,'tail','left');

%get CorEx (CORDA) Pie charts (Total and alternative networks with equal or
%greater metabolic test score result)
Pielabels = {'Active Non-Core','Inactive Non-Core','Variable Non-Core'};
LiverHCore=length(HCrecon1red);
LiverCORDAFixed=length(find(LiverCORDAtotnoC==round(max(LiverCORDAtotnoC))));
LiverCORDAVar=length(find(LiverCORDAtotnoC>0 & LiverCORDAtotnoC<round(max(LiverCORDAtotnoC))));
LiverCORDATestFixed=length(find(LiverCORDATesttotnoC==round(max(LiverCORDATesttotnoC))));
LiverCORDATestVar=length(find(LiverCORDATesttotnoC>0 & LiverCORDATesttotnoC<round(max(LiverCORDATesttotnoC))));

figure(5);pie([LiverCORDAFixed,size(Recon1red.S,2)-(LiverCORDAFixed+LiverCORDAVar+LiverHCore),LiverCORDAVar]);title('LiverCORDA');legend(Pielabels);colormap([51/256,153/256,255/256;0.1,0.85,0.15;1,1,0])
figure(6);pie([LiverCORDATestFixed,size(Recon1red.S,2)-(LiverCORDATestFixed+LiverCORDATestVar+LiverHCore),LiverCORDATestVar]);title('LiverCORDAtest');legend(Pielabels);colormap([51/256,153/256,255/256;0.1,0.85,0.15;1,1,0])

%get Network-centered Table:
ColNames={'Case study','P','#models','MRmax','MRmean','MRCV','p-value'};
RowNames={'LeafCorex';'LeafFastCORE';'LiverCorEx';'LiverFastCORE';'LiverCORDA';'LiverCORDAtest'};
C1=[LeafCorEx.AddedRxns;PleafFast;LiverCorEx.AddedRxns;PliverFast;369+11+1147;369+11+1147];
C2=[size(LeafCorExAO.Modmatrix,2);size(LeafFastCOREAO.Modmatrix,2);size(LiverCorExAO.Modmatrix,2);size(LiverFastCOREAO.Modmatrix,2);size(LiverCORDAAO.Modmatrix,2);size(LiverCORDAAO.ModmatrixTest,2)];
C3=[max(HDisLeafCorEx);max(HDisLeafFastCORE);max(HDisLiverCorEx);max(HDisLiverFastCORE);max(HDisLiverCORDA);max(HDisLiverCORDAtest)];
C4=[mean(HDisLeafCorEx);mean(HDisLeafFastCORE);mean(HDisLiverCorEx);mean(HDisLiverFastCORE);mean(HDisLiverCORDA);mean(HDisLiverCORDAtest)];
C5=[std(HDisLeafCorEx);std(HDisLeafFastCORE);std(HDisLiverCorEx);std(HDisLiverFastCORE);std(HDisLiverCORDA);std(HDisLiverCORDAtest)]./C4;
C6={num2str(HpLeaf);'';num2str(HpLiver);'';num2str(HpLiverCorda);''};
NetworkTable=[ColNames;[RowNames,num2cell([C1,C2,C3,C4,C5]),C6]];

%**************************************************************************
%***************************iMAT AO Analysis*****************************
%**************************************************************************

%Obtain iMAT model and analyze alternative optima (iMAT FVA-like approach)
LeafiMAT=iMAT(AraCOREred,[],CLeaf,1e-6,'LeafiMAT');

%Approach proposed in the original iMAT publication
LeafiMATFVA=iMATFVA(AraCOREred,[],CLeaf,1e-6);
LeafiMATFVA(:,3)=LeafiMATFVA(:,1)-LeafiMATFVA(:,2);
LeafiMATFVA(abs(LeafiMAT.Flux)>=1e-6,4)=1;

%Approach proposed in this study (AO sample)
LeafiMATAO=iMATAO(AraCOREred,LeafiMAT.Zopt,[],CLeaf,1e-6,2000,Options00.FVArange);
%calculate hamming distance
HmatLeafiMATAO=HammingMat(LeafiMATAO.Modmatrix);
%get distribution of hamming distances 
HDisLeafiMATAO = HmatLeafiMATAO(triu(logical(ones(size(HmatLeafiMATAO,1))),1));[mean(HDisLeafiMATAO),std(HDisLeafiMATAO)/mean(HDisLeafiMATAO)]
%get the fixed and variable fractions of the total number 
%reactions among alternative networks and rank reactions according to their
%frequency of occurrence
[FreqLeafiMATAO(:,1),FreqLeafiMATAO(:,2)]=sort(sum(LeafiMATAO.Modmatrix')'/max(sum(LeafiMATAO.Modmatrix')),'descend');

%Get table of ranked reactions according to their frequency in the AO space
%of iMAT
LeafiMATRxnRanking=[{'Frequency in AO','Rxn index','Rxn Name','Rxn SubSystem'};[num2cell(FreqLeafiMATAO(:,1)),num2cell(FreqLeafiMATAO(:,2)),AraCOREred.rxnNames(FreqLeafiMATAO(:,2)),AraCOREred.subSystems(FreqLeafiMATAO(:,2))]];

length(intersect(find(LeafiMATFVA(:,3)<0),FreqLeafiMATAO(FreqLeafiMATAO(:,1)==1,2)))/length(find(LeafiMATFVA(:,3)<0))%fraction of active reactions per iMATFVA in iMATAO sample
length(intersect(find(LeafiMATFVA(:,3)>0),FreqLeafiMATAO(FreqLeafiMATAO(:,1)==0,2)))/length(find(LeafiMATFVA(:,3)>0))%fraction of inactive reactions per iMATFVA in iMATAO sample
length(intersect(find(LeafiMATFVA(:,3)==0),FreqLeafiMATAO(FreqLeafiMATAO(:,1)<1 & FreqLeafiMATAO(:,1)>0,2)))/length(find(LeafiMATFVA(:,3)==0))%fraction of undetermined reactions per iMATFVA in iMATAO sample

%Obtain iMAT model and analyze alternative optima
LiveriMAT=iMAT(Recon1red,[],CLiver,1e-6,'LiveriMAT');
%Approach proposed in the original iMAT publication
LiveriMATFVA=iMATFVA(Recon1red,[],CLiver,1e-6);
LiveriMATFVA(:,3)=LiveriMATFVA(:,1)-LiveriMATFVA(:,2);
LiveriMATFVA(abs(LiveriMAT.Flux)>=1e-6,4)=1;
length(find(LiveriMATFVA(:,3)==0)) %number of undetermined across AO space
length(find(LiveriMATFVA(:,3)>0)) %number of inactive across AO space
length(find(LiveriMATFVA(:,3)<0)) %number of active across AO space

%Approach proposed in this study (AO sample)
LiveriMATAO=iMATAO(Recon1red,LiveriMAT.Zopt,[],CLiver,1e-6,2000,Recon1FVArange);
%calculate hamming distance
HmatLiveriMATAO=HammingMat(LiveriMATAO.Modmatrix);
%get distribution of hamming distances 
HDisLiveriMATAO = HmatLiveriMATAO(triu(logical(ones(size(HmatLiveriMATAO,1))),1));[mean(HDisLiveriMATAO),std(HDisLiveriMATAO)/mean(HDisLiveriMATAO)]
%get the fixed and variable fractions of the total number 
%reactions among alternative networks and rank reactions according to their
%frequency of occurrence
[FreqLiveriMATAO(:,1),FreqLiveriMATAO(:,2)]=sort(sum(LiveriMATAO.Modmatrix')'/max(sum(LiveriMATAO.Modmatrix')),'descend');

length(intersect(find(LiveriMATFVA(:,3)<0),FreqLiveriMATAO(FreqLiveriMATAO(:,1)==1,2)))/length(find(LiveriMATFVA(:,3)<0))%fraction of active reactions per iMATFVA in iMATAO sample
length(intersect(find(LiveriMATFVA(:,3)>0),FreqLiveriMATAO(FreqLiveriMATAO(:,1)==0,2)))/length(find(LiveriMATFVA(:,3)>0))%fraction of inactive reactions per iMATFVA in iMATAO sample
length(intersect(find(LiveriMATFVA(:,3)==0),FreqLiveriMATAO(FreqLiveriMATAO(:,1)<1 & FreqLiveriMATAO(:,1)>0,2)))/length(find(LiveriMATFVA(:,3)==0))%fraction of undetermined reactions per iMATFVA in iMATAO sample

 %ranksum test to compare Leaf and Liver Hamming distance distributions (normalized by
 %total number of reactions)
[pHDisiMAT,~]=ranksum(HDisLiveriMATAO/size(Recon1red.S,2),HDisLeafiMATAO/size(AraCOREred.S,2),'tail','right')

%**************************************************************************
%***************************RegrEx AO Analysis*****************************
%**************************************************************************

%******************************Leaf case***********************************

%get RegrExLAD solutions
[Options.FVArange(:,1),Options.FVArange(:,2)]=FVA(AraCOREred);
Options.Nsamples=10;Options.samplesize=200;

Options00=Options;Options00.L_min=0;Options00.L_int=0;Options00.L_max=0;
Options01=Options;Options01.L_min=0.1;Options01.L_int=0.1;Options01.L_max=0.1;
Options03=Options;Options03.L_min=0.3;Options03.L_int=0.3;Options03.L_max=0.3;
Options05=Options;Options05.L_min=0.5;Options05.L_int=0.5;Options05.L_max=0.5;

Leaf00=RegrExLAD(AraCOREred,LeafData(:,1),Options00)
Leaf01=RegrExLAD(AraCOREred,LeafData(:,1),Options01)
Leaf03=RegrExLAD(AraCOREred,LeafData(:,1),Options03)
Leaf05=RegrExLAD(AraCOREred,LeafData(:,1),Options05)

%explore RegrExLAD alternative optima
Leaf00AO=RegrExAO(AraCOREred,LeafData(:,1),Leaf00,Options00)
Leaf01AO=RegrExAO(AraCOREred,LeafData(:,1),Leaf01,Options01)
Leaf03AO=RegrExAO(AraCOREred,LeafData(:,1),Leaf03,Options03)
Leaf05AO=RegrExAO(AraCOREred,LeafData(:,1),Leaf05,Options05)

%extract alternative RegrEx networks
Leaf00Anet=zeros(size(Leaf00AO.Vsample));Leaf00Anet(Leaf00AO.Vsample~=0)=1;Leaf00Anet=unique(Leaf00Anet','rows')';
Leaf01Anet=zeros(size(Leaf01AO.Vsample));Leaf01Anet(Leaf01AO.Vsample~=0)=1;Leaf01Anet=unique(Leaf01Anet','rows')';
Leaf03Anet=zeros(size(Leaf03AO.Vsample));Leaf03Anet(Leaf03AO.Vsample~=0)=1;Leaf03Anet=unique(Leaf03Anet','rows')';
Leaf05Anet=zeros(size(Leaf05AO.Vsample));Leaf05Anet(Leaf05AO.Vsample~=0)=1;Leaf05Anet=unique(Leaf05Anet','rows')';

%get the fixed and variable fractions of the total number
%reactions among alternative networks and rank reactions according to their
%frequency of occurrence
[FreqLeaf00Anet(:,1),FreqLeaf00Anet(:,2)]=sort(sum(Leaf00Anet')'/max(sum(Leaf00Anet')),'descend');
[FreqLeaf01Anet(:,1),FreqLeaf01Anet(:,2)]=sort(sum(Leaf01Anet')'/max(sum(Leaf01Anet')),'descend');
[FreqLeaf03Anet(:,1),FreqLeaf03Anet(:,2)]=sort(sum(Leaf03Anet')'/max(sum(Leaf03Anet')),'descend');
[FreqLeaf05Anet(:,1),FreqLeaf05Anet(:,2)]=sort(sum(Leaf05Anet')'/max(sum(Leaf05Anet')),'descend');

%get s2table
s2Table=[{'Frequency in AO space', 'Rxn Idx', 'Rxn Name', 'Subsystem'};[num2cell(FreqLeaf00Anet),AraCOREred.rxnNames(FreqLeaf00Anet(:,2)),AraCOREred.subSystems(FreqLeaf00Anet(:,2))]];
xlswrite('S2Table.xlsx',s2Table)

%calculate hamming distance
HmatLeaf00=HammingMat(Leaf00Anet);
HmatLeaf01=HammingMat(Leaf01Anet);
HmatLeaf03=HammingMat(Leaf03Anet);
HmatLeaf05=HammingMat(Leaf05Anet);

%get distribution of hamming distances
HDis00 = HmatLeaf00(triu(logical(ones(size(HmatLeaf00,1))),1));mean(HDis00)
HDis01 = HmatLeaf01(triu(logical(ones(size(HmatLeaf01,1))),1));
HDis03 = HmatLeaf03(triu(logical(ones(size(HmatLeaf03,1))),1));
HDis05 = HmatLeaf05(triu(logical(ones(size(HmatLeaf05,1))),1));
meanHDis=[mean(HDis00),mean(HDis01),mean(HDis03),mean(HDis05)];
CVHDis=[std(HDis00),std(HDis01),std(HDis03),std(HDis05)]./meanHDis;
DataRxns=find(LeafData(:,1)>0);
DataOrphan=setdiff(1:length(AraCOREred.rxns),DataRxns);

%calculate entropies of RegrEx AO solution
LeafEntro(:,1)=DistEntropy(Leaf00AO.Vsample,[],'F','F');
LeafEntro(:,2)=DistEntropy(Leaf01AO.Vsample,[],'F','F');
LeafEntro(:,3)=DistEntropy(Leaf03AO.Vsample,[],'F','F');
LeafEntro(:,4)=DistEntropy(Leaf05AO.Vsample,[],'F','F');
EntroDataOrphan=LeafEntro(intersect(DataOrphan,find(sum(abs(Leaf00AO.Vsample)')>0)),:);
EntroDataBounded=LeafEntro(intersect(DataRxns,find(sum(abs(Leaf00AO.Vsample)')>0)),:);
mean(EntroDataOrphan(:,1))  %mean entropy in DataOrhpan
mean(EntroDataBounded(:,1)) %mean entropy in DataBounded
[pLeafOrphanBounded,~]=ranksum(EntroDataOrphan(:,1),EntroDataBounded(:,1),'tail','right')

%find reversible reactions with fixed direction
% 0 -> fixed reverse direction
% Inf -> fixed forward direction
% NaN -> inactive reaction
% r -> forward/reverse ratio
RevRxns=find(AraCOREred.rev==1);
for i=1:length(RevRxns),
    FixedRev(i,1)=length(find(Leaf00AO.Vsample(RevRxns(i),:)>eps))/length(find(Leaf00AO.Vsample(RevRxns(i),:)<-eps));
    FixedRev(i,2)=length(find(Leaf01AO.Vsample(RevRxns(i),:)>eps))/length(find(Leaf01AO.Vsample(RevRxns(i),:)<-eps));
    FixedRev(i,3)=length(find(Leaf03AO.Vsample(RevRxns(i),:)>eps))/length(find(Leaf03AO.Vsample(RevRxns(i),:)<-eps));
    FixedRev(i,4)=length(find(Leaf05AO.Vsample(RevRxns(i),:)>eps))/length(find(Leaf05AO.Vsample(RevRxns(i),:)<-eps));
end
for i=1:4,
ProportionFixedRev(i)=length(find(FixedRev(:,i)==0 | isinf(FixedRev(:,i)))) / length(find(~isnan(FixedRev(:,i))));
end

%calculate AraCORE reference entropies without data integration
AraCOREredsample=coneSampling(AraCOREred,'F',Options00.FVArange,[],2000,[],'gurobi')
AraCOREredEntro=DistEntropy(AraCOREredsample.Vsample,[],'F','F');

%compare entropies
colnames={'Entropy(RegrEx)','Entropy(flux cone)','Rxn Idx','Data status','Rxn Name','Subsystem','Mechanism'};
[sortentro,sortidx]=sort(LeafEntro(:,1),'descend');
sortidx((sortentro==0))=[];sortentro(sortentro==0)=[];
datastatus=zeros(length(AraCOREred.rxns),1);datastatus(DataRxns)=1;
s1Table=[colnames;[num2cell([sortentro,AraCOREredEntro(sortidx),sortidx, datastatus(sortidx)]),AraCOREred.rxnNames(sortidx),AraCOREred.subSystems(sortidx),AraCOREred.rxnMechanisms(sortidx)]];
xlswrite('S1Table.xlsx',s1Table)

%get RegrEx alternative optima table
colnames={'L0','L0.1','L0.3','L0.5'};
rownames={'RegrExAO';'H_D';'H_nD';'Htotal';'median H_D';'median H_nD';'median Htotal';'F_rev';'D_mean';'D_CV'};
RegrExAOtable=[rownames,[colnames;num2cell([sum(EntroDataBounded);sum(EntroDataOrphan);sum(EntroDataBounded)+sum(EntroDataOrphan);mean(EntroDataBounded);mean(EntroDataOrphan);mean([EntroDataBounded;EntroDataOrphan]);ProportionFixedRev;meanHDis;CVHDis],3)]];

%get figures
figure(7)
X=sort(LeafEntro(DataOrphan,:),'descend');
xcutoff=find(X<1e-3,1);
stairs(X(1:xcutoff,:),'LineWidth',1.5)
ylabel('H_s')
xlabel('rxns')
legend({'\lambda = 0','\lambda = 0.1','\lambda = 0.2','\lambda = 0.5'})
title('Data-orphan reactions')
figure(8)
X=sort(LeafEntro(DataRxns,:),'descend');
xcutoff=find(X<1e-3,1);
stairs(X(1:xcutoff,:),'LineWidth',1.5)
ylabel('H_s')
xlabel('rxns')
legend({'\lambda = 0','\lambda = 0.1','\lambda = 0.2','\lambda = 0.5'})
title('Data-bounded reactions')
figure(9)
boxplot(EntroDataOrphan)
set(gca,'XTickLabel',{'0','0.1','0.2','0.5'},'fontsize',16)
title('Data-orphan reactions')
figure(10)
boxplot(EntroDataBounded)
set(gca,'XTickLabel',{'0','0.1','0.2','0.5'},'fontsize',16)
title('Data-bounded reactions')

%rank-sum tests
[pDataBounded(1,1),~] = ranksum(EntroDataBounded(:,1),EntroDataBounded(:,2),'tail','right');
[pDataBounded(2,1),~] = ranksum(EntroDataBounded(:,2),EntroDataBounded(:,3),'tail','right');
[pDataBounded(3,1),~] = ranksum(EntroDataBounded(:,3),EntroDataBounded(:,4),'tail','right');
[pDataOrphan(1,1),~] = ranksum(EntroDataOrphan(:,1),EntroDataOrphan(:,2),'tail','right');
[pDataOrphan(2,1),~] = ranksum(EntroDataOrphan(:,2),EntroDataOrphan(:,3),'tail','right');
[pDataOrphan(3,1),~] = ranksum(EntroDataOrphan(:,3),EntroDataOrphan(:,4),'tail','right');

%RegrExOLS vs RegrExLAD on data contaminated with outliers
OLSLADtestReconGaussian=RegrExOLSLADstudy(AraCOREred,LeafData(:,1),1e4)

%*****************************Liver case***********************************

%get RegrExLAD solutions
[Options.FVArange(:,1),Options.FVArange(:,2)]=FVA(Recon1red);
Options.Nsamples=10;Options.samplesize=200;

Options00=Options;Options00.L_min=0;Options00.L_int=0;Options00.L_max=0;
Options01=Options;Options01.L_min=0.1;Options01.L_int=0.1;Options01.L_max=0.1;
Options03=Options;Options03.L_min=0.3;Options03.L_int=0.3;Options03.L_max=0.3;
Options05=Options;Options05.L_min=0.5;Options05.L_int=0.5;Options05.L_max=0.5;

Liver00=RegrExLAD(Recon1red,LiverData(:,1),Options00)
Liver01=RegrExLAD(Recon1red,LiverData(:,1),Options01)
Liver03=RegrExLAD(Recon1red,LiverData(:,1),Options03)
Liver05=RegrExLAD(Recon1red,LiverData(:,1),Options05)

%explore RegrExLAD alternative optima
Liver00AO=RegrExAO(Recon1red,LiverData(:,1),Liver00,Options00)
Liver01AO=RegrExAO(Recon1red,LiverData(:,1),Liver01,Options01)
Liver03AO=RegrExAO(Recon1red,LiverData(:,1),Liver03,Options03)
Liver05AO=RegrExAO(Recon1red,LiverData(:,1),Liver05,Options05)

%extract alternative RegrEx networks
Liver00Anet=zeros(size(Liver00AO.Vsample));Liver00Anet(Liver00AO.Vsample~=0)=1;Liver00Anet=unique(Liver00Anet','rows')';
Liver01Anet=zeros(size(Liver01AO.Vsample));Liver01Anet(Liver01AO.Vsample~=0)=1;Liver01Anet=unique(Liver01Anet','rows')';
Liver03Anet=zeros(size(Liver03AO.Vsample));Liver03Anet(Liver03AO.Vsample~=0)=1;Liver03Anet=unique(Liver03Anet','rows')';
Liver05Anet=zeros(size(Liver05AO.Vsample));Liver05Anet(Liver05AO.Vsample~=0)=1;Liver05Anet=unique(Liver05Anet','rows')';

%get the fixed and variable fractions of the total number
%reactions among alternative networks and rank reactions according to their
%frequency of occurrence
[FreqLiver00Anet(:,1),FreqLiver00Anet(:,2)]=sort(sum(Liver00Anet')'/max(sum(Liver00Anet')),'descend');
[FreqLiver01Anet(:,1),FreqLiver01Anet(:,2)]=sort(sum(Liver01Anet')'/max(sum(Liver01Anet')),'descend');
[FreqLiver03Anet(:,1),FreqLiver03Anet(:,2)]=sort(sum(Liver03Anet')'/max(sum(Liver03Anet')),'descend');
[FreqLiver05Anet(:,1),FreqLiver05Anet(:,2)]=sort(sum(Liver05Anet')'/max(sum(Liver05Anet')),'descend');

% %get s2table
 s2LiverTable=[{'Frequency in AO space', 'Rxn Idx', 'Rxn Name', 'Subsystem'};[num2cell(FreqLiver00Anet),Recon1red.rxnNames(FreqLiver00Anet(:,2)),Recon1red.subSystems(FreqLiver00Anet(:,2))]];

%calculate hamming distance
HmatLiver00=HammingMat(Liver00Anet);
HmatLiver01=HammingMat(Liver01Anet);
HmatLiver03=HammingMat(Liver03Anet);
HmatLiver05=HammingMat(Liver05Anet);

%get distribution of hamming distances
HDis00 = HmatLiver00(triu(logical(ones(size(HmatLiver00,1))),1));mean(HDis00)
HDis01 = HmatLiver01(triu(logical(ones(size(HmatLiver01,1))),1));
HDis03 = HmatLiver03(triu(logical(ones(size(HmatLiver03,1))),1));
HDis05 = HmatLiver05(triu(logical(ones(size(HmatLiver05,1))),1));
meanHDis=[mean(HDis00),mean(HDis01),mean(HDis03),mean(HDis05)];
CVHDis=[std(HDis00),std(HDis01),std(HDis03),std(HDis05)]./meanHDis;
DataRxns=find(LiverData(:,1)>0);
DataOrphan=setdiff(1:length(Recon1red.rxns),DataRxns);

%calculate entropies of RegrEx AO solution
LiverEntro(:,1)=DistEntropy(Liver00AO.Vsample,[],'F','F');
LiverEntro(:,2)=DistEntropy(Liver01AO.Vsample,[],'F','F');
LiverEntro(:,3)=DistEntropy(Liver03AO.Vsample,[],'F','F');
LiverEntro(:,4)=DistEntropy(Liver05AO.Vsample,[],'F','F');
EntroDataOrphan=LiverEntro(intersect(DataOrphan,find(sum(abs(Liver00AO.Vsample)')>0)),:);
EntroDataBounded=LiverEntro(intersect(DataRxns,find(sum(abs(Liver00AO.Vsample)')>0)),:);
mean(EntroDataOrphan(:,1))  %mean entropy in DataOrhpan
mean(EntroDataBounded(:,1)) %mean entropy in DataBounded
[pLiverOrphanBounded,~]=ranksum(EntroDataOrphan(:,1),EntroDataBounded(:,1),'tail','right')

%find reversible reactions with fixed direction
% 0 -> fixed reverse direction
% Inf -> fixed forward direction
% NaN -> inactive reaction
% r -> forward/reverse ratio
RevRxns=find(Recon1red.rev==1);
for i=1:length(RevRxns),
    FixedRev(i,1)=length(find(Liver00AO.Vsample(RevRxns(i),:)>eps))/length(find(Liver00AO.Vsample(RevRxns(i),:)<-eps));
    FixedRev(i,2)=length(find(Liver01AO.Vsample(RevRxns(i),:)>eps))/length(find(Liver01AO.Vsample(RevRxns(i),:)<-eps));
    FixedRev(i,3)=length(find(Liver03AO.Vsample(RevRxns(i),:)>eps))/length(find(Liver03AO.Vsample(RevRxns(i),:)<-eps));
    FixedRev(i,4)=length(find(Liver05AO.Vsample(RevRxns(i),:)>eps))/length(find(Liver05AO.Vsample(RevRxns(i),:)<-eps));
end
for i=1:4,
ProportionFixedRev(i)=length(find(FixedRev(:,i)==0 | isinf(FixedRev(:,i)))) / length(find(~isnan(FixedRev(:,i))));
end

%calculate Recon1 reference entropies without data integration
Recon1redsample=coneSampling(Recon1red,'F',Options00.FVArange,[],2000,[],'gurobi')
Recon1redEntro=DistEntropy(Recon1redsample.Vsample,[],'F','F');

%compare entropies
colnames={'Entropy(RegrEx)','Entropy(flux cone)','Rxn Idx','Data status','Rxn Name','Subsystem','Mechanism'};
[sortentro,sortidx]=sort(LiverEntro(:,1),'descend');
sortidx((sortentro==0))=[];sortentro(sortentro==0)=[];
datastatus=zeros(length(Recon1red.rxns),1);datastatus(DataRxns)=1;
s1Table=[colnames;[num2cell([sortentro,Recon1redEntro(sortidx),sortidx, datastatus(sortidx)]),Recon1red.rxnNames(sortidx),Recon1red.subSystems(sortidx),Recon1red.rxnMechanisms(sortidx)]];
% xlswrite('S1Table.xlsx',s1Table)

%get RegrEx alternative optima table
colnames={'L0','L0.1','L0.3','L0.5'};
rownames={'RegrExAO';'H_D';'H_nD';'Htotal';'median H_D';'median H_nD';'median Htotal';'F_rev';'D_mean';'D_CV'};
LiverRegrExAOtable=[rownames,[colnames;num2cell([sum(EntroDataBounded);sum(EntroDataOrphan);sum(EntroDataBounded)+sum(EntroDataOrphan);mean(EntroDataBounded);mean(EntroDataOrphan);mean([EntroDataBounded;EntroDataOrphan]);ProportionFixedRev;meanHDis;CVHDis],3)]];

%get figures
figure(7)
X=sort(LiverEntro(DataOrphan,:),'descend');
xcutoff=find(X<1e-3,1);
stairs(X(1:xcutoff,:),'LineWidth',1.5)
ylabel('H_s')
xlabel('rxns')
legend({'\lambda = 0','\lambda = 0.1','\lambda = 0.3','\lambda = 0.5'})
title('Data-orphan reactions')
figure(8)
X=sort(LiverEntro(DataRxns,:),'descend');
xcutoff=find(X<1e-3,1);
stairs(X(1:xcutoff,:),'LineWidth',1.5)
ylabel('H_s')
xlabel('rxns')
legend({'\lambda = 0','\lambda = 0.1','\lambda = 0.2','\lambda = 0.5'})
title('Data-bounded reactions')
figure(9)
boxplot(EntroDataOrphan)
set(gca,'XTickLabel',{'0','0.1','0.3','0.5'},'fontsize',16)
title('Data-orphan reactions')
figure(10)
boxplot(EntroDataBounded)
set(gca,'XTickLabel',{'0','0.1','0.3','0.5'},'fontsize',16)
title('Data-bounded reactions')

%rank-sum tests
[pDataBounded(1,1),~] = ranksum(EntroDataBounded(:,1),EntroDataBounded(:,2),'tail','right');
[pDataBounded(2,1),~] = ranksum(EntroDataBounded(:,2),EntroDataBounded(:,3),'tail','right');
[pDataBounded(3,1),~] = ranksum(EntroDataBounded(:,3),EntroDataBounded(:,4),'tail','right');
[pDataOrphan(1,1),~] = ranksum(EntroDataOrphan(:,1),EntroDataOrphan(:,2),'tail','right');
[pDataOrphan(2,1),~] = ranksum(EntroDataOrphan(:,2),EntroDataOrphan(:,3),'tail','right');
[pDataOrphan(3,1),~] = ranksum(EntroDataOrphan(:,3),EntroDataOrphan(:,4),'tail','right');
