function Sol=RegrExOLSLADstudy(GEM,Data,samplesize)

%Create sample of perturbed data vectors
Rxns=size(GEM.S,2);
% Data=Data/max(Data);
DataRxns=find(Data>0);
DataSample=zeros(Rxns,samplesize);
epsmax=0.01*mean(Data(DataRxns));epsmin=-epsmax;
for i=1:samplesize,
    DataSample(DataRxns,i)=Data(DataRxns)+(epsmax-epsmin).*rand(length(DataRxns),1)+epsmin;
    DataSample(:,i)=DataSample(:,i)/max(DataSample(:,i));
end

%Introduce outlier(s) in original data
[~,sortedDrxn]=sort(Data(DataRxns),'ascend');
Options.L_int=0;Options.L_min=0;Options.L_max=0;Options.time_limit=60;

DataOut=Data;DataOut(DataRxns(sortedDrxn(1)))=5*max(Data(DataRxns));
DataOut=DataOut/max(DataOut);
Data=Data/max(Data);
%Obtain context-specific flux distributions
%with outliers
VolsOut=RegrExOLS(GEM,DataOut,Options);
VladOut=RegrExLAD(GEM,DataOut,Options);
%without outliers
Vols=RegrExOLS(GEM,Data,Options);
Vlad=RegrExLAD(GEM,Data,Options);

%Obtain mean squared error over the sample data
for i=1:samplesize,
   OLSerrorOut(i,1)=sum(abs(DataSample(DataRxns,i)-VolsOut.Flux(DataRxns)));
   LADerrorOut(i,1)=sum(abs(DataSample(DataRxns,i)-VladOut.Flux(DataRxns)));
   OLSerror(i,1)=sum(abs(DataSample(DataRxns,i)-Vols.Flux(DataRxns)));
   LADerror(i,1)=sum(abs(DataSample(DataRxns,i)-Vlad.Flux(DataRxns)));
end
[pvalueOut,~]=ranksum(LADerrorOut(:,1),OLSerrorOut(:,1),'tail','left');
[pvalue,~]=ranksum(LADerror(:,1),OLSerror(:,1),'tail','left');
OLSmsreOut=mean(OLSerrorOut(:,1));
LADmsreOut=mean(LADerrorOut(:,1));
OLSmsre=mean(OLSerror(:,1));
LADmsre=mean(LADerror(:,1));

Sol.pvalueOut=pvalueOut;
Sol.pvalue=pvalue;
Sol.OLSLADmsreOut=[OLSmsreOut/length(DataRxns),LADmsreOut/length(DataRxns)];
Sol.OLSLADmsre=[OLSmsre/length(DataRxns),LADmsre/length(DataRxns)];

end