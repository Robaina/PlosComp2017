function Entropy=DistEntropy(X,Nbins,scaled,relative)

%Calculates the (relative) entropy of a distribution
%Semidan, April 2015
%**************************************************************************

if nargin<2,
   Nbins=20;
end
if isempty(Nbins),
   Nbins=20;
end

if nargin<3,
   scaled='F';
end
if isempty(scaled),
   scaled='F';
end  

if nargin<4,
   relative='F';
end
if isempty(relative),
   relative='F';
end  

E=zeros(Nbins,1); 
p=size(X,2)/(Nbins*size(X,2));
Entropy=zeros(size(X,1),1);
for i=1:size(X,1),
    [counts,~]=hist(X(i,:),Nbins);
    counts=counts/sum(counts);
    for j=1:Nbins, 
        E(j)=counts(j)*log(counts(j));
        if isnan(E(j)),
            E(j)=0;
        end
    end
    Entropy(i)=-sum(E);
    if scaled=='T',
       Entropy(i)=abs(max(X(i,:))-min(X(i,:)))*Entropy(i);
    end
    if relative=='T',
        Entropy(i)=Entropy(i)/(-log(p));
    end

end