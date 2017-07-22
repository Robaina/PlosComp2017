function m=patternfind(cellarray,pattern,CaseIns,first)

%Find patterns in cell array
%cellarray is a cell array...
%pattern is a character string
%CaseIns, either 0 if the search should be case sensitive or 1 if case
%insensitive
%first, either 0 if all instances of the pattern are required or 1 if only
%the first instance in the cellarray is required

%Semidán, September 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n=1;k=1;
    m=[];
    if nargin<3,
        CaseIns='F';
    end
    if isempty(CaseIns),
        CaseIns='F';
    end
    if nargin<4,
        first='F';
    end
    if isempty(first),
        first='F';
    end
    if ~iscellstr(cellarray(1)),
        for i=1:length(cellarray),
            cellarray(i)=cellarray{i};
        end
    end
    if CaseIns=='T',
        cellarray=lower(cellarray);
        pattern=lower(pattern);
    end

    l=strfind(cellarray,pattern);
    for i=1:size(l,1),
        for j=1:size(l,2),
            if l{i,j}~=0,
               m(n)=i;
               h(k)=j;
               n=n+1;
               k=k+1;
            end
        end
    end
    
    if first=='T',
        m=m(1);
    end

end
    