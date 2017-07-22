function [maprxnvalues]=mapgene2rxn(genenames,genevalues,modelgenes,modelgrules,modellogic)


% c=0;
% %Modellogic=1 => A | B, A+B
% %Modellogic=0 => A | B, max(A,B)
%Maps gene value to gene name in model genes list
%genenames must be a cell of char entries     

 
modeldata=[];loc=[];geneloc=[];
for n=1:length(modelgenes),
    for p=1:length(genenames),
    loc(p)=strcmp(modelgenes{n},genenames{p});
    end
    try
        geneloc(n)=find(loc==1,1);
        modeldata(n)=genevalues(geneloc(n));
    catch
        if modellogic==1,
        modeldata(n)=0;
        elseif modellogic==0,
        modeldata(n)=nan;
        end
     
    end
        
end
for k=1:length(modelgenes),
    VAMOS{k}=num2str(modeldata(k));
end

%Asings value to gene name in reaction rules

VAMOS=VAMOS';
mrm=modelgrules;
for k=1:length(modelgenes),
    mrm=strrep(mrm,modelgenes{k},VAMOS{k});   
    mrm=strrep(mrm,'OR','|');
    mrm=strrep(mrm,'AND','&');
    %mrm=strrep(mrm,modelgenedata(n,1),modelgenedata(n,3)); 
end

for n=1:length(mrm)
    
    %Selects substring within last pair of brackets
    
while isempty(strfind(mrm{n},'('))==0,

    lastbracket=strfind(mrm{n},'(');
    firstbracket=strfind(mrm{n},')');
    lastsubstr=mrm{n}(max(lastbracket):firstbracket(find(firstbracket>max(lastbracket),1)));

if isempty(strfind(lastsubstr,'&'))==0,
    lastsubstr=strrep(lastsubstr,'&',',');
    lastsubstr=strrep(lastsubstr,'(','[');lastsubstr=strrep(lastsubstr,')',']');lastsubstr=strcat('min','(',lastsubstr,')');
    
elseif isempty(strfind(lastsubstr,'|'))==0,
    if modellogic==0,
        lastsubstr=strrep(lastsubstr,'|',',');
        lastsubstr=strrep(lastsubstr,'(','[');lastsubstr=strrep(lastsubstr,')',']');lastsubstr=strcat('max','(',lastsubstr,')');
    elseif modellogic==1,
        lastsubstr=strrep(lastsubstr,'|','+');
    end
end

     %Evaluates expression within last brackets & substitutes value in logical rule

    try
       lastsubeval{n}=eval(lastsubstr); 
       mrm{n}=strrep(mrm{n},mrm{n}(max(lastbracket):firstbracket(find(firstbracket>max(lastbracket),1))),sprintf('%d',lastsubeval{n}));
    catch
        lastsubeval{n}=0;
        mrm{n}=strrep(mrm{n},'(','');mrm{n}=strrep(mrm{n},')','');
    end
     
end

    % Perfomrs substitution when no brackets left 

if isempty(strfind(mrm{n},'('))==1 && (isempty(strfind(mrm{n},'&'))==0) && (isempty(strfind(mrm{n},'|'))==1),
   mrm{n}=strrep(mrm{n},'&',',');mrm{n}=strcat('min','(','[',mrm{n},']',')');   
    
elseif isempty(strfind(mrm{n},'('))==1 && (isempty(strfind(mrm{n},'|'))==0) && (isempty(strfind(mrm{n},'&'))==1),
    if modellogic==0,
        mrm{n}=strrep(mrm{n},'|',',');mrm{n}=strcat('max','(','[',mrm{n},']',')');
    elseif modellogic==1,
        mrm{n}=strrep(mrm{n},'|','+');
    end
elseif isempty(strfind(mrm{n},'('))==1 && (isempty(strfind(mrm{n},'&'))==0) && (isempty(strfind(mrm{n},'|'))==0),
    mrm{n}=strrep(mrm{n},'&',',');isor=strfind(mrm{n},'|');
    if length(isor)==1,
    mrm{n}=strcat('min','(','[',mrm{n}(1:isor(1)-1),']',')',mrm{n}(isor(1)),'min','(','[',mrm{n}(isor(1)+1:length(mrm{n})),']',')');
    elseif length(isor)>1,
        mrm{n}=strcat('min','(','[',mrm{n}(1:isor(1)-1),']',')',mrm{n}(isor(1):length(mrm{n})));
        for q=2:length(isor),
            if q==length(isor),
               isor=strfind(mrm{n},'|');
               mrm{n}=strcat(mrm{n}(1:isor(q-1)),'min','(','[',mrm{n}(isor(q-1)+1:isor(q)-1),']',')',mrm{n}(isor(q)),'min','(','[',mrm{n}(isor(q)+1:length(mrm{n})),']',')');
            else
              isor=strfind(mrm{n},'|');
              mrm{n}=strcat(mrm{n}(1:isor(q-1)),'min','(','[',mrm{n}(isor(q-1)+1:isor(q)-1),']',')',mrm{n}(isor(q):length(mrm{n})));
            end            
        end
        
    end

end
end
  
for n=1:length(mrm),
    if isempty(mrm{n})==0,
        if modellogic==0,
        mrm{n}=strrep(mrm{n},'|',',');mrm{n}=strcat('max','(','[',mrm{n},']',')');
        elseif modellogic==1,
        mrm{n}=strrep(mrm{n},'|','+');
        end
        try
        mrm{n}=strrep(mrm{n},mrm{n},sprintf('%d',eval(mrm{n})));
        catch
            mrm{n}=modelgrules{n};
        end
    end
end

    %Pinpoints mistakes in model.rules   

errormrm=[];
for t=1:length(mrm),
   errormrm(t)=isempty(strfind(mrm{t},'T'));
end
 numerrors=find(errormrm==0); 
if numerrors>0,
   disp('');
   fprintf('%g reaction rules have been found to be incomplete',length(numerrors));
   disp('');
   disp('Please, check following reaction rules: ');
   disp('');
   mrm{[numerrors]}
    
end

maprxnvalues=[];
for i=1:length(mrm),
    if isempty(mrm{i})==1,
        mrm{i}='0';
    end
   maprxnvalues(i)=str2num(mrm{i});
end
   maprxnvalues=maprxnvalues';
   maprxnvalues(isnan(maprxnvalues)==1)=0;
  

end
        
    
    