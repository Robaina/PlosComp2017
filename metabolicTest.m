function Sol=metabolicTest(GEM,A)

%This function implements the metabolic test employed in:

%Schultz A, Qutub AA. Reconstruction of Tissue-Specific Metabolic Networks
%Using CORDA. PLoS Comput Biol [Internet]. Public Library of Science; 2016 %

%Get context-specif model COBRA structure
if nargin<2,
    A=[];
end
if ~isempty(A),
    Abinary=zeros(length(GEM.rxns),1);Abinary(A)=1;
    GEM=contextmodel2COBRA(Abinary,GEM,'',0,0);
end
BasalEx={'EX_co2(e)';'EX_h2o(e)';'EX_h(e)';'EX_o2(e)';'EX_pi(e)';'EX_h2o2(e)';'EX_o2s(e)';'EX_hco3(e)';'EX_co(e)'};
AminoacidTest={'EX_glu_L(e)';'EX_gly(e)';'EX_ala_L(e)';'EX_lys_L(e)';'EX_asp_L(e)';'EX_arg_L(e)';'EX_gln_L(e)';'EX_ser_L(e)';'EX_met_L(e)';'EX_trp_L(e)';'EX_phe_L(e)';'EX_tyr_L(e)';'EX_cys_L(e)';'EX_leu_L(e)';'EX_his_L(e)';'EX_pro_L(e)';'EX_asn_L(e)';'EX_val_L(e)';'EX_ile_L(e)'};
GlucoseTest={'EX_glu_L(e)';'EX_gly(e)';'EX_ala_L(e)';'EX_asp_L(e)';'EX_arg_L(e)';'EX_gln_L(e)';'EX_ser_L(e)';'EX_met_L(e)';'EX_trp_L(e)';'EX_phe_L(e)';'EX_tyr_L(e)';'EX_cys_L(e)';'EX_his_L(e)';'EX_pro_L(e)';'EX_asn_L(e)';'EX_val_L(e)';'EX_thr_L(e)';'EX_ile_L(e)';'EX_lac_L(e)';'EX_glyc(e)';'EX_pyr(e)'};
NucleotideTest = {'DM_atp[c]' 'atp[c] -> ';'DM_ctp[c]' 'ctp[c] -> ';'DM_datp[c]' 'datp[c] -> ';'DM_dctp[c]' 'dctp[c] -> ';'DM_dgtp[c]' 'dgtp[c] -> ';'DM_dttp[c]' 'dttp[c] -> ';'DM_gtp[c]' 'gtp[c] -> ';'DM_utp[c]' 'utp[c] -> '};

GEM.rxnMechanisms=createRxnMechanism(GEM);
ExchangeRxns=patternfind(GEM.rxnMechanisms,'Extracellular');

%Aminoacid Test
mod=GEM;
TestAminoResult=zeros(length(AminoacidTest),2);
%Block intake of exchange reactions (minus basal reactions).Block excretion
%of ammonia
mod.lb(ExchangeRxns)=0;
mod.ub(strcmp(mod.rxns,'EX_nh4(e)'))=0;
%Allow exchange of basal reactions and glucose intake
mod.lb(ismember(mod.rxns,BasalEx)==1)=-1000;
mod.lb(strcmp(mod.rxns,'EX_glc(e)'))=-1000;
%Test aminoacid production
mod.c=zeros(length(mod.rxns),1);
%Set fobj to urea exchange
mod.c(strcmp(mod.rxns,'EX_urea(e)'))=1;

for i=1:length(AminoacidTest),
    mod2=mod;
    mod2.lb(strcmp(mod2.rxns,AminoacidTest{i}))=-1000;
    %maximize urea excretion
    Val=optimizeCbModel(mod2);
    TestAminoResult(i,1)=Val.f;
    TestAminoResult(i,2)=Val.stat;
end
TestAminoResult(TestAminoResult<1e-6)=0;
TestAminoResultT=[{'Metabolite','Fobj','Status'};[AminoacidTest,num2cell(TestAminoResult)]];

%Glucose Test
mod=GEM;
TestGlucoResult=zeros(length(GlucoseTest),2);
%Block intake of exchange reactions (minus basal reactions)
mod.lb(ExchangeRxns)=0;
mod.ub(strcmp(mod.rxns,'EX_nh4(e)'))=0;
%Allow exchange of basal reactions and glucose intake
mod.lb(ismember(mod.rxns,BasalEx)==1)=-1000;
%Test glucose production
mod.c=zeros(length(mod.rxns),1);
mod.c(strcmp(mod.rxns,'EX_glc(e)'))=1;

for i=1:length(GlucoseTest),
    mod2=mod;
    mod2.lb(strcmp(mod2.rxns,GlucoseTest{i}))=-1000;
    %Maximize glucose excretion
    Val=optimizeCbModel(mod2);
    TestGlucoResult(i,1)=Val.f;
    TestGlucoResult(i,2)=Val.stat;
end
TestGlucoResult(TestGlucoResult<1e-6)=0;
TestGlucoResultT=[{'Metabolite','Fobj','Status'};[GlucoseTest,num2cell(TestGlucoResult)]];

%Nucleotide Test
mod=GEM;
TestNucleoResult=zeros(length(NucleotideTest),2);
%Block intake of exchange reactions (minus basal reactions)
mod.lb(ExchangeRxns)=0;
mod.ub(strcmp(mod.rxns,'EX_nh4(e)'))=0;
%Allow exchange of basal reactions and glucose intake
mod.lb(ismember(mod.rxns,BasalEx)==1)=-1000;
mod.lb(strcmp(mod.rxns,'EX_glc(e)'))=-1000;
mod.lb(strcmp(mod.rxns,'EX_nh4(e)'))=-1000;
mod.c=zeros(length(mod.rxns),1);

for i=1:length(NucleotideTest(:,1)),
    mod2=mod;
    %Add sink reaction
    mod2=addReaction(mod2,NucleotideTest{i,1},NucleotideTest{i,2});
    mod2.c(strcmp(mod2.rxns,NucleotideTest{i,1}))=1;
    %Maximize nucleotide sink reaction
    Val=optimizeCbModel(mod2);
    TestNucleoResult(i,1)=Val.f;
    TestNucleoResult(i,2)=Val.stat;
end
TestNucleoResult(TestNucleoResult<1e-6)=0;
TestNucleoResultT=[{'Nucleotide Sink','Fobj','Status'};[NucleotideTest(:,2),num2cell(TestNucleoResult)]];

Sol.TestAmino=TestAminoResultT;
Sol.TestGluco=TestGlucoResultT;
Sol.TestNucleo=TestNucleoResultT;
Sol.TestScore=(length(find(TestAminoResult(:,1)>0))+length(find(TestGlucoResult(:,1)>0))+length(find(TestNucleoResult(:,1)>0)))/(length(AminoacidTest)+length(GlucoseTest)+length(NucleotideTest(:,1)));
end
    
    