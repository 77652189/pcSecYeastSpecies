%% buildModel
% Timing: ~  s

% Before run this script, please 1) check the annotation file is correct
% for all sce proteins
addpath('../../Enzymedata/')
addpath('../../Model/')
addpath('../../pcSecPichia/');

[~,~,protein_info] = xlsread('Protein_information_PP.xlsx');
tic;

load('iMT1026modif.mat');
%% Modify the original model
model_updated = modifiMTv3(model);
model = model_updated;

%% Reformulate the orginal model
%model2 = ravenCobraWrapper(model);
%model2.metChEBIID(1:length(model2.mets),1) = {''};
%model2.metKEGGID(1:length(model2.mets),1) = {''};
%model_updated = model2;
%Split Reactions
model_split = splitModel_iMT(model_updated);
model=model_split;%

load('Protein_Sequence.mat');
load('Protein_stoichiometry.mat');
model = addComplexRxns(model,Protein_stoichiometry,protein_info,ProteinSequence);

%% reformulate the Sec part
[~,~,proteins] = xlsread('TableS1_PP.xlsx','Secretory');
model = addMachineryComplex(model,proteins,protein_info,Protein_stoichiometry,ProteinSequence);

%% reformulate the ribosome/assembly factor/proteasome part
[~,~,protein_machinery]=xlsread('TableS1_PP.xlsx','Machinery');
model = addMachineryComplex(model,protein_machinery,protein_info,Protein_stoichiometry,ProteinSequence);

% manually update complex formation reactions for some complexes.
% model = updateComplexformation(model,protein_info);
%% Collect kcats for enzymes

% use the one for with sec parameter
enzymedata = collectkcats(model);
% load('Kcatenzymedata.mat');
list = endsWith(model.rxns,'_translation');
list = model.rxns(list);
list = strrep(list,'r_','');
enzymedata.proteins = strrep(list,'_peptide_translation',''); % all protein

% get kdeg info from petri's proteome data
enzymedata = getkdeg(enzymedata);

% Calculate molecular weight for each enzyme
enzymedata = calculateMW(enzymedata,ProteinSequence,protein_info);

% which is based on Sec protein abundance and the sec machinery abundance
[enzymedataSEC,modeled_ratio,meanprotein_info,missingsecP_ratio] = SimulateSecParam(model,protein_info,ProteinSequence);
enzymedataSEC = calculateMW(enzymedataSEC,ProteinSequence,protein_info);
enzymedata = SimulateRxnKcatCoef(model,enzymedataSEC,enzymedata);

% calculate the enzyme molecule weight of other machinery proteins such as
% ribosome and proteasome
enzymedataMachine = SimulateMachineParam(model);
enzymedataMachine = calculateMW(enzymedataMachine,ProteinSequence,protein_info);

%% Add dummy complex reactions
% Dummy complex is assumed to be a part of metabolic protein pool.
% Note that the dummy complex is synthesized or diluted in the unit of
% mmol/gCDW/h.
[protein_info,~] = SimulateDummyERParam(model,meanprotein_info,protein_info,enzymedataSEC);
model = addDummyER(model,'PROTEINS','PROTEIN[c]',protein_info); 
[~,enzymedataDummyER] = SimulateDummyERParam(model,meanprotein_info,protein_info,enzymedataSEC);

model = addDummy(model,'PROTEINS','PROTEIN[c]'); 
% close refolding and dilute for misfolding protein those reactions are for
% later usage
refold_list = model.rxns(contains(model.rxns,'refolding_'));
model = changeRxnBounds(model,refold_list,0,'b');
% misfold_dilute_list = model.rxns(contains(model.rxns,'dilution_misfolding'));
% model = changeRxnBounds(model,misfold_dilute_list,0,'b');

%% Change the original biomass equation
% Estimate modeled proteome
f_modeled_protein = estimateModeledprotein(model);

f_modeled_protein = f_modeled_protein + missingsecP_ratio; % small amount of protein has been added into ER
% Change the biomass equation
model = changeBiomass(model,f_modeled_protein,'BIOMASS','PROTEIN[c]');%%

% change GAM and NGAM
model = fitGAMforGlc(model);

%% Save model
save('pcSecPichia.mat','model');
save('enzymedata_PP.mat','enzymedata');
save('enzymedataSEC_PP.mat','enzymedataSEC');
save('enzymedataMachine_PP.mat','enzymedataMachine');
save('enzymedataDummyER_PP.mat','enzymedataDummyER');

%% Save model to Excel
model_excel = model;
model_excel.subSystems = cell(length(model_excel.rxns),1);
model_excel.mets = model_excel.metNames;
writeCbModel(model_excel,'xls','test\pcSecPichia.xls');
clear model_excel;
toc;
