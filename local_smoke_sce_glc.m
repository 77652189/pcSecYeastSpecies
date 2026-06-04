function local_smoke_sce_glc(mu)
%LOCAL_SMOKE_SCE_GLC Generate one S. cerevisiae glucose LP for SoPlex.

if nargin < 1 || isempty(mu)
    mu = 0.01;
end

repoRoot = fileparts(mfilename('fullpath'));
if isempty(repoRoot)
    repoRoot = pwd;
end

startup_pcsec_local('yeast');

runDir = fullfile(repoRoot, 'local_runs', 'SCE_GLC_smoke');
if ~exist(runDir, 'dir')
    mkdir(runDir);
end

modelData = load(fullfile(repoRoot, 'Model', 'pcSecYeast.mat'), 'model');
enzData = load(fullfile(repoRoot, 'Enzymedata', 'pcSecYeast', 'enzymedata_SCE.mat'), 'enzymedata');
secData = load(fullfile(repoRoot, 'Enzymedata', 'pcSecYeast', 'enzymedataSEC_SCE.mat'), 'enzymedataSEC');
dummyERData = load(fullfile(repoRoot, 'Enzymedata', 'pcSecYeast', 'enzymedataDummyER_SCE.mat'), 'enzymedataDummyER');
machineData = load(fullfile(repoRoot, 'Enzymedata', 'pcSecYeast', 'enzymedataMachine_SCE.mat'), 'enzymedataMachine');

model = modelData.model;
enzymedata = enzData.enzymedata;
enzymedataSEC = secData.enzymedataSEC;
enzymedataDummyER = dummyERData.enzymedataDummyER;
enzymedataMachine = machineData.enzymedataMachine;

model = setMedia(model, 2);
model = changeRxnBounds(model, 'r_1714', -1000, 'l');
model = changeRxnBounds(model, 'r_1709', 0, 'l');
model = changeRxnBounds(model, 'r_2058', 0, 'l');
model = changeRxnBounds(model, 'r_1931', 0, 'l');
model = changeRxnBounds(model, 'r_1710', 0, 'l');
model = changeRxnBounds(model, 'r_1992', -1000, 'l');
model = blockRxns(model);
model = changeRxnBounds(model, 'r_1634', 0, 'b');
model = changeRxnBounds(model, 'r_1631', 0, 'b');
model = changeRxnBounds(model, 'r_1810', 0, 'b');
model = changeRxnBounds(model, 'r_2033', 0, 'b');

rxnID = 'r_1714';
osenseStr = 'Maximize';
tot_protein = 0.46;
f_modeled_protein = extractModeledprotein(model, 'r_4041', 's_3717[c]');
f = tot_protein * f_modeled_protein;
f_unmodelER = tot_protein * 0.046;
factor_k = 1;

enzymedata_all = CombineEnzymedata(enzymedata, enzymedataSEC, enzymedataMachine, enzymedataDummyER);
model.ub(contains(model.rxns, 'dilution_misfolding')) = 0;
model_tmp = changeRxnBounds(model, 'r_2111', mu, 'b');

oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(runDir);

muTag = strrep(sprintf('%.2f', mu), '.', 'p');
fileName = writeLPSCE(model_tmp, mu, f, f_unmodelER, osenseStr, rxnID, enzymedata_all, factor_k, [muTag, '_GLC_SCE_smoke']);
writeclusterfileLP({fileName}, 'sub_1');

fprintf('SCE glucose smoke LP generated: %s\n', fullfile(runDir, fileName));
fprintf('SoPlex runner generated: %s\n', fullfile(runDir, 'sub_1.sh'));
fprintf('Run in WSL from this directory: bash sub_1.sh\n');

clear cleanupObj;

end
