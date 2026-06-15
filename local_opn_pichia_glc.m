function local_opn_pichia_glc(mu, productionRatio, misfoldRatioOverride, options, candidateId)
%LOCAL_OPN_PICHIA_GLC Generate a pcSecPichia LP for secreted human OPN.
%
% OPN candidates are described in
% Data/pcSecPichia/TargetProtein_OPN_candidates.csv. Each candidate uses the
% mature human osteopontin sequence supplied by the user and swaps the
% N-terminal secretory leader.
%
% Usage:
%   local_opn_pichia_glc
%   local_opn_pichia_glc(0.10, 1e-6)
%   local_opn_pichia_glc(0.10, 1e-8, 0)
%   local_opn_pichia_glc(0.10, 1e-8, [], struct('mediaType', 4, 'blockMisfoldDilution', false))
%   local_opn_pichia_glc(0.10, 1e-8, [], struct('mediaType', 4), 'OPN_PPA_DDDK18')

if nargin < 1 || isempty(mu)
    mu = 0.10;
end
if nargin < 2 || isempty(productionRatio)
    productionRatio = 1e-6;
end
if nargin < 3
    misfoldRatioOverride = [];
end
if nargin < 4 || isempty(options)
    options = struct();
end
if nargin < 5 || isempty(candidateId)
    candidateId = getOption(options, 'candidateId', 'OPN_ALPHA_FULL_PROJECT');
end
candidateId = char(candidateId);

repoRoot = fileparts(mfilename('fullpath'));
if isempty(repoRoot)
    repoRoot = pwd;
end

startup_pcsec_local('pichia');

outDir = fullfile(repoRoot, 'local_runs', 'OPN_PPA_glc_smoke');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

targetCsv = fullfile(repoRoot, 'Data', 'pcSecPichia', 'TargetProtein_OPN_candidates.csv');
if ~exist(targetCsv, 'file')
    targetCsv = fullfile(repoRoot, 'Data', 'pcSecPichia', 'TargetProtein_OPN.csv');
end
fakeProteinInfo = readOpnTargetInfo(targetCsv, candidateId);
targetId = char(fakeProteinInfo{2});

load(fullfile(repoRoot, 'Enzymedata', 'pcSecPichia', 'enzymedataMachine_PP.mat'));
load(fullfile(repoRoot, 'Enzymedata', 'pcSecPichia', 'enzymedataSEC_PP.mat'));
load(fullfile(repoRoot, 'Enzymedata', 'pcSecPichia', 'enzymedataDummyER_PP.mat'));
load(resolveFirstExisting({ ...
    fullfile(repoRoot, 'Enzymedata', 'pcSecPichia', 'enzymedata_PP.mat'), ...
    fullfile(repoRoot, 'Enzymedata', 'pcSecPichia', 'enzymedat_PP.mat') ...
}));
load(fullfile(repoRoot, 'Model', 'pcSecPichia.mat'));

if ~isempty(misfoldRatioOverride)
    enzymedata.kdeg(:) = misfoldRatioOverride;
end

mediaType = getOption(options, 'mediaType', 2);
model = setMediaPP(model, mediaType);
model = changeRxnBounds(model, 'Ex_glc_D', -1000, 'l');
model = changeRxnBounds(model, 'BIOMASS', 1000, 'u');
model = changeRxnBounds(model, 'Ex_glyc', 0, 'l');
model = changeRxnBounds(model, 'BIOMASS_glyc', 0, 'b');
model = changeRxnBounds(model, 'Ex_meoh', 0, 'l');
model = changeRxnBounds(model, 'BIOMASS_meoh', 0, 'b');
model = changeRxnBounds(model, 'Ex_o2', -1000, 'l');
if getOption(options, 'blockMisfoldDilution', true)
    model.ub(contains(model.rxns, 'dilution_misfolding')) = 0;
end

tot_protein = 0.37;
f_modeled_protein = extractModeledprotein(model, 'BIOMASS', 'PROTEIN[c]');
f = tot_protein * f_modeled_protein;
f_unmodelER = 0.040;
factor_k = 1;

[model_tmp, enzymedataTP] = addTargetProtein(model, {targetId}, false, fakeProteinInfo);
enzymedataTP = SimulateRxnKcatCoef(model_tmp, enzymedataSEC, enzymedataTP);

enzymedata_new = enzymedata;
enzymedata_new.proteins = [enzymedata_new.proteins; enzymedataTP.proteins];
enzymedata_new.proteinMWs = [enzymedata_new.proteinMWs; enzymedataTP.proteinMWs];
enzymedata_new.proteinLength = [enzymedata_new.proteinLength; enzymedataTP.proteinLength];
enzymedata_new.proteinExtraMW = [enzymedata_new.proteinExtraMW; enzymedataTP.proteinExtraMW];
enzymedata_new.kdeg = [enzymedata_new.kdeg; enzymedataTP.kdeg];
enzymedata_new.proteinPST = [enzymedata_new.proteinPST; enzymedataTP.proteinPST];
enzymedata_new.rxns = [enzymedata_new.rxns; enzymedataTP.rxns];
enzymedata_new.rxnscoef = [enzymedata_new.rxnscoef; enzymedataTP.rxnscoef];
enzymedata_new = CombineEnzymedata(enzymedata_new, enzymedataSEC, enzymedataMachine, enzymedataDummyER);

model_tmp = changeRxnBounds(model_tmp, 'BIOMASS', mu, 'b');
model_tmp = changeRxnBounds(model_tmp, [targetId, ' exchange'], productionRatio, 'b');

cd(outDir);
rxnID = 'Ex_glc_D';
osenseStr = 'Maximize';
modeTag = optionTag(options);
candidateTag = safeText(targetId);
name = sprintf('%s_mu%s_media%d_ratio%s_misfold%s%s_PP', candidateTag, safeNum(mu), mediaType, safeNum(productionRatio), safeOptionalNum(misfoldRatioOverride), modeTag);
fileName = writeLPGlc(model_tmp, mu, f, f_unmodelER, osenseStr, rxnID, enzymedata_new, factor_k, name, 0, options);
writeclusterfileLP({fileName}, sprintf('sub_%s_PPA_glc_smoke_media%d_misfold%s%s', candidateTag, mediaType, safeOptionalNum(misfoldRatioOverride), modeTag));
save(sprintf('%s_pichia_glc_smoke_media%d_misfold%s%s_setup.mat', candidateTag, mediaType, safeOptionalNum(misfoldRatioOverride), modeTag), ...
    'mu', 'mediaType', 'productionRatio', 'misfoldRatioOverride', 'options', 'candidateId', 'targetCsv', 'fakeProteinInfo', 'enzymedataTP', 'fileName');
cd(repoRoot);

fprintf('OPN pcSecPichia LP generated for %s: %s\n', targetId, fullfile(outDir, fileName));

end

function fakeProteinInfo = readOpnTargetInfo(targetCsv, candidateId)
raw = readcell(targetCsv, 'Delimiter', ',');
if size(raw, 1) < 2
    error('No OPN target row found in %s.', targetCsv);
end
data = raw(2:end, 1:15);
candidateId = char(candidateId);
rowIdx = find(strcmp(data(:, 1), candidateId) | strcmp(data(:, 2), candidateId), 1);
if isempty(rowIdx)
    if size(data, 1) == 1
        rowIdx = 1;
    else
        available = strjoin(string(data(:, 1)), ', ');
        error('Candidate %s was not found in %s. Available candidates: %s', candidateId, targetCsv, available);
    end
end
fakeProteinInfo = data(rowIdx, 1:15);
numericCols = [3 4 5 6 7 8 9 12 14 15];
for idx = numericCols
    value = fakeProteinInfo{idx};
    if ischar(value) || isstring(value)
        fakeProteinInfo{idx} = str2double(value);
    end
end
end

function pathName = resolveFirstExisting(candidates)
for i = 1:numel(candidates)
    if exist(candidates{i}, 'file')
        pathName = candidates{i};
        return;
    end
end
error('None of the candidate files exist.');
end

function text = safeNum(value)
text = strrep(sprintf('%.2g', value), '.', 'p');
text = strrep(text, '-', 'm');
end

function text = safeOptionalNum(value)
if isempty(value)
    text = 'default';
else
    text = safeNum(value);
end
end

function text = safeText(value)
text = regexprep(char(value), '[^A-Za-z0-9]+', '_');
text = regexprep(text, '^_+|_+$', '');
end

function value = getOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName)
    value = options.(fieldName);
else
    value = defaultValue;
end
end

function tag = optionTag(options)
parts = {};
if ~getOption(options, 'blockMisfoldDilution', true)
    parts{end+1} = 'openMisfoldDilution';
end
if ~getOption(options, 'writeMisfoldingConstraints', true)
    parts{end+1} = 'noMisfoldEq';
end
if ~getOption(options, 'writeRibosomeConstraint', true)
    parts{end+1} = 'noRiboEq';
end
if isempty(parts)
    tag = '';
else
    tag = ['_', strjoin(parts, '_')];
end
end
