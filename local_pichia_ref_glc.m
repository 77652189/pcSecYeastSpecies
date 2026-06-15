function local_pichia_ref_glc(mu, mediaType, misfoldRatioOverride, options)
%LOCAL_PICHIA_REF_GLC Generate one no-target pcSecPichia glucose LP.
%
% Usage:
%   local_pichia_ref_glc
%   local_pichia_ref_glc(0.10, 2)
%   local_pichia_ref_glc(0.10, 4)
%   local_pichia_ref_glc(0.10, 4, 0)
%   local_pichia_ref_glc(0.10, 2, [], struct('blockMisfoldDilution', false))

if nargin < 1 || isempty(mu)
    mu = 0.10;
end
if nargin < 2 || isempty(mediaType)
    mediaType = 2;
end
if nargin < 3
    misfoldRatioOverride = [];
end
if nargin < 4 || isempty(options)
    options = struct();
end

repoRoot = fileparts(mfilename('fullpath'));
if isempty(repoRoot)
    repoRoot = pwd;
end

startup_pcsec_local('pichia');

runDir = fullfile(repoRoot, 'local_runs', 'PPA_GLC_ref_smoke');
if ~exist(runDir, 'dir')
    mkdir(runDir);
end

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

model = setMediaPP(model, mediaType);
model = changeRxnBounds(model, 'Ex_glc_D', -1000, 'l');
model = changeRxnBounds(model, 'BIOMASS', 1000, 'u');
model = changeRxnBoundsIfPresent(model, 'LIPIDS', 1000, 'u');
model = changeRxnBoundsIfPresent(model, 'PROTEINS', 1000, 'u');
model = changeRxnBoundsIfPresent(model, 'STEROLS', 1000, 'u');
model = changeRxnBounds(model, 'Ex_glyc', 0, 'l');
model = changeRxnBoundsIfPresent(model, 'BIOMASS_glyc', 0, 'b');
model = changeRxnBoundsIfPresent(model, 'LIPIDS_glyc', 0, 'b');
model = changeRxnBoundsIfPresent(model, 'PROTEINS_glyc', 0, 'b');
model = changeRxnBoundsIfPresent(model, 'STEROLS_glyc', 0, 'b');
model = changeRxnBounds(model, 'Ex_meoh', 0, 'l');
model = changeRxnBoundsIfPresent(model, 'BIOMASS_meoh', 0, 'b');
model = changeRxnBoundsIfPresent(model, 'LIPIDS_meoh', 0, 'b');
model = changeRxnBoundsIfPresent(model, 'PROTEINS_meoh', 0, 'b');
model = changeRxnBoundsIfPresent(model, 'STEROLS_meoh', 0, 'b');
model = changeRxnBounds(model, 'Ex_o2', -1000, 'l');
if getOption(options, 'blockMisfoldDilution', true)
    model.ub(contains(model.rxns, 'dilution_misfolding')) = 0;
end

tot_protein = 0.37;
f_modeled_protein = extractModeledprotein(model, 'BIOMASS', 'PROTEIN[c]');
f = tot_protein * f_modeled_protein;
f_unmodelER = tot_protein * 0.04;
factor_k = 1;
rxnID = 'Ex_glc_D';
osenseStr = 'Maximize';
enzymedata_all = CombineEnzymedata(enzymedata, enzymedataSEC, enzymedataMachine, enzymedataDummyER);

model_tmp = changeRxnBounds(model, 'BIOMASS', mu, 'b');

oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(runDir);

modeTag = optionTag(options);
name = sprintf('ref_mu%s_media%d_misfold%s%s_PP', safeNum(mu), mediaType, safeOptionalNum(misfoldRatioOverride), modeTag);
fileName = writeLPGlc(model_tmp, mu, f, f_unmodelER, osenseStr, rxnID, enzymedata_all, factor_k, name, 0, options);
writeclusterfileLP({fileName}, sprintf('sub_ref_media%d_misfold%s%s', mediaType, safeOptionalNum(misfoldRatioOverride), modeTag));
save(sprintf('ref_mu%s_media%d_misfold%s%s_setup.mat', safeNum(mu), mediaType, safeOptionalNum(misfoldRatioOverride), modeTag), ...
    'mu', 'mediaType', 'misfoldRatioOverride', 'options', 'fileName');
fprintf('pcSecPichia reference LP generated: %s\n', fullfile(runDir, fileName));

clear cleanupObj;

end

function model = changeRxnBoundsIfPresent(model, rxnID, value, boundType)
if any(strcmp(model.rxns, rxnID))
    model = changeRxnBounds(model, rxnID, value, boundType);
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
