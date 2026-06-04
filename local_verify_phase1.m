function local_verify_phase1()
%LOCAL_VERIFY_PHASE1 Verify model loading and one processed-result figure.

repoRoot = fileparts(mfilename('fullpath'));
if isempty(repoRoot)
    repoRoot = pwd;
end

startup_pcsec_local('yeast');

runDir = fullfile(repoRoot, 'local_runs', 'phase1');
if ~exist(runDir, 'dir')
    mkdir(runDir);
end

modelFile = fullfile(repoRoot, 'Model', 'pcSecYeast.mat');
assert(exist(modelFile, 'file') == 2, 'Missing model file: %s', modelFile);

loadedModel = load(modelFile, 'model');
assert(isfield(loadedModel, 'model'), 'Model file does not contain variable "model".');
assert(isfield(loadedModel.model, 'rxns') && ~isempty(loadedModel.model.rxns), 'model.rxns is empty.');
assert(isfield(loadedModel.model, 'mets') && ~isempty(loadedModel.model.mets), 'model.mets is empty.');

figScript = fullfile(repoRoot, 'Code', 'Figures', 'Fig1b_ModelComparisonPpa.m');
assert(exist(figScript, 'file') == 2, 'Missing figure script: %s', figScript);

oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(fileparts(figScript));
run(figScript);

outputPng = fullfile(runDir, 'Fig1b_ModelComparisonPpa.png');
saveas(gcf, outputPng);
fprintf('Phase 1 verification passed. Figure saved to: %s\n', outputPng);

clear cleanupObj;

end
