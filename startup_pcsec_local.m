function startup_pcsec_local(species)
%STARTUP_PCSEC_LOCAL Configure local MATLAB paths for this repository.
%
% Usage:
%   startup_pcsec_local
%   startup_pcsec_local('yeast')
%   startup_pcsec_local('pichia')
%   startup_pcsec_local('kmarx')

if nargin < 1 || isempty(species)
    species = 'all';
end

repoRoot = fileparts(mfilename('fullpath'));
if isempty(repoRoot)
    repoRoot = pwd;
end
cd(repoRoot);

addpath(repoRoot, '-begin');

toolboxRoot = fullfile(getenv('USERPROFILE'), 'MATLAB', 'toolboxes');
addExistingPath(fullfile(toolboxRoot, 'cobratoolbox'));
addExistingPath(fullfile(toolboxRoot, 'RAVEN'));

addExistingPath(fullfile(repoRoot, 'Model'));
addExistingPath(fullfile(repoRoot, 'Results'));
addExistingPath(fullfile(repoRoot, 'Code', 'Figures'));

switch lower(species)
    case {'sce', 'yeast', 'pcsecyeast'}
        addSpeciesPath(repoRoot, 'pcSecYeast');
    case {'ppa', 'pichia', 'pcsecpichia'}
        addSpeciesPath(repoRoot, 'pcSecPichia');
    case {'kmx', 'kmarx', 'pcseckmarx'}
        addSpeciesPath(repoRoot, 'pcSecKmarx');
    case 'all'
        addSpeciesPath(repoRoot, 'pcSecYeast');
        addSpeciesPath(repoRoot, 'pcSecPichia');
        addSpeciesPath(repoRoot, 'pcSecKmarx');
    otherwise
        error('Unknown species "%s". Use yeast, pichia, kmarx, or all.', species);
end

fprintf('pcSecYeastSpecies local paths configured: %s\n', repoRoot);

end

function addSpeciesPath(repoRoot, speciesDir)
addExistingPath(fullfile(repoRoot, 'Code', speciesDir));
addExistingPath(fullfile(repoRoot, 'Data', speciesDir));
addExistingPath(fullfile(repoRoot, 'Enzymedata', speciesDir));
end

function addExistingPath(pathName)
if exist(pathName, 'dir')
    addpath(genpath(pathName), '-begin');
end
end
