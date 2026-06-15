%% setMediaPP
function model = setMediaPP(model,type)
%SETMEDIAPP Configure exchange bounds for pcSecPichia media.
%
% type  = 1: minimal medium
%       = 2: yeast nitrogen base without amino acids
%       = 3: YNB + CSM-Ura-like amino acid supplement
%       = 4: minimal + amino acids + uracil
%       = 5: YNB with all amino acids

if nargin < 2 || isempty(type)
    type = 1;
end

exchangeRxns = findExcRxns(model);
model.lb(exchangeRxns) = 0;
model.ub(exchangeRxns) = 1000;

minimal = {'Ex_nh4'; ...   % ammonium
           'Ex_o2'; ...    % oxygen
           'Ex_pi'; ...    % phosphate
           'Ex_so4'; ...   % sulphate
           'Ex_fe2'; ...   % iron
           'Ex_h'; ...     % hydrogen
           'Ex_h2o'; ...   % water
           'Ex_na1'; ...   % sodium
           'Ex_k'; ...     % potassium
           'Ex_co2'};      % carbon dioxide

ynbVitamins = {'Ex_btn'; ...     % biotin
               'Ex_thm'; ...     % thiamine
               'Ex_4abz'; ...    % 4-aminobenzoate
               'Ex_pnto_R'; ...  % pantothenate
               'Ex_inost'; ...   % myo-inositol
               'Ex_nac'; ...     % nicotinate
               'Ex_ribflv'};     % riboflavin

coreAminoAcids = {'Ex_arg_L'; ...
                  'Ex_asp_L'; ...
                  'Ex_glu_L'; ...
                  'Ex_gly'; ...
                  'Ex_his_L'; ...
                  'Ex_ile_L'; ...
                  'Ex_leu_L'; ...
                  'Ex_lys_L'; ...
                  'Ex_met_L'; ...
                  'Ex_phe_L'; ...
                  'Ex_thr_L'; ...
                  'Ex_trp_L'; ...
                  'Ex_tyr_L'; ...
                  'Ex_val_L'; ...
                  'Ex_ura'};

allAminoAcids = [coreAminoAcids; ...
                 {'Ex_ala_L'; ...
                  'Ex_asn_L'; ...
                  'Ex_cys_L'; ...
                  'Ex_gln_L'; ...
                  'Ex_pro_L'; ...
                  'Ex_ser_L'}];

setLowerBounds(model, minimal, -1000);

if type == 2
    model = setLowerBounds(model, ynbVitamins, -2);
elseif type == 3
    model = setLowerBounds(model, ynbVitamins, -2);
    model = setLowerBounds(model, coreAminoAcids, -0.08);
elseif type == 4
    model = setLowerBounds(model, ynbVitamins, -2);
    model = setLowerBounds(model, coreAminoAcids, -0.08);
elseif type == 5
    model = setLowerBounds(model, ynbVitamins, -2);
    model = setLowerBounds(model, allAminoAcids, -0.08);
end

end

function model = setLowerBounds(model, rxns, lowerBound)
idx = findRxnIDs(model, rxns);
missing = rxns(idx == 0);
if ~isempty(missing)
    warning('setMediaPP:MissingExchange', ...
        'Some exchange reactions were not found: %s', strjoin(missing, ', '));
end
idx = idx(idx ~= 0);
model.lb(idx) = lowerBound;
end
