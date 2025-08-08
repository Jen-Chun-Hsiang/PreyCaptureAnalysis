% Script to split gauss_est and Gauss_TF_est into 4 groups: ON-temporal, ON-nasal, OFF-temporal, OFF-nasal
% Assumes gauss_est, Gauss_TF_est, cell_type_numeric, location_type_numeric are in workspace

% Logical indices for each group
isON = cell_type_numeric == 1;
isOFF = cell_type_numeric == 0;
isTemporal = location_type_numeric == 1;
isNasal = location_type_numeric == 0;

% Group indices
ON_temporal_idx = find(isON & isTemporal);
ON_nasal_idx    = find(isON & isNasal);
OFF_temporal_idx= find(isOFF & isTemporal);
OFF_nasal_idx   = find(isOFF & isNasal);

% Split gauss_est
gauss_est_ON_temporal     = gauss_est(ON_temporal_idx, :);
gauss_est_ON_nasal        = gauss_est(ON_nasal_idx, :);
gauss_est_OFF_temporal    = gauss_est(OFF_temporal_idx, :);
gauss_est_OFF_nasal       = gauss_est(OFF_nasal_idx, :);

% Split Gauss_TF_est
Gauss_TF_est_ON_temporal  = Gauss_TF_est(ON_temporal_idx, :);
Gauss_TF_est_ON_nasal     = Gauss_TF_est(ON_nasal_idx, :);
Gauss_TF_est_OFF_temporal = Gauss_TF_est(OFF_temporal_idx, :);
Gauss_TF_est_OFF_nasal    = Gauss_TF_est(OFF_nasal_idx, :);

% Check accuracy of split
assert(size(gauss_est_ON_temporal,1) == sum(isON & isTemporal), 'ON-temporal split mismatch');
assert(size(gauss_est_ON_nasal,1)    == sum(isON & isNasal),    'ON-nasal split mismatch');
assert(size(gauss_est_OFF_temporal,1)== sum(isOFF & isTemporal),'OFF-temporal split mismatch');
assert(size(gauss_est_OFF_nasal,1)   == sum(isOFF & isNasal),   'OFF-nasal split mismatch');

assert(size(Gauss_TF_est_ON_temporal,1) == sum(isON & isTemporal), 'ON-temporal TF split mismatch');
assert(size(Gauss_TF_est_ON_nasal,1)    == sum(isON & isNasal),    'ON-nasal TF split mismatch');
assert(size(Gauss_TF_est_OFF_temporal,1)== sum(isOFF & isTemporal),'OFF-temporal TF split mismatch');
assert(size(Gauss_TF_est_OFF_nasal,1)   == sum(isOFF & isNasal),   'OFF-nasal TF split mismatch');

% Optional: print summary
fprintf('ON-temporal: %d\n', numel(ON_temporal_idx));
fprintf('ON-nasal:    %d\n', numel(ON_nasal_idx));
fprintf('OFF-temporal: %d\n', numel(OFF_temporal_idx));
fprintf('OFF-nasal:   %d\n', numel(OFF_nasal_idx));

