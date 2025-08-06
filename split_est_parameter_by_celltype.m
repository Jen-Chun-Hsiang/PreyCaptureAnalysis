% after loading and vectorize data from
% "WhiteNoise_ONOFFalpha_Comparison.m"
ON_nasal_ids = cell_type_numeric == 1 & location_type_numeric == 0;
ON_temporal_ids = cell_type_numeric == 1 & location_type_numeric == 1;
OFF_nasal_ids = cell_type_numeric == 0 & location_type_numeric == 0;
OFF_temporal_ids = cell_type_numeric == 0 & location_type_numeric == 1;

ON_nasal_params = gauss_est(ON_nasal_ids, :);
ON_temporal_params = gauss_est(ON_temporal_ids, :);
OFF_nasal_params = gauss_est(OFF_nasal_ids, :);
OFF_temporal_params = gauss_est(OFF_temporal_ids, :);

