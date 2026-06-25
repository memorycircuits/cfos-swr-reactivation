function result = compute_gmm_cell_properties(custom_settings)
% Compute the three cell properties required by the GMM classifier.

if nargin < 1 || isempty(custom_settings)
    custom_settings = struct();
end

result = struct();
result.spatialCoverage = compute_spatial_coverage_for_classification(custom_settings);
result.acgMean = compute_acg_mean_for_classification(custom_settings);
result.classificationFiringRate = compute_classification_firing_rate(custom_settings);

fprintf('Computed GMM cell properties: spatial_coverage_meanOFs, acg_mean, classific_firingRate.\n');
end
