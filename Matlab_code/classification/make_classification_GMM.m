function [GMM_based_classification ] = make_classification_GMM(X_all, feature_names, CCGbased_CellClassification, outfolder)
%% Gaussian Mixture Model (GMM) Classification for Cell Types

% Ground truth labels: 1 = Excitatory, 2 = Inhibitory
labels = CCGbased_CellClassification;

if nargin < 4 || isempty(outfolder)
    outfolder = fullfile(pwd, 'CellClassification', 'figures');
end
if ~exist(outfolder, 'dir')
    mkdir(outfolder);
end

% Only include labeled samples
labeled_idx = (labels == 1 | labels == 2);
X_train = X_all(labeled_idx, :);
y_labeled = labels(labeled_idx);

if isempty(X_train) || ~any(y_labeled == 1) || ~any(y_labeled == 2)
    error('GMM classification requires CCG ground-truth labels for both classes: 1 principal/excitatory and 2 interneuron/inhibitory.');
end

% Z-score standardization using labeled data
mu = mean(X_train, 1);
sigma = std(X_train, 0, 1);
sigma(sigma == 0 | ~isfinite(sigma)) = 1;
X_all_z = (X_all - mu) ./ sigma;

% Fit 1-component GMMs to each class
X_exc = X_all_z(labels == 1, :);
X_inh = X_all_z(labels == 2, :);
gmm_exc = fitgmdist(X_exc, 1, 'RegularizationValue', 1e-6);
gmm_inh = fitgmdist(X_inh, 1, 'RegularizationValue', 1e-6);

disp('Convergence status (Excitatory):'); disp(gmm_exc.Converged)
disp('Convergence status (Inhibitory):'); disp(gmm_inh.Converged)

% Compute likelihoods for full data
likelihood_exc = pdf(gmm_exc, X_all_z);
likelihood_inh = pdf(gmm_inh, X_all_z);

% Use empirical priors from labeled data
prior_exc = mean(y_labeled == 1);
prior_inh = mean(y_labeled == 2);

% Bayes rule for posterior probability
numerator_exc = likelihood_exc * prior_exc;
numerator_inh = likelihood_inh * prior_inh;
total = numerator_exc + numerator_inh;

P_exc_given_x = numerator_exc ./ total;
P_inh_given_x = numerator_inh ./ total;
P_exc_given_x(~isfinite(P_exc_given_x)) = NaN;
P_inh_given_x(~isfinite(P_inh_given_x)) = NaN;

% Hard classification
predicted_labels = ones(size(P_exc_given_x));
predicted_labels(P_inh_given_x > P_exc_given_x) = 2;

%% Combined accuracy using ALL features (evaluate on labeled cells only)
acc_all = mean(predicted_labels(labeled_idx) == labels(labeled_idx));
fprintf('\nCombined GMM accuracy (all %d features): %.2f%% (n=%d labeled)\n', ...
    size(X_all,2), 100*acc_all, sum(labeled_idx));

% Confusion matrix (rows = GT, cols = Pred)
cm = confusionmat(labels(labeled_idx), predicted_labels(labeled_idx), 'Order', [1 2]);
disp('Confusion matrix (rows=GT [Exc; Inh], cols=Pred [Exc Inh]):');
disp(cm);

% Balanced accuracy (handles class imbalance)
rec_exc = cm(1,1) / sum(cm(1,:));
rec_inh = cm(2,2) / sum(cm(2,:));
bal_acc = mean([rec_exc, rec_inh]);
fprintf('Balanced accuracy: %.2f%% (Exc recall %.2f%%, Inh recall %.2f%%)\n', ...
    100*bal_acc, 100*rec_exc, 100*rec_inh);

% Optional: visual confusion chart (if available in your MATLAB)
%{
if exist('confusionchart','file') == 2
    figure;
    confusionchart(labels(labeled_idx), predicted_labels(labeled_idx), ...
        'Order', [1 2], 'RowSummary','row-normalized', 'ColumnSummary','column-normalized');
    title(sprintf('GMM (all %d features) confusion chart', size(X_all,2)));
end
%}


%% Export
GMM_based_classification = predicted_labels;

%% 3D Plot of Posterior Probability
%{
figure;
scatter3(X_all(:,1), X_all(:,2), X_all(:,3), 30, P_exc_given_x, 'filled');
xlabel(feature_names{1}); ylabel(feature_names{2}); zlabel(feature_names{3});
title('Posterior Probability of Excitatory (GMM)');
colorbar;

% Custom colormap
nColors = 256;
red_half = [linspace(1, 0.85, nColors/2)', linspace(0.3, 0.3, nColors/2)', linspace(0.0, 0.85, nColors/2)'];
blue_half = [linspace(0.85, 0.0, nColors/2)', linspace(0.3, 0.347, nColors/2)', linspace(0.85, 0.9, nColors/2)'];
redToBlueNeutral = [red_half; blue_half];
colormap(redToBlueNeutral);
caxis([0 1]);
    filename = fullfile(outfolder, sprintf('Posterior Probability of GMM E and I classification.pdf')) ;
    exportgraphics(gcf, filename, 'ContentType', 'vector'); % High-quality PDF output
%}

%% 3D posterior plot with hybrid vector/raster export
f = figure('Color','w', ...
    'DefaultAxesFontName','Arial', ...
    'DefaultAxesFontSize',20, ...
    'DefaultAxesFontWeight','bold', ...
    'DefaultTextFontName','Arial', ...
    'DefaultTextFontSize',20,...
    'DefaultTextFontWeight','bold');

ax = axes(f, ...
    'FontName','Arial', ...
    'DefaultAxesFontWeight','bold', ...
    'FontSize',20, ... 
    'LineWidth', 3);
hold(ax, 'on');
%ax.LineWidth = 3;
% Base posterior scatter (this should become raster in PDF)
s = scatter3(ax, ...
    X_all(:,1), X_all(:,2), X_all(:,3), ...
    20, P_exc_given_x, 'filled', ...
    'MarkerFaceAlpha', 0.65, ...   % encourages rasterization
    'MarkerEdgeAlpha', 0.15);

xlabel(ax, feature_names{1});
ylabel(ax, feature_names{2});
zlabel(ax, feature_names{3});
title(ax, 'Posterior Probability of Excitatory (GMM)');

% Custom colormap
nColors = 256;
red_half = [linspace(1, 0.85, nColors/2)', ...
            linspace(0.3, 0.3, nColors/2)', ...
            linspace(0.0, 0.85, nColors/2)'];
blue_half = [linspace(0.85, 0.0, nColors/2)', ...
             linspace(0.3, 0.347, nColors/2)', ...
             linspace(0.85, 0.9, nColors/2)'];
redToBlueNeutral = [red_half; blue_half];
colormap(ax, redToBlueNeutral);
clim(ax, [0 1]);

%cb = colorbar(ax);
%cb.Label.String = 'P(excitatory | x)';

grid(ax, 'on');
axis(ax, 'tight');
view(ax, 3);

% Optional: set log scales if relevant for the chosen 3 features
for k = 1:3
    if strcmp(feature_names{k}, 'acg tau rise') || strcmp(feature_names{k}, 'acg tau decay')
        switch k
            case 1, set(ax, 'XScale', 'log');
            case 2, set(ax, 'YScale', 'log');
            case 3, set(ax, 'ZScale', 'log');
        end
    end
end

%% Ground-truth annotation on top of posterior cloud
% labels == 1 : GT excitatory
% labels == 2 : GT inhibitory

% Use only a subset if too crowded
idxE = find(labels == 1);
idxI = find(labels == 2);

% Example: plot all GT points, but larger and unfilled so they remain visible
scatter3(ax, ...
    X_all(idxE,1), X_all(idxE,2), X_all(idxE,3), ...
    70, '^', ...
    'MarkerEdgeColor', [0 0 0], ...
    'LineWidth', 0.9, ...
    'MarkerFaceColor', 'none');

scatter3(ax, ...
    X_all(idxI,1), X_all(idxI,2), X_all(idxI,3), ...
    70, 'o', ...
    'MarkerEdgeColor', [0 0 0], ...
    'LineWidth', 0.9, ...
    'MarkerFaceColor', 'none');

%legend(ax, ...
%    {'Posterior probability', 'GT excitatory', 'GT inhibitory'}, ...
%    'Location', 'northeastoutside');

f.Position = [100 100 700 600];

%% Export
filename = fullfile(outfolder, 'Posterior_Probability_GMM_EI_classification.pdf');
exportgraphics(f, filename, ...
    'ContentType', 'auto', ...
    'BackgroundColor', 'white', ...
    'Resolution', 300);

cbFig = figure('Color','w', 'Units','pixels', 'Position',[100 100 180 700]);

cbAx = axes('Parent', cbFig, ...
    'Position',[0.10 0.05 0.10 0.90], ...
    'Visible','off');

% Dummy image just to host the colormap / clim
imagesc(cbAx, [0 1; 0 1]);
set(cbAx, 'Visible', 'off');

colormap(cbAx, redToBlueNeutral);
clim(cbAx, [0 1]);

cb = colorbar(cbAx, 'Position',[0.38 0.08 0.22 0.84]);
cb.Label.String = 'P(excitatory | x)';
cb.TickDirection = 'out';
cb.Box = 'off';
cb.FontSize = 30 ;
cb.Label.FontSize = 30;

filename_cb = fullfile(outfolder, 'Posterior_GMM_colorbar_only.pdf');
exportgraphics(cbFig, filename_cb, 'ContentType', 'vector');


%% Pairwise 2D Feature Plots with Predictions
%figure;
%plot_num = 1;
num_features = size(X_all, 2);
gt_exc_color = [0.00 0.16 0.55];
gt_inh_color = [0.62 0.00 0.00];
for i = 1:num_features
    for j = i+1:num_features
        %subplot(2, 3, plot_num);
        f1 = figure;
        h_pred = gscatter(X_all(:,i), X_all(:,j), predicted_labels, 'br', '..'); hold on
        xlabel(feature_names{i});
        ylabel(feature_names{j});
        if strcmp('acg tau rise', feature_names{i}) || strcmp('acg tau decay', feature_names{i})
            set(gca, 'XScale', 'log');
        elseif strcmp('acg tau rise', feature_names{j})  || strcmp('acg tau decay', feature_names{j})
            set(gca, 'YScale', 'log');
        end
        
        title(sprintf('GMM Classification: %s vs %s', feature_names{i}, feature_names{j}));
       
        f1.Position = [100 100 500 500];
        axis square; grid on;

        h_gt_exc = scatter(X_all(labels == 1, i), X_all(labels == 1, j), ...
            55, ...
            'Marker', '^', ...
            'MarkerEdgeColor', gt_exc_color, ...
            'MarkerFaceColor', 'none', ...
            'LineWidth', 1.3);
        h_gt_inh = scatter(X_all(labels == 2, i), X_all(labels == 2, j), ...
            55, ...
            'Marker', 'o', ...
            'MarkerEdgeColor', gt_inh_color, ...
            'MarkerFaceColor', 'none', ...
            'LineWidth', 1.3);
        uistack([h_gt_exc h_gt_inh], 'top');

        legend([h_pred(:); h_gt_exc; h_gt_inh], ...
            {'excitatory','inhibitory', 'GT excitatory','GT inhibitory'}, ...
            'Location','northwestoutside')
        
        filename = fullfile(outfolder, sprintf('GMM_Classification_%s_%s.pdf', feature_names{i}, feature_names{j}));
        exportgraphics(gcf, filename); 
        
        %legend('Location','eastoutside')
        %plot_num = plot_num + 1;
    end
end

%% Feature Importance (1D GMM Accuracy)
fprintf('\nFeature importance (1D GMM classification accuracy):\n');
for f = 1:num_features
    Xf = X_all_z(:, f);
    Xf_exc = Xf(labels == 1);
    Xf_inh = Xf(labels == 2);

    gmm_f_exc = fitgmdist(Xf_exc, 1, 'RegularizationValue', 1e-6);
    gmm_f_inh = fitgmdist(Xf_inh, 1, 'RegularizationValue', 1e-6);

    l_exc = pdf(gmm_f_exc, Xf);
    l_inh = pdf(gmm_f_inh, Xf);

    pred_f = ones(size(Xf));
    pred_f(l_inh > l_exc) = 2;

    acc = mean(pred_f(labeled_idx) == labels(labeled_idx));
    fprintf('%s: %.2f%%\n', feature_names{f}, acc * 100);
end

%% 2-Feature GMM Accuracy
fprintf('\nFeature combination GMM classification accuracy (2D):\n');
combs = nchoosek(1:num_features, 2);  % All 2-feature combinations

for i = 1:size(combs, 1)
    idx = combs(i, :);
    X_sub = X_all_z(:, idx);
    X_sub_exc = X_sub(labels == 1, :);
    X_sub_inh = X_sub(labels == 2, :);

    gmm_exc_sub = fitgmdist(X_sub_exc, 1, 'RegularizationValue', 1e-6);
    gmm_inh_sub = fitgmdist(X_sub_inh, 1, 'RegularizationValue', 1e-6);

    l_exc = pdf(gmm_exc_sub, X_sub(labeled_idx, :));
    l_inh = pdf(gmm_inh_sub, X_sub(labeled_idx, :));

    pred = ones(sum(labeled_idx), 1);
    pred(l_inh > l_exc) = 2;

    acc = mean(pred == labels(labeled_idx));
    fprintf('%s + %s: %.2f%%\n', feature_names{idx(1)}, feature_names{idx(2)}, acc * 100);
end



end

