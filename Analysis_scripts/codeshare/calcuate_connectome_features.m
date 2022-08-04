
%% Calculate individual and group-level connectome features
% ===================================================================
% ===================================================================
% input:
% % fc_flat = Dataset with nRun cells where each cell contains the lower
% triangle of functional connectivity for all subjects
% size: nSubjects x ((nNodes^2 - nNodes)/2)

% behav_scores = nSubjects x 1 matrix of performance scores (to calculate
% optimal FC pattern)

% % top_percent = decimal (0.05 = 5%) indicating what percentage of
% performers should be used to calculate the optimal FC pattern (default is top 5%)

% % optimal_task_runs = which runs should be used to create the "optimal"
% FC pattern? Vector of size [1 x nOptimalityRuns] If empty, will use the average across all runs in fc_flat

% ===================================================================
% output:
% % stability = nRun x nRun cells of Pearson's correlation of functional connectivity between
% run pairs, each cell size nSubjects x 1

% % typicality = nRun x nRun cells of Pearson's correlation of a subjects
% functional connectivity with the group average FC pattern

% % optimality = nRun x nRun cells of Pearson's correlation of a subjects
% functional connectivity with FC of the highest performer(s)

% % discriminability = nRun x nRun cells of subjects' stability /
% typicality scores
% ===================================================================
% ===================================================================




function [stability, typicality, optimality, discriminability] = calcuate_connectome_features(fc_flat, behav_scores, top_percent, optimal_task_runs)

if isempty(top_percent)
    top_percent = .05;
end

% identify top performers
topK = ceil(top_percent * length(behav_scores(~isnan(behav_scores))));
top_scores = maxk(behav_scores,topK);
top_behav = zeros([1 length(behav_scores)]);
for t = 1:length(top_scores)
    where = find(behav_scores == top_scores(t));
    top_behav(where) = 1;
end

top_log = logical(top_behav);

opt_fc = [];

if ~isempty(optimal_task_runs)
    for r = 1:length(optimal_task_runs)
        opt_fc = vertcat(opt_fc, fc_flat{optimal_task_runs(r)}(top_log,:));
    end
else
    for r = 1:length(fc_flat)
        opt_fc = vertcat(opt_fc, fc_flat{r}(top_log,:));
    end
end
fc_top = tanh(nanmean(opt_fc));

% create dataset with optimal performers removed
fc_flat_minustop = fc_flat;
for r = 1:length(fc_flat)
    fc_flat_minustop{r}(top_log,:) = NaN;
end


for p = 1:size(fc_flat{1},1)
    disp(p);
    for r1 = 1:size(fc_flat,2)
        for r2 = 1:size(fc_flat,2)

            fc_flat{r1}(isinf(fc_flat{r1})) = NaN;
            fc_flat{r2}(isinf(fc_flat{r2})) = NaN;
            fc_flat_minustop{r1}(isinf(fc_flat_minustop{r1})) = NaN;
            fc_flat_minustop{r2}(isinf(fc_flat_minustop{r2})) = NaN;


            run1 = fc_flat{r1}(p,:);
            run2 = fc_flat{r2}(p,:);

            self_corr = corr(run1', run2', 'rows','pairwise');

            self_FC = tanh(nanmean([run1; run2])');

            others1 = fc_flat{r1};
            others1(p,:) = [];
            others2 = fc_flat{r2};
            others2(p,:) = [];

            cat_others = cat(1, others1, others2);

            others_FC = tanh(nanmean(cat_others)');

            typicality{r1,r2}(p,1) = corr(self_FC, others_FC,'rows','pairwise');
            stability{r1,r2}(p,:) = self_corr;

            tmpstab(r1,r2) = self_corr;

            % Optimality analysis
            run1_top = fc_flat_minustop{r1}(p,:);
            run2_top = fc_flat_minustop{r2}(p,:);
            self_FC_top = tanh(nanmean([run1_top; run2_top])');
            optimality{r1,r2}(p,1) = corr(self_FC_top, fc_top' ,'rows','pairwise');



        end % for r1 = 1:size(all_runs,2)
    end % for r2 = 1:size(all_runs,2)
    %
end


for r1 = 1:size(fc_flat,2)
    for r2 = 1:size(fc_flat,2)
        discriminability{r1,r2} = atanh(stability{r1,r2})./atanh(typicality{r1,r2});
    end
end


