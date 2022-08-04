% Calculate stability and typicality within networks to understand which
% networks are more stable and typical than others

% AC March 2022

clear
% =======================================================================
%% Modify this section
% =======================================================================

% Which dataset do you want to calculate? 1, 2, or 3
which_DS = [1];

% What proportion of top performers will be averaged to create "optimal" connectome
top_percent = .05; 

% =======================================================================
%% Load data
% =======================================================================
% Loads dataset determined by which_DS

% fc_flat variable is a 1 x nRun cell where each cell contains matrix of size
% nSub x nEdges (35778)
    % in DS1, sustained attention runs are runs 1-3
    % in DS2 and DS3, 
    %   sustained attention runs are runs 1-2
    %   working memory runs are runs 3-4

% behav_ATTN (and behav_WM for Datasets 2 and 3) are vectors reflecting
% performance on sustained attention and working memory tasks

% mot variable is a 1 x nRn cell where each cell contains vector of size
% nSub x 1 reflecting in-scanner frame displacement


if which_DS == 1
    load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS1_data/gradCPT_flatmats_motion_behav_data']);
elseif which_DS == 2
    load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS2_data/Chun_flatmats_motion_behav']);
elseif which_DS == 3
    load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS3_data/HCP_flatmats_motion_behav_data']);
end


% =======================================================================
%% load 10-network assignments
% =======================================================================

load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/map268_subnetwork_062521.mat');
pre_networks(:,1) = map.oldroi;
pre_networks(:,2) = map.category;

[~,inds] = sort(pre_networks(:,1));
networks = pre_networks(inds,:);
Shen_network_labels = networks(:,2);
network_names = {'Medial Frontal', 'Frontoparietal', 'Default Mode', 'Motor', 'Visual I', 'Visual II', 'Visual Association','Limbic','Basal Ganglia','Cerebellum'};

% =======================================================================
%% Identify top performers
% =======================================================================

% identify top performers
topK = ceil(top_percent * length(behav_ATTN(~isnan(behav_ATTN))));
top_scores = maxk(behav_ATTN,topK);
top_behav_ATTN = zeros([1 length(behav_ATTN)]);
for t = 1:length(top_scores)
    where = find(behav_ATTN == top_scores(t));
    top_behav_ATTN(where) = 1;
end

topSA_log = logical(top_behav_ATTN);
disp(['SA opt N = ' num2str(sum(top_behav_ATTN))])

if which_DS > 1
    fc_topSA = tanh(nanmean([fc_flat{1}(topSA_log,:); fc_flat{2}(topSA_log,:)])); % Calculate "optimal" sustained attention (runs 1 + 2) FC matrix 

    % identify top performers WM
    topK_WM = ceil(top_percent * length(behav_WM(~isnan(behav_WM))));
    top_scores_WM = maxk(behav_WM,topK_WM);
    top_behav_WM = zeros([1 length(behav_WM)]);
    for t = 1:length(top_scores_WM)
        where_WM = find(behav_WM == top_scores_WM(t));
        top_behav_WM(where_WM) = 1;
    end

    topWM_log = logical(top_behav_WM);
    disp(['WM opt N = ' num2str(sum(top_behav_WM))])

    fc_topWM = tanh(nanmean([fc_flat{3}(topWM_log,:); fc_flat{4}(topWM_log,:)])); % Calculate optimal working memory (runs 3 + 4) FC matrix

    fc_flat_minustop = fc_flat;
    fc_flat_minustop{1}(topSA_log,:) = NaN; % remove top sustained attention performers from runs 1 and 2
    fc_flat_minustop{2}(topSA_log,:) = NaN;
    fc_flat_minustop{3}(topWM_log,:) = NaN; % remove top working memory perfomers from runs 3 + 4
    fc_flat_minustop{4}(topWM_log,:) = NaN;

else % if calculating optimal FC for DS1
    
    % create optimal sustained attention FC pattern for DS1 (runs 1-3)
    fc_topSA = tanh(nanmean([fc_flat{1}(topSA_log,:); fc_flat{2}(topSA_log,:); fc_flat{3}(topSA_log,:)])); % Do not take atanh before averaging, these are already fisher-z scored
    
    % remove top sustained attention performers from data so we don't
    % double-dip
    fc_flat_minustop = fc_flat;
    fc_flat_minustop{1}(topSA_log,:) = NaN;
    fc_flat_minustop{2}(topSA_log,:) = NaN;
    fc_flat_minustop{3}(topSA_log,:) = NaN;

end

% =======================================================================
%% Loop through networks
% =======================================================================

blankmat = zeros([268]); % create empty matrix of size

trilmat = ones([268]);
trlind = find(tril(trilmat, -1));

for net1 = 1:max(Shen_network_labels)
    disp(['Calculating for network ' num2str(net1)])
    for net2 = 1:max(Shen_network_labels)

%         num = 1:length(Shen_network_labels);
        which1 = Shen_network_labels == net1; % find nodes in relevant networks
        which2 = Shen_network_labels == net2;

        % indicate which nodes belong to relevant networks
        tmpmat = blankmat;
        tmpmat(which1,which2) = 1;
        tmpmat(which2,which1) = 1;

        % create logical for network indices in lower triangle
        net_log = tmpmat(trlind) == 1;

        % loop through participants
        for p = 1:size(fc_flat{1},1)

            % loop through runs
            for r1 = 1:4 % We are only interested in first 4 runs (SA and WM) so save computational power
                for r2 = 1:4 %

                    % remove any infinite values
                    fc_flat{r1}(isinf(fc_flat{r1})) = NaN;
                    fc_flat{r2}(isinf(fc_flat{r2})) = NaN;
                    fc_flat_minustop{r1}(isinf(fc_flat_minustop{r1})) = NaN;
                    fc_flat_minustop{r2}(isinf(fc_flat_minustop{r2})) = NaN;

                    % find network-relevant edges in participant p's run1+2
                    run1 = fc_flat{r1}(p,net_log);
                    run2 = fc_flat{r2}(p,net_log);

                    % Calculate stability between network edges across runs
                    self_corr = corr(run1', run2', 'rows','pairwise');
                    
                    % find participant's mean FC in relevant edges across runs
                    self_FC = tanh(nanmean([run1; run2])');

                    % find mean FC across all other participants
                    others1 = fc_flat{r1}(:,net_log);
                    others1(p,:) = [];
                    others2 = fc_flat{r2}(:,net_log);
                    others2(p,:) = [];
                    cat_others = cat(1, others1, others2);
                    others_FC = tanh(nanmean(cat_others)');

                    % Calculate stability and typicality in network
                    Stab{net1,net2}{r1,r2}(p,:) = self_corr;
                    Typ{net1,net2}{r1,r2}(p,1) = corr(self_FC, others_FC,'rows','pairwise');

                    % Calculate similarity to optimal network FC
                    run1_top = fc_flat_minustop{r1}(p,net_log);
                    run2_top = fc_flat_minustop{r2}(p,net_log);
                    self_FC_top = tanh(nanmean([run1_top; run2_top])');

                    fc_topSA_net = fc_topSA(net_log); % Find network indices of optimal FC pattern
                    Opt_SA{net1,net2}{r1,r2}(p,1) = corr(self_FC_top, fc_topSA_net' ,'rows','pairwise');

                    if which_DS > 1
                        fc_topWM_net = fc_topWM(net_log); % Find network indeces of optimal FC pattern
                        Opt_WM{net1,net2}{r1,r2}(p,1) = corr(self_FC_top, fc_topWM_net' ,'rows','pairwise');
                    end

                end % for r1 = 1:size(all_runs,2)
            end % for r2 = 1:size(all_runs,2)
    
        end


        % Calculate discriminability in each network
        for r1 = 1:4 %
            for r2 = 1:4 %

                Disc{net1,net2}{r1,r2} = atanh(Stab{net1,net2}{r1,r2})./atanh(Typ{net1,net2}{r1,r2});

            end
        end

        % Find average connectome feature values across participants
        if which_DS ~= 1
            Net_stab_SA(net1,net2) = tanh(nanmean(atanh(Stab{net1,net2}{1,2})));
            Net_typ_SA(net1,net2) = tanh(nanmean(atanh(Typ{net1,net2}{1,2})));
            Net_disc_SA(net1,net2) = nanmean(Disc{net1,net2}{1,2});
            Net_opt_SA(net1,net2) = tanh(nanmean(atanh(Opt_SA{net1,net2}{1,2})));

            Net_stab_WM(net1,net2) = tanh(nanmean(atanh(Stab{net1,net2}{3,4})));
            Net_typ_WM(net1,net2) = tanh(nanmean(atanh(Typ{net1,net2}{3,4})));
            Net_disc_WM(net1,net2) = nanmean(Disc{net1,net2}{3,4});
            Net_opt_WM(net1,net2) = tanh(nanmean(atanh(Opt_WM{net1,net2}{3,4})));
        end

        % Since Dataset 1 has 3 gradCPT runs, we will average values across runs
        if which_DS == 1
            Net_stab_SA(net1,net2) = tanh(nanmean(nanmean([atanh(Stab{net1,net2}{1,2}), atanh(Stab{net1,net2}{1,3}), atanh(Stab{net1,net2}{2,3})],2)));
            Net_typ_SA(net1,net2) = tanh(nanmean(nanmean([atanh(Typ{net1,net2}{1,2}), atanh(Typ{net1,net2}{1,3}), atanh(Typ{net1,net2}{2,3})],2)));
            Net_disc_SA(net1,net2) = nanmean(nanmean([Disc{net1,net2}{1,2}, Disc{net1,net2}{1,3}, Disc{net1,net2}{2,3}],2));
            Net_opt_SA(net1,net2) = tanh(nanmean(nanmean([atanh(Opt_SA{net1,net2}{1,2}), atanh(Opt_SA{net1,net2}{1,3}), atanh(Opt_SA{net1,net2}{2,3})],2)));
        end

    end % for net2 =

end % for net1

% =======================================================================
%% Save data
% =======================================================================

if which_DS == 1
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_NetworkStabTypOptDisc_forReview'],'Net_stab_SA','Net_typ_SA','Net_disc_SA','Net_opt_SA');
elseif which_DS == 2
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_NetworkStabTypOptDisc_forReview'],'Net_stab_SA','Net_typ_SA','Net_stab_WM','Net_typ_WM','Net_disc_SA','Net_opt_SA','Net_disc_WM','Net_opt_WM');
elseif which_DS == 3
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_NetworkStabTypOptDisc_forReview'],'Net_stab_SA','Net_typ_SA','Net_stab_WM','Net_typ_WM','Net_disc_SA','Net_opt_SA','Net_disc_WM','Net_opt_WM');
end





