% Iterate null relationships between connectome stability and performance
% Computationally lesion networks in the brain and calculate correlation btw stability and performance
% Goal is to determine whether canonical networks significantly contribute
% to relationship observed

% AC March 2022

clear
% =======================================================================
%% Modify this section
% =======================================================================

% Which dataset do you want to calculate? 1, 2, or 3
which_DS = [2];

% How many iterations?
IT = 250;

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
%% Shuffle lesioned network
% =======================================================================

blankmat = ones([268]);
trlind = find(tril(blankmat, -1));

for iter = 1:IT % iterate ___ times (set at top)

    for net = 1:max(Shen_network_labels)
        disp(['Iteration ' num2str(iter) ' network ' num2str(net)]) ;

        num = 1:length(Shen_network_labels);
        which = find(Shen_network_labels == net); % which indices are actually the network(want to make sure we don't include these in MC stats)

        num(which) = NaN; % delete those as possible options
        inds = num(~isnan(num));

        inds = inds(randperm(length(inds))); %shuffle so different indices end up at the beginning of the vector

        log = inds(1:length(which)); % grab the same number of indices as the network


        tmpmat = blankmat;
        tmpmat(log,:) = NaN;
        tmpmat(:,log) = NaN;

        lesion_log = tmpmat(trlind) == 1;

        for p = 1:size(fc_flat{1},1) % loop over participants

            for r1 = 1:4 % We are only interested in first 4 runs (SA and WM) so save computational power
                for r2 = 1:4 %

                    fc_flat{r1}(isinf(fc_flat{r1})) = NaN;
                    fc_flat{r2}(isinf(fc_flat{r2})) = NaN;
                    run1 = fc_flat{r1}(p,lesion_log);
                    run2 = fc_flat{r2}(p,lesion_log);

                    self_corr = corr(run1', run2', 'rows','pairwise');

                    self_FC = nanmean([run1; run2])';


                    others1 = fc_flat{r1}(:,lesion_log);
                    others1(p,:) = [];
                    others2 = fc_flat{r2}(:,lesion_log);
                    others2(p,:) = [];
                    cat_others = cat(1, others1, others2);
                    others_FC = nanmean(cat_others)';

                    Stab{r1,r2}(p,:) = self_corr;
                    Typ{r1,r2}(p,1) = corr(self_FC, others_FC,'rows','pairwise');

                end % for r1 = 1:size(all_runs,2)
            end % for r2 = 1:size(all_runs,2)
        end % for p =

        if which_DS ~= 1

            for r1 = 1%4%size(Stab,2)
                for r2 = 2%:4%size(Stab,2)

                    stab = Stab{r1,r2};
                    typ = Typ{r1,r2};

                    [corr_stab_ATTN{net,iter}(r1,r2), sig_stab_ATTN{net,iter}(r1,r2)] = partialcorr(stab, behav_ATTN, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
                    [corr_typ_ATTN{net,iter}(r1,r2), sig_typ_ATTN{net,iter}(r1,r2)] = partialcorr(typ, behav_ATTN, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
                end
            end


            for r1 = 3
                for r2 = 4
                    stab = Stab{r1,r2};
                    typ = Typ{r1,r2};

                    [corr_stab_WM{net,iter}(r1,r2), sig_stab_WM{net,iter}(r1,r2)] = partialcorr(stab, behav_WM, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
                    [corr_typ_WM{net,iter}(r1,r2), sig_typ_WM{net,iter}(r1,r2)] = partialcorr(typ, behav_WM, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');

                end
            end
        end

        % Since Dataset 1 has 3 gradCPT runs, we will average values across runs
        if which_DS == 1
            Avg_Stab = tanh(nanmean([atanh(Stab{1,2}), atanh(Stab{1,3}), atanh(Stab{2,3})],2));
            Avg_Typ = tanh(nanmean([atanh(Typ{1,2}), atanh(Typ{1,3}), atanh(Typ{2,3})],2));

            [r_stab_ds1{net,iter},p_stab_ds1{net,iter}] = partialcorr(Avg_Stab, behav_ATTN, [nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{1} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');
            [r_typ_ds1{net,iter},p_typ_ds1{net,iter}] = partialcorr(Avg_Typ, behav_ATTN, [nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{1} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');

        end



        % =======================================================================
        %% Save
        % =======================================================================

        if which_DS == 1
            save(['//Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_NullStabTypAnat'],'Stab','r_stab_ds1','p_stab_ds1','Avg_Stab','r_typ_ds1','p_typ_ds1','Avg_Typ');
        elseif which_DS == 2
            save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_NullStabTypAnat'],'Stab','Typ','corr_typ_ATTN','sig_typ_ATTN','corr_typ_WM','sig_typ_WM','corr_stab_ATTN', 'sig_stab_ATTN', 'corr_stab_WM', 'sig_stab_WM');
        elseif which_DS == 3
            save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_NullStabTypAnat'],'Stab','Typ','corr_typ_ATTN','sig_typ_ATTN','corr_typ_WM','sig_typ_WM','corr_stab_ATTN', 'sig_stab_ATTN', 'corr_stab_WM', 'sig_stab_WM');
        end



    end % for net =
end % for iter = IT