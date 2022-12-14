%% Calculate stability, typicality, optimality, and discriminability for each dataset
% Final version written July 2022


% =======================================================================
%% Modify this section
% =======================================================================
clear
% Which dataset do you want to calculate? 1, 2, or 3
which_DS = [2];

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
% edges are Fisher-z values

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
    fc_topSA = nanmean([fc_flat{1}(topSA_log,:); fc_flat{2}(topSA_log,:)]); % Calculate "optimal" sustained attention (runs 1 + 2) FC matrix 

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

    fc_topWM = nanmean([fc_flat{3}(topWM_log,:); fc_flat{4}(topWM_log,:)]); % Calculate optimal working memory (runs 3 + 4) FC matrix

    fc_flat_minustop = fc_flat;
    fc_flat_minustop{1}(topSA_log,:) = NaN; % remove top sustained attention performers from runs 1 and 2
    fc_flat_minustop{2}(topSA_log,:) = NaN;
    fc_flat_minustop{3}(topWM_log,:) = NaN; % remove top working memory perfomers from runs 3 + 4
    fc_flat_minustop{4}(topWM_log,:) = NaN;

    % calculate difference from mean performance (performance "typicality")
    behav_diff_ATTN = abs(behav_ATTN - nanmean(behav_ATTN));
    behav_diff_WM = abs(behav_WM - nanmean(behav_WM));

else % if calculating optimal FC for DS1
    
    % create optimal sustained attention FC pattern for DS1 (runs 1-3)
    fc_topSA = nanmean([fc_flat{1}(topSA_log,:); fc_flat{2}(topSA_log,:); fc_flat{3}(topSA_log,:)]); % Do not take atanh before averaging, these are already fisher-z scored
    
    % remove top sustained attention performers from data so we don't
    % double-dip
    fc_flat_minustop = fc_flat;
    fc_flat_minustop{1}(topSA_log,:) = NaN;
    fc_flat_minustop{2}(topSA_log,:) = NaN;
    fc_flat_minustop{3}(topSA_log,:) = NaN;

    % calculate difference from mean performance (performance "typicality")
    behav_diff_ATTN = abs(behav_ATTN - nanmean(behav_ATTN));
end



% =======================================================================
%% Calculate individual and group-level connectome features
% =======================================================================

trlinds = find(tril(ones(size(fc_flat,2)),-1)); % calculate indices of lower triangle in nRun x nRun matrix (used to calculate overall stability later)

for p = 1:size(fc_flat{1},1) % p represents 1 participant
    disp(['running participant ' num2str(p)]);

    for r1 = 1:size(fc_flat,2) % loop through runs to calculate connectome features 
        for r2 = 1:size(fc_flat,2) 

            % remove any infinite values
            fc_flat{r1}(isinf(fc_flat{r1})) = NaN; 
            fc_flat{r2}(isinf(fc_flat{r2})) = NaN;
            fc_flat_minustop{r1}(isinf(fc_flat_minustop{r1})) = NaN;
            fc_flat_minustop{r2}(isinf(fc_flat_minustop{r2})) = NaN;

            % participant p's run r1 and run r2 FC pattern
            run1 = fc_flat{r1}(p,:);
            run2 = fc_flat{r2}(p,:);

            % correlate participant p's FC patterns (this is stability)
            self_corr = corr(run1', run2', 'rows','pairwise');

            % calculate mean FC pattern between runs
            self_FC = nanmean([run1; run2])';

            % calculate average FC pattern of all participants except ppt p
            others1 = fc_flat{r1};
            others1(p,:) = [];
            others2 = fc_flat{r2};
            others2(p,:) = [];
            cat_others = cat(1, others1, others2);
            others_FC = nanmean(cat_others)';

            % Calculate connectome features stability and typicality
            Stab{r1,r2}(p,:) = self_corr;
            Typ{r1,r2}(p,1) = corr(self_FC, others_FC,'rows','pairwise');
            
            % record between run stability to calculate overall stability later
            tmpstab(r1,r2) = self_corr;

            % Calculate optimality using variables with top performers removed
            run1_top = fc_flat_minustop{r1}(p,:);
            run2_top = fc_flat_minustop{r2}(p,:);
            self_FC_top = nanmean([run1_top; run2_top])';
            Opt_SA{r1,r2}(p,1) = corr(self_FC_top, fc_topSA' ,'rows','pairwise');

            if which_DS > 1 % if not using DS1, calculate optimality to working memory FC
                Opt_WM{r1,r2}(p,1) = corr(self_FC_top, fc_topWM' ,'rows','pairwise');
            end

        end % for r1 = 1:size(all_runs,2)
    end % for r2 = 1:size(all_runs,2)

    % Calculate overall stability 
    stab_fishz = atanh(tmpstab);
    OverStab(p,:) = tanh(nanmean(stab_fishz(trlinds))); % average stability overall

    if which_DS == 3 % If using DS3 which has 4 rest runs, calculate mean rest stability and typicality (DS1 and DS2 have 2 only rest runs, therefore only one value for rest stability and typicality) 
        RestStab(p,:) = tanh(nanmean([stab_fishz(17,18),stab_fishz(17,19),stab_fishz(17,20),stab_fishz(18,19),stab_fishz(18,20),stab_fishz(19,20)]));
        RestTyp(p,:) = tanh(nanmean(atanh([Typ{17,18}(p), Typ{17,19}(p), Typ{17,20}(p), Typ{18,19}(p), Typ{18,20}(p), Typ{19,20}(p)]),2));
    end % if which_DS == 3
end % for p = 1:size(fc_flat{1},1)


% Calculate discriminability 
for r1 = 1:size(fc_flat,2)
    for r2 = 1:size(fc_flat,2)
        Disc{r1,r2} = atanh(Stab{r1,r2})./atanh(Typ{r1,r2});
    end
end

% =======================================================================
%% Relate connectome features to performance
% =======================================================================
% controlling for motion 

for r1 = 1:size(Disc,2)
    for r2 = 1:size(Disc,2)

        stab = Stab{r1,r2};
        typ = Typ{r1,r2};
        disc = Disc{r1,r2};
        optSA = Opt_SA{r1,r2};


        [corr_stab_ATTN(r1,r2), sig_stab_ATTN(r1,r2)] = partialcorr(stab, behav_ATTN, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
        [corr_typ_ATTN(r1,r2), sig_typ_ATTN(r1,r2)] = partialcorr(typ, behav_ATTN, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
        [corr_disc_ATTN(r1,r2), sig_disc_ATTN(r1,r2)] = partialcorr(disc, behav_ATTN, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
        [corr_opt_ATTN(r1,r2), sig_opt_ATTN(r1,r2)] = partialcorr(optSA, behav_ATTN, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');


        if which_DS ~= 1 % 

            optWM = Opt_WM{r1,r2};

            [corr_stab_WM(r1,r2), sig_stab_WM(r1,r2)] = partialcorr(stab, behav_WM, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
            [corr_typ_WM(r1,r2), sig_typ_WM(r1,r2)] = partialcorr(typ, behav_WM, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
            [corr_disc_WM(r1,r2), sig_disc_WM(r1,r2)] = partialcorr(disc, behav_WM, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
            [corr_opt_WM(r1,r2), sig_opt_WM(r1,r2)] = partialcorr(optWM, behav_WM, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');

            % Currently, only have run-specific data for Dataset 2, so only do this consistency analysis for DS2
            if which_DS == 2
                [corr_run1_ATTN(r1,r2), sig_run1_ATTN(r1,r2)] = partialcorr(stab, behav_ATTN_run1, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
                [corr_run2_ATTN(r1,r2), sig_run2_ATTN(r1,r2)] = partialcorr(stab, behav_ATTN_run2, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
                [corr_run1_WM(r1,r2), sig_run1_WM(r1,r2)] = partialcorr(stab, behav_WM_run1, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
                [corr_run2_WM(r1,r2), sig_run2_WM(r1,r2)] = partialcorr(stab, behav_WM_run2, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');

                [corr_consistencycontrol_ATTN(r1,r2), sig_consistencycontrol_ATTN(r1,r2)] = partialcorr(stab, behav_ATTN, [abs(behav_ATTN_run1 - behav_ATTN_run2), nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
                [corr_consistencycontrol_WM(r1,r2), sig_consistencycontrol_WM(r1,r2)] = partialcorr(stab, behav_WM, [abs(behav_WM_run1 - behav_WM_run2), nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');

            end

            % calculate similarity between sustained attention and working
            % memory stability
            [corr_attn_wm, sig_attn_wm] = partialcorr(Stab{1,2},Stab{3,4},[nanmean([mot{1}, mot{2}, mot{3}, mot{4}], 2), abs(nanmean([mot{1}, mot{2}],2) - nanmean([mot{3},mot{4}],2))],'rows','pairwise');
            

            % calculate correlation between typicality and performance
            % typicality
            [corr_typ_diff_ATTN(r1,r2), sig_typ_diff_ATTN(r1,r2)] = partialcorr(typ, behav_diff_ATTN, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
            [corr_typ_diff_WM(r1,r2), sig_typ_diff_WM(r1,r2)] = partialcorr(typ, behav_diff_WM, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
          
        end % if which_DS ~= 1

    end % for r2
end % for r1

% Since Dataset 1 has 3 gradCPT runs, we will average values across runs
if which_DS == 1
    Avg_Stab = tanh(nanmean([atanh(Stab{1,2}), atanh(Stab{1,3}), atanh(Stab{2,3})],2));
    Avg_Typ = tanh(nanmean([atanh(Typ{1,2}), atanh(Typ{1,3}), atanh(Typ{2,3})],2));
    Avg_Disc = nanmean([Disc{1,2}, Disc{1,3}, Disc{2,3}],2);
    Avg_Opt = tanh(nanmean([atanh(Opt_SA{1,2}), atanh(Opt_SA{1,3}), atanh(Opt_SA{2,3})],2));

    [r_stab_ds1,p_stab_ds1] = partialcorr(Avg_Stab, behav_ATTN, [nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{1} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');
    [r_typ_ds1,p_typ_ds1] = partialcorr(Avg_Typ, behav_ATTN, [nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{1} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');
    [r_disc_ds1,p_disc_ds1] = partialcorr(Avg_Disc, behav_ATTN, [nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{1} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');
    [r_opt_ds1,p_opt_ds1] = partialcorr(Avg_Opt, behav_ATTN, [nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{1} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');


    % Does Avg stability correlate with individual run performance in dataset 1?
    [corr_run1_ATTN, sig_run1_ATTN] = partialcorr(Avg_Stab, behav_ATTN_run1, [nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{1} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');
    [corr_run2_ATTN, sig_run2_ATTN] = partialcorr(Avg_Stab, behav_ATTN_run2, [nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{1} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');
    [corr_run3_ATTN, sig_run3_ATTN] = partialcorr(Avg_Stab, behav_ATTN_run3, [nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{1} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');

    % Does Avg Typicality correlate with performance typicality? *difference from mean
    % performance*
    [r_diffmean_ds1,p_diffmean_ds1] = partialcorr(Avg_Typ, behav_diff_ATTN, [nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{2} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');

end

% =======================================================================
%% Relate connectome features to performance
% =======================================================================

if which_DS == 3
    [corr_rest_ATTN, sig_rest_ATTN] = partialcorr(RestStab, behav_ATTN,[nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
    [corr_rest_WM, sig_rest_WM] = partialcorr(RestStab, behav_WM,[nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');

    [corr_resttyp_ATTN, sig_resttyp_ATTN] = partialcorr(RestTyp, behav_ATTN,[nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');
    [corr_resttyp_WM, sig_resttyp_WM] = partialcorr(RestTyp, behav_WM,[nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');

end % if which_DS == 3


% =======================================================================
%% Relate connectome features to each other
% =======================================================================
if which_DS == 1

    for m = 1:size(mot,2)
        avg_mot_pre(:,m) = mot{m};
    end
    avg_mot = nanmean(avg_mot_pre,2);

    % calculate relationship between overall stability and performance
    [OverCorr, sig_OverCorr] = partialcorr(OverStab, behav_ATTN, [avg_mot],'rows','pairwise','type','Spearman');

    % control for mean difference in performance
    abs_diff_perf = nanmean([abs(behav_ATTN_run1 - behav_ATTN_run2), abs(behav_ATTN_run1 - behav_ATTN_run3), abs(behav_ATTN_run3 - behav_ATTN_run2)],2);
    [corr_consistencycontrol, p_consistencycontrol] = partialcorr(Avg_Stab, behav_ATTN, [abs_diff_perf, nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{1} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');

    % correlation between typicality and performance typicality (difference
    % from mean performance)
    behav_diff_ATTN = abs(behav_ATTN - nanmean(behav_ATTN));
    [r_diffmean_ds1,p_diffmean_ds1] = partialcorr(Avg_Typ, behav_diff_ATTN, [nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{2} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');

    % Spearman correlation between connectome features
    [rST_ATTN,pST_ATTN] = corr(Avg_Stab,Avg_Typ,'rows','pairwise','type','Spearman');
    [rSO_ATTN,pSO_ATTN] = corr(Avg_Stab,Avg_Opt,'rows','pairwise','type','Spearman');
    [rSD_ATTN,pSD_ATTN] = corr(Avg_Stab,Avg_Disc,'rows','pairwise','type','Spearman');
    [rTO_ATTN,pTO_ATTN] = corr(Avg_Typ,Avg_Opt,'rows','pairwise','type','Spearman');
    [rTD_ATTN,pTD_ATTN] = corr(Avg_Typ,Avg_Disc,'rows','pairwise','type','Spearman');
    [rOD_ATTN,pOD_ATTN] = corr(Avg_Opt,Avg_Disc,'rows','pairwise','type','Spearman');

elseif which_DS == 2

    for m = 1:size(mot,2)
        avg_mot_pre(:,m) = mot{m};
    end
    avg_mot = nanmean(avg_mot_pre,2);

    [OverCorr_ATTN, sig_OverCorr_ATTN] = partialcorr(OverStab, behav_ATTN, [avg_mot],'rows','pairwise','type','Spearman');
    [OverCorr_WM, sig_OverCorr_WM] = partialcorr(OverStab, behav_WM, [avg_mot],'rows','pairwise','type','Spearman');

    % how much do WM and SA stability values correlate?
    [r_SAWM, p_SAWM] = partialcorr(Stab{1,2},Stab{3,4},[nanmean([mot{1} mot{2}],2), nanmean([mot{3} mot{4}],2), abs(nanmean([mot{1} mot{2}],2) - nanmean([mot{3} mot{4}],2))],'rows','pairwise','type','Spearman');

    % Correlation between typicality and performance typicality
    behav_diff_ATTN = abs(behav_ATTN - nanmean(behav_ATTN));
    [r_diffmean_ATTN_ds2,p_diffmean_ATTN_ds2] = partialcorr(Typ{1,2}, behav_diff_ATTN, [nanmean([mot{1}, mot{2}],2), abs(mot{1}-mot{2})],'rows','pairwise','type','Spearman');
    behav_diff_WM = abs(behav_WM - nanmean(behav_WM));
    [r_diffmean_WM_ds2,p_diffmean_WM_ds2] = partialcorr(Typ{3,4}, behav_diff_WM, [nanmean([mot{3}, mot{4}],2), abs(mot{3}-mot{4})],'rows','pairwise','type','Spearman');

    % Spearman correlation between connectome features
    [rST_ATTN,pST_ATTN] = corr(Stab{1,2},Typ{1,2},'rows','pairwise','type','Spearman');
    [rST_WM,pST_WM] = corr(Stab{3,4},Typ{3,4},'rows','pairwise','type','Spearman');

    [rSO_ATTN,pSO_ATTN] = corr(Stab{1,2},Opt_SA{1,2},'rows','pairwise','type','Spearman');
    [rSO_WM,pSO_WM] = corr(Stab{3,4},Opt_WM{3,4},'rows','pairwise','type','Spearman');

    [rSD_ATTN,pSD_ATTN] = corr(Stab{1,2},Disc{1,2},'rows','pairwise','type','Spearman');
    [rSD_WM,pSD_WM] = corr(Stab{3,4},Disc{3,4},'rows','pairwise','type','Spearman');

    [rTO_ATTN,pTO_ATTN] = corr(Typ{1,2},Opt_SA{1,2},'rows','pairwise','type','Spearman');
    [rTO_WM,pTO_WM] = corr(Typ{3,4},Opt_WM{3,4},'rows','pairwise','type','Spearman');

    [rTD_ATTN,pTD_ATTN] = corr(Typ{1,2},Disc{1,2},'rows','pairwise','type','Spearman');
    [rTD_WM,pTD_WM] = corr(Typ{3,4},Disc{3,4},'rows','pairwise','type','Spearman');

    [rOD_ATTN,pOD_ATTN] = corr(Opt_SA{1,2},Disc{1,2},'rows','pairwise','type','Spearman');
    [rOD_WM,pOD_WM] = corr(Opt_WM{3,4},Disc{3,4},'rows','pairwise','type','Spearman');


elseif which_DS == 3

    for m = 1:size(mot,2)
        avg_mot_pre(:,m) = mot{m};
    end
    avg_mot = nanmean(avg_mot_pre,2);

    % correlation between overall stability and performance
    [OverCorr_ATTN, sig_OverCorr_ATTN] = partialcorr(OverStab, behav_ATTN, [avg_mot],'rows','pairwise','type','Spearman');
    [OverCorr_WM, sig_OverCorr_WM] = partialcorr(OverStab, behav_WM, [avg_mot],'rows','pairwise','type','Spearman');

    % correlation between rest stability and performance
    [OverRest_ATTN, sig_OverRest_ATTN] = partialcorr(RestStab, behav_ATTN, [avg_mot],'rows','pairwise','type','Spearman');
    [OverRest_typ_ATTN, sig_OverRest_typ_ATTN] = partialcorr(RestTyp, behav_ATTN, [avg_mot],'rows','pairwise','type','Spearman');

    [OverRest_WM, sig_OverRest_WM] = partialcorr(RestStab, behav_WM, [avg_mot],'rows','pairwise','type','Spearman');
    [OverRest_typ_WM, sig_OverRest_typ_WM] = partialcorr(RestTyp, behav_WM, [avg_mot],'rows','pairwise','type','Spearman');

   
    % how much do WM and SA stability values correlate?
    [r_SAWM, p_SAWM] = partialcorr(Stab{1,2},Stab{3,4},[nanmean([mot{1} mot{2}],2), nanmean([mot{3} mot{4}],2), abs(nanmean([mot{1} mot{2}],2) - nanmean([mot{3} mot{4}],2))],'rows','pairwise','type','Spearman');

    % Correlation between typicality and performance typicality
    behav_diff_ATTN = abs(behav_ATTN - nanmean(behav_ATTN));
    [r_diffmean_ATTN_ds3,p_diffmean_ATTN_ds3] = partialcorr(Typ{1,2}, behav_diff_ATTN, [nanmean([mot{1}, mot{2}],2), abs(mot{1}-mot{2})],'rows','pairwise','type','Spearman');
    behav_diff_WM = abs(behav_WM - nanmean(behav_WM));
    [r_diffmean_WM_ds3,p_diffmean_WM_ds3] = partialcorr(Typ{3,4}, behav_diff_WM, [nanmean([mot{3}, mot{4}],2), abs(mot{3}-mot{4})],'rows','pairwise','type','Spearman');

    % Spearman correlation between connectome features
    [rST_ATTN,pST_ATTN] = corr(Stab{1,2},Typ{1,2},'rows','pairwise','type','Spearman');
    [rST_WM,pST_WM] = corr(Stab{3,4},Typ{3,4},'rows','pairwise','type','Spearman');

    [rSO_ATTN,pSO_ATTN] = corr(Stab{1,2},Opt_SA{1,2},'rows','pairwise','type','Spearman');
    [rSO_WM,pSO_WM] = corr(Stab{3,4},Opt_WM{3,4},'rows','pairwise','type','Spearman');

    [rSD_ATTN,pSD_ATTN] = corr(Stab{1,2},Disc{1,2},'rows','pairwise','type','Spearman');
    [rSD_WM,pSD_WM] = corr(Stab{3,4},Disc{3,4},'rows','pairwise','type','Spearman');

    [rTO_ATTN,pTO_ATTN] = corr(Typ{1,2},Opt_SA{1,2},'rows','pairwise','type','Spearman');
    [rTO_WM,pTO_WM] = corr(Typ{3,4},Opt_WM{3,4},'rows','pairwise','type','Spearman');

    [rTD_ATTN,pTD_ATTN] = corr(Typ{1,2},Disc{1,2},'rows','pairwise','type','Spearman');
    [rTD_WM,pTD_WM] = corr(Typ{3,4},Disc{3,4},'rows','pairwise','type','Spearman');

    [rOD_ATTN,pOD_ATTN] = corr(Opt_SA{1,2},Disc{1,2},'rows','pairwise','type','Spearman');
    [rOD_WM,pOD_WM] = corr(Opt_WM{3,4},Disc{3,4},'rows','pairwise','type','Spearman');

end % if which_DS == 


% =======================================================================
%% Save data
% =======================================================================

if which_DS == 1
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_StabTypDisc'],'Stab','Typ','Disc','Opt_SA','behav_ATTN','corr_stab_ATTN', 'corr_typ_ATTN', 'corr_disc_ATTN','sig_stab_ATTN', 'sig_typ_ATTN', 'sig_disc_ATTN','r_stab_ds1', 'p_stab_ds1', 'r_typ_ds1', 'p_typ_ds1', 'r_disc_ds1', 'p_disc_ds1', 'r_opt_ds1', 'p_opt_ds1', 'Avg_Stab', 'Avg_Typ', 'Avg_Disc', 'Avg_Opt','OverStab','corr_run1_ATTN','corr_run2_ATTN','corr_run3_ATTN','sig_run1_ATTN','sig_run2_ATTN','sig_run3_ATTN','OverStab');
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_OptFC'],'fc_topSA');
elseif which_DS == 2
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_StabTypDisc'],'corr_run1_ATTN','corr_run2_ATTN','corr_run1_WM','corr_run2_WM','corr_consistencycontrol_ATTN','corr_consistencycontrol_WM','sig_consistencycontrol_ATTN','sig_consistencycontrol_WM','corr_attn_wm','sig_attn_wm','Stab','Typ','Disc','Opt_SA','Opt_WM','behav_ATTN','behav_WM','corr_stab_ATTN', 'corr_typ_ATTN', 'corr_disc_ATTN', 'corr_opt_ATTN', 'sig_stab_ATTN', 'sig_typ_ATTN', 'sig_disc_ATTN', 'sig_opt_ATTN', 'corr_stab_WM', 'corr_typ_WM', 'corr_disc_WM', 'corr_opt_WM', 'sig_stab_WM', 'sig_typ_WM', 'sig_disc_WM', 'sig_opt_WM','OverStab','corr_run1_ATTN','corr_run2_ATTN','corr_run1_WM','corr_run2_WM','sig_run1_ATTN','sig_run2_ATTN','sig_run1_WM','sig_run2_WM');
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_OptFC'],'fc_topSA','fc_topWM');
elseif which_DS == 3
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_StabTypDisc'],'corr_attn_wm','sig_attn_wm','corr_rest_ATTN','sig_rest_ATTN','sig_rest_WM','corr_rest_WM','RestStab','Stab','Typ','Disc','Opt_SA','Opt_WM','behav_ATTN','behav_WM','corr_stab_ATTN', 'corr_typ_ATTN', 'corr_disc_ATTN', 'corr_opt_ATTN', 'sig_stab_ATTN', 'sig_typ_ATTN', 'sig_disc_ATTN', 'sig_opt_ATTN', 'corr_stab_WM', 'corr_typ_WM', 'corr_disc_WM', 'corr_opt_WM', 'sig_stab_WM', 'sig_typ_WM', 'sig_disc_WM', 'sig_opt_WM','OverStab','corr_resttyp_ATTN', 'sig_resttyp_ATTN','corr_resttyp_WM', 'sig_resttyp_WM','RestTyp','OverRest_typ_ATTN', 'sig_OverRest_typ_ATTN' ,'OverRest_typ_WM', 'sig_OverRest_typ_WM');
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_OptFC'],'fc_topSA','fc_topWM');
end


