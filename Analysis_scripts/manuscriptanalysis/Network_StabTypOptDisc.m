% Calculate stability and typicality within networks to understand which
% networks are more stable and typical than others
% AC March 2022

clear
% Which dataset do you want to calculate?
which_DS = [2];


if which_DS == 1
    load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS1_data/gradCPT_flatmats_motion_behav_data']);
elseif which_DS == 2
    load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS2_data/Chun_flatmats_motion_behav']);
    run_names = {'gradCPT','VSTM','MOT','movie','rest'};     
elseif which_DS == 3
    load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS3_data/HCP_flatmats_motion_behav_data']);
end

top_percent = .05;

%trlinds = find(tril(ones(size(fc_flat,2)),-1));


% load 10-network assignments:
        load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/map268_subnetwork_062521.mat');
        pre_networks(:,1) = map.oldroi;
        pre_networks(:,2) = map.category;

        [~,inds] = sort(pre_networks(:,1));
        networks = pre_networks(inds,:);
        Shen_network_labels = networks(:,2);
        network_names = {'Medial Frontal', 'Frontoparietal', 'Default Mode', 'Motor', 'Visual I', 'Visual II', 'Visual Association','Limbic','Basal Ganglia','Cerebellum'};


blankmat = zeros([268]);

trilmat = ones([268]);
trlind = find(tril(trilmat, -1));

% % identify top performers
% topSA_log = behav_ATTN == max(behav_ATTN);
% 
% if which_DS > 1
% 
% fc_topSA = tanh(nanmean([fc_flat{1}(topSA_log,:); fc_flat{2}(topSA_log,:)]));
% topWM_log = behav_WM == max(behav_WM);
% fc_topWM = tanh(nanmean([fc_flat{3}(topWM_log,:); fc_flat{4}(topWM_log,:)]));
% 
% fc_flat_minustop = fc_flat;
% fc_flat_minustop{1}(topSA_log,:) = NaN;
% fc_flat_minustop{2}(topSA_log,:) = NaN;
% fc_flat_minustop{3}(topWM_log,:) = NaN;
% fc_flat_minustop{4}(topWM_log,:) = NaN;
% 
% else
% fc_topSA = tanh(nanmean([fc_flat{1}(topSA_log,:); fc_flat{2}(topSA_log,:); fc_flat{3}(topSA_log,:)])); % Do not take atanh before averaging, these are already fisher-z scored
% fc_flat_minustop = fc_flat;
% fc_flat_minustop{1}(topSA_log,:) = NaN;
% fc_flat_minustop{2}(topSA_log,:) = NaN;
% fc_flat_minustop{3}(topSA_log,:) = NaN;
% end

% identify top performers
topK = ceil(top_percent * length(behav_ATTN(~isnan(behav_ATTN))));
% top_scores = mink(behav_ATTN,topK);
top_scores = maxk(behav_ATTN,topK);
top_behav_ATTN = zeros([1 length(behav_ATTN)]);
for t = 1:length(top_scores)
where = find(behav_ATTN == top_scores(t));
top_behav_ATTN(where) = 1;
end

% new
top_score_onlySA = behav_ATTN == max(behav_ATTN);
%if is new
if sum(top_score_onlySA) > sum(top_behav_ATTN)
    %topSA_log = behav_ATTN == max(behav_ATTN);
    topSA_log = top_score_onlySA;
    disp(['SA opt N = ' num2str(sum(top_score_onlySA))])
else

%this is the only old part
topSA_log = logical(top_behav_ATTN);
disp(['SA opt N = ' num2str(sum(top_behav_ATTN))])
 end



if which_DS > 1

fc_topSA = tanh(nanmean([fc_flat{1}(topSA_log,:); fc_flat{2}(topSA_log,:)]));

% identify top performers WM
topK_WM = ceil(top_percent * length(behav_WM(~isnan(behav_WM))));
top_scores_WM = maxk(behav_WM,topK_WM);
%top_scores_WM = mink(behav_WM,topK_WM);
top_behav_WM = zeros([1 length(behav_WM)]);
for t = 1:length(top_scores_WM)
where_WM = find(behav_WM == top_scores_WM(t));
top_behav_WM(where_WM) = 1;
end


% new
top_score_onlyWM = behav_WM == max(behav_WM);
%if is new
if sum(top_score_onlyWM) > sum(top_behav_WM)
  
    topWM_log = top_score_onlyWM;
    disp(['WM opt N = ' num2str(sum(top_score_onlyWM))])
else

%this is the only old part
topWM_log = logical(top_behav_WM);
 disp(['WM opt N = ' num2str(sum(top_behav_WM))])

 end


%topWM_log = behav_WM == max(behav_WM);
fc_topWM = tanh(nanmean([fc_flat{3}(topWM_log,:); fc_flat{4}(topWM_log,:)]));

fc_flat_minustop = fc_flat;
fc_flat_minustop{1}(topSA_log,:) = NaN;
fc_flat_minustop{2}(topSA_log,:) = NaN;
fc_flat_minustop{3}(topWM_log,:) = NaN;
fc_flat_minustop{4}(topWM_log,:) = NaN;

else
fc_topSA = tanh(nanmean([fc_flat{1}(topSA_log,:); fc_flat{2}(topSA_log,:); fc_flat{3}(topSA_log,:)])); % Do not take atanh before averaging, these are already fisher-z scored
fc_flat_minustop = fc_flat;
fc_flat_minustop{1}(topSA_log,:) = NaN;
fc_flat_minustop{2}(topSA_log,:) = NaN;
fc_flat_minustop{3}(topSA_log,:) = NaN;
end


    
for net1 = 1:max(Shen_network_labels)
    disp(net1)
    for net2 = 1:max(Shen_network_labels)

    num = 1:length(Shen_network_labels);
    which1 = Shen_network_labels == net1; % which indices are actually the network(want to make sure we don't include these in MC stats)
   which2 = Shen_network_labels == net2;
   
   
    tmpmat = blankmat;
    tmpmat(which1,which2) = 1;
    tmpmat(which2,which1) = 1;
%    
    net_log = tmpmat(trlind) == 1;
% 
% 

for p = 1:size(fc_flat{1},1)

   for r1 = 1:4%size(fc_flat,2) % We are only interested in first 4 runs (SA and WM) so save computational power
         for r2 = 1:4%size(fc_flat,2)


fc_flat{r1}(isinf(fc_flat{r1})) = NaN;
fc_flat{r2}(isinf(fc_flat{r2})) = NaN;

 fc_flat_minustop{r1}(isinf(fc_flat_minustop{r1})) = NaN;
 fc_flat_minustop{r2}(isinf(fc_flat_minustop{r2})) = NaN;
% 
               
          run1 = fc_flat{r1}(p,net_log);
          run2 = fc_flat{r2}(p,net_log);
          
         %run1(isinf(run1)) = NaN;
         %run2(isinf(run2)) = NaN;

          self_corr = corr(run1', run2', 'rows','pairwise');
            
          self_FC = tanh(nanmean([run1; run2])');
            

          others1 = fc_flat{r1}(:,net_log);;
               others1(p,:) = [];
               %others1(isinf(others1)) = NaN;
          others2 = fc_flat{r2}(:,net_log);;
               others2(p,:) = [];
               %others2(isinf(others2)) = NaN;

           cat_others = cat(1, others1, others2);
%            for c = 1:size(cat_others, 1)
%                cat_others(c,:) = atanh(cat_others(c,:));
%            end

           others_FC = tanh(nanmean(cat_others)');

           Typ{net1,net2}{r1,r2}(p,1) = corr(self_FC, others_FC,'rows','pairwise');
            Stab{net1,net2}{r1,r2}(p,:) = self_corr;
%            
%            tmpstab(r1,r2) = self_corr;
%      
           % Optimality analysis
            run1_top = fc_flat_minustop{r1}(p,net_log);
            run2_top = fc_flat_minustop{r2}(p,net_log);
            self_FC_top = tanh(nanmean([run1_top; run2_top])');

            fc_topSA_net = fc_topSA(net_log);
            Opt_SA{net1,net2}{r1,r2}(p,1) = corr(self_FC_top, fc_topSA_net' ,'rows','pairwise');

            if which_DS > 1
                fc_topWM_net = fc_topWM(net_log);
            Opt_WM{net1,net2}{r1,r2}(p,1) = corr(self_FC_top, fc_topWM_net' ,'rows','pairwise');
            end

  
end % for r1 = 1:size(all_runs,2)
end % for r2 = 1:size(all_runs,2)
%    stab_fishz = atanh(tmpstab);
%    OverStab(p,:) = tanh(nanmean(stab_fishz(trlinds)));
end


% 
 for r1 = 1:4%size(fc_flat,2)
         for r2 = 1:4%size(fc_flat,2)
             
             Disc{net1,net2}{r1,r2} = atanh(Stab{net1,net2}{r1,r2})./atanh(Typ{net1,net2}{r1,r2});
             
         end
 end


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

if which_DS == 1
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_NetworkStabTypOptDisc_forReview'],'Net_stab_SA','Net_typ_SA','Net_disc_SA','Net_opt_SA');
elseif which_DS == 2
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_NetworkStabTypOptDisc_forReview'],'Net_stab_SA','Net_typ_SA','Net_stab_WM','Net_typ_WM','Net_disc_SA','Net_opt_SA','Net_disc_WM','Net_opt_WM');
elseif which_DS == 3
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_NetworkStabTypOptDisc_forReview'],'Net_stab_SA','Net_typ_SA','Net_stab_WM','Net_typ_WM','Net_disc_SA','Net_opt_SA','Net_disc_WM','Net_opt_WM');
end





