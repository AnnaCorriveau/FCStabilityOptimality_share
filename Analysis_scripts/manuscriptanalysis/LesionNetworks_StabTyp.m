
% Computationally lesion networks in the brain and calculate correlation btw stability and performance
% Goal is to determine whether canonical networks significantly contribute
% to relationship observed

% AC April 2022

clear
% Which dataset do you want to calculate?
which_DS = [3];


if which_DS == 1
    load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS1_data/gradCPT_flatmats_motion_behav_data']);
elseif which_DS == 2
    load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS2_data/Chun_flatmats_motion_behav']);
    run_names = {'gradCPT','VSTM','MOT','movie','rest'};     
elseif which_DS == 3
    load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS3_data/HCP_flatmats_motion_behav_data']);
end


% load 10-network assignments:
        load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/map268_subnetwork_062521.mat');
        pre_networks(:,1) = map.oldroi;
        pre_networks(:,2) = map.category;

        [~,inds] = sort(pre_networks(:,1));
        networks = pre_networks(inds,:);
        Shen_network_labels = networks(:,2);
        network_names = {'Medial Frontal', 'Frontoparietal', 'Default Mode', 'Motor', 'Visual I', 'Visual II', 'Visual Association','Limbic','Basal Ganglia','Cerebellum'};


blankmat = ones([268]);
trlind = find(tril(blankmat, -1));

    
for net = 1:max(Shen_network_labels)
disp(net);

    num = 1:length(Shen_network_labels);
    which = find(Shen_network_labels == net); % which indices are actually the network(want to make sure we don't include these in MC stats)
   
    num(which) = NaN; % delete those as possible options
    inds = num(~isnan(num));
        
   tmpmat = blankmat;
   tmpmat(which,:) = NaN;
   tmpmat(:,which) = NaN;
   
   lesion_log = tmpmat(trlind) == 1;



for p = 1:size(fc_flat{1},1)

   for r1 = 1:4%size(fc_flat,2) % We are only interested in first 4 runs (SA and WM) so save computational power
         for r2 = 1:4%size(fc_flat,2)


fc_flat{r1}(isinf(fc_flat{r1})) = NaN;
fc_flat{r2}(isinf(fc_flat{r2})) = NaN;
               
          run1 = fc_flat{r1}(p,lesion_log);
          run2 = fc_flat{r2}(p,lesion_log);
          
          self_corr = corr(run1', run2', 'rows','pairwise');
            
          self_FC = tanh(nanmean([run1; run2])');
            

          others1 = fc_flat{r1}(:,lesion_log);;
               others1(p,:) = [];
          others2 = fc_flat{r2}(:,lesion_log);;
               others2(p,:) = [];

           cat_others = cat(1, others1, others2);


           others_FC = tanh(nanmean(cat_others)');

           Typ{r1,r2}(p,1) = corr(self_FC, others_FC,'rows','pairwise');
            Stab{r1,r2}(p,:) = self_corr;
%            


  
end % for r1 = 1:size(all_runs,2)
end % for r2 = 1:size(all_runs,2)

end
            

for r1 = 1:4%size(Stab,2)
         for r2 = 1:4%size(Stab,2)

              stab = Stab{r1,r2};
             typ = Typ{r1,r2};

 
         [corr_stab_ATTN{net}(r1,r2), sig_stab_ATTN{net}(r1,r2)] = partialcorr(stab, behav_ATTN, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');     
        [corr_typ_ATTN{net}(r1,r2), sig_typ_ATTN{net}(r1,r2)] = partialcorr(typ, behav_ATTN, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');        
        

     if which_DS ~= 1

     [corr_stab_WM{net}(r1,r2), sig_stab_WM{net}(r1,r2)] = partialcorr(stab, behav_WM, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');     
     [corr_typ_WM{net}(r1,r2), sig_typ_WM{net}(r1,r2)] = partialcorr(typ, behav_WM, [nanmean([mot{r1}, mot{r2}],2), abs(mot{r1}-mot{r2})],'rows','pairwise','type','Spearman');     

     end       

         end
end

% Since Dataset 1 has 3 gradCPT runs, we will average values across runs
if which_DS == 1
 Avg_Stab = tanh(nanmean([atanh(Stab{1,2}), atanh(Stab{1,3}), atanh(Stab{2,3})],2));
Avg_Typ = tanh(nanmean([atanh(Typ{1,2}), atanh(Typ{1,3}), atanh(Typ{2,3})],2));

 [r_stab_ds1{net},p_stab_ds1{net}] = partialcorr(Avg_Stab, behav_ATTN, [nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{1} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');
[r_typ_ds1{net},p_typ_ds1{net}] = partialcorr(Avg_Typ, behav_ATTN, [nanmean([mot{1},mot{2},mot{3}],2),nanmean([abs(mot{1}-mot{2}), abs(mot{1} - mot{3}), abs(mot{2}-mot{3})],2)],'rows','pairwise','type','Spearman');

end

end % for net = 


if which_DS == 1
    save(['//Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_StabAnat'],'Stab','corr_stab_ATTN','sig_stab_ATTN','r_stab_ds1', 'p_stab_ds1', 'Avg_Stab','r_typ_ds1','Avg_Typ');
elseif which_DS == 2
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_StabAnat'],'Stab','Typ','corr_typ_ATTN','corr_typ_WM','corr_stab_ATTN', 'sig_stab_ATTN', 'corr_stab_WM', 'sig_stab_WM');
elseif which_DS == 3
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_StabAnat'],'Stab','Typ','corr_typ_ATTN','corr_typ_WM','corr_stab_ATTN', 'sig_stab_ATTN', 'corr_stab_WM', 'sig_stab_WM');
end