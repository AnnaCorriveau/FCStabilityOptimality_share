% Calculate the saCPM and wmCPM summary features for each ppt in each
% dataset. Output these values in a data format that can be called in R.

% Rewritten March 2022 AC

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

trlinds = find(tril(ones(size(fc_flat,2)),-1));


% Load in CPM masks
    mat = ones([268]);
    inds = find(tril(mat,-1));

    % load saCPM
    load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/CPM_masks/saCPM.mat')
    hi_mask = high_attention_mask(inds);
    lo_mask = low_attention_mask(inds);

    % load wmCPM
    load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/CPM_masks/HCP_WMmasks_final.mat')
    hi_wm = tril(pos_overlap_wm_task.',-1) + triu(pos_overlap_wm_task);
    lo_wm = tril(neg_overlap_wm_task.',-1) + triu(neg_overlap_wm_task);

    hi_wm_mask = hi_wm(inds);
    lo_wm_mask = lo_wm(inds);





for p = 1:size(fc_flat{1},1)


if which_DS > 1 % If we are working with DS 2 or 3, use {1,2} for SA and {3,4} for WM

   for r1 = 1:4%size(fc_flat,2)
         for r2 = 1:4%size(fc_flat,2)

    
            fc_flat{r1}(isinf(fc_flat{r1})) = NaN;
            fc_flat{r2}(isinf(fc_flat{r2})) = NaN;

            run1 = fc_flat{r1}(p,:);
            run2 = fc_flat{r2}(p,:);
          
            
             self_FC = nanmean([run1; run2])';
            
            % calculate sa and wm summary features
             if r1 == 1 & r2 == 2
              sa_strength(p)=nanmean(self_FC(hi_mask>0))-nanmean(self_FC(lo_mask>0));
             elseif r1 == 3 & r2 == 4
              wm_strength(p)=nanmean(self_FC(hi_wm_mask>0))-nanmean(self_FC(lo_wm_mask>0));
             end


            

  
end % for r1 = 1:size(all_runs,2)
end % for r2 = 1:size(all_runs,2)


            
elseif which_DS == 1

            fc_flat{1}(isinf(fc_flat{1})) = NaN;
            fc_flat{2}(isinf(fc_flat{2})) = NaN;
            fc_flat{3}(isinf(fc_flat{3})) = NaN;

            run1 = fc_flat{1}(p,:);
            run2 = fc_flat{2}(p,:);
            run3 = fc_flat{3}(p,:);
          
            
             self_FC = nanmean([run1; run2; run3])';
            
            % calculate sa and wm summary features
              sa_strength(p)=nanmean(self_FC(hi_mask>0))-nanmean(self_FC(lo_mask>0));


end % if which_DS = 

end % for p = 



if which_DS == 1
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_CPMout'],'sa_strength');
else
    save(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(which_DS) '_CPMout'],'sa_strength','wm_strength');
end


%% Create a table to output into the R GLM code
% load each dataset and combine into a table

clear

ct = 0;

mot_ds{1} = '/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS1_data/gradCPT_flatmats_motion_behav_data';
mot_ds{2} = '/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS2_data/Chun_flatmats_motion_behav';
mot_ds{3} = '/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS3_data/HCP_flatmats_motion_behav_data';


for d = 1:3
    
 
 load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(d) '_StabTypDisc_OptCutoff']);
 new_ct = ct + size(Stab{1,2},1);

 load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(d) '_CPMout']); % CPM model vals

if d == 1

  
 df_out_SA(ct+1:new_ct,1) = normalize(behav_ATTN); % normalized behavior
 df_out_SA(ct+1:new_ct,2) = Avg_Stab; % normalized stability
 df_out_SA(ct+1:new_ct,3) = Avg_Typ; % normalized typicality
 df_out_SA(ct+1:new_ct,4) = Avg_Disc; % normalized discriminability
 df_out_SA(ct+1:new_ct,5) = Avg_Opt; % normalized optimality
 df_out_SA(ct+1:new_ct,6) = repelem(d,length(behav_ATTN));
 df_out_SA(ct+1:new_ct,7) = sa_strength;% normalized saCPM strength summary feature

 load([mot_ds{d}]);
 df_out_SA(ct+1:new_ct,8) = nanmean([mot{1},mot{2},mot{3}],2);
 df_out_SA(ct+1:new_ct,9) = nanmean([abs(mot{1}-mot{2}), abs(mot{1} - mot{3}), abs(mot{2}-mot{3})],2);

else
 
    df_out_SA(ct+1:new_ct,1) = normalize(behav_ATTN); % normalized behavior
  df_out_SA(ct+1:new_ct,2) = Stab{1,2}; % normalized stability
  df_out_SA(ct+1:new_ct,3) = Typ{1,2}; % normalized typicality
  df_out_SA(ct+1:new_ct,4) = Disc{1,2}; % normalized discriminability
  df_out_SA(ct+1:new_ct,5) = Opt_SA{1,2}; % normalized optimality
 df_out_SA(ct+1:new_ct,6) = repelem(d,length(behav_ATTN));
 df_out_SA(ct+1:new_ct,7) = sa_strength;% normalized saCPM strength summary feature

 load([mot_ds{d}]);
 df_out_SA(ct+1:new_ct,8) = nanmean([mot{1}, mot{2}],2);
 df_out_SA(ct+1:new_ct,9) = abs(mot{1}-mot{2});

end

 ct = new_ct;

end


% repeat for WM
ct = 0;

for d = 2:3
    
 
 load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(d) '_StabTypDisc_OptCutoff']);
 load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS' num2str(d) '_CPMout']); % CPM model vals

 new_ct = ct + size(Stab{1,2},1);


 df_out_WM(ct+1:new_ct,1) = normalize(behav_WM); % normalized behavior
 df_out_WM(ct+1:new_ct,2) = Stab{3,4}; % normalized stability
 df_out_WM(ct+1:new_ct,3) = Typ{3,4}; % normalized typicality
 df_out_WM(ct+1:new_ct,4) = Disc{3,4}; % normalized discriminability
 df_out_WM(ct+1:new_ct,5) = Opt_WM{3,4}; % normalized optimality
 df_out_WM(ct+1:new_ct,6) = repelem(d,length(behav_WM));
 df_out_WM(ct+1:new_ct,7) = wm_strength; % normalized wm CPM strength summary feature

 load([mot_ds{d}]);
 df_out_WM(ct+1:new_ct,8) = nanmean([mot{3}, mot{4}],2);
 df_out_WM(ct+1:new_ct,9) = abs(mot{3}-mot{4});

 ct = new_ct;

end









writematrix(df_out_SA,['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/FC_DOTS_GLMinput_SA_NoNorm.txt'],'FileType','text','Delimiter','\t');
writematrix(df_out_WM,['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/FC_DOTS_GLMinput_WM_NoNorm.txt'],'FileType','text','Delimiter','\t');




