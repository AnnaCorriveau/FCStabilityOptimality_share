% Create figures 
clear

%% Figure 2.
% ----------------------------------------------------------------------------------------------
% ------------------    SETUP FIGURE 2     -----------------------------------------------------
% ----------------------------------------------------------------------------------------------

% scatter plot of relationship between performance and stability (red) and
% typicality (blue)

% DS1
load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS1_StabTypDisc_OptCutoff']);
fig2_stab_SA_ds1 = Avg_Stab;
fig2_behav_SA_ds1 = behav_ATTN;
fig3_stabmat_SA_ds1 = corr_stab_ATTN;
    fig3_sig_SA_ds1 = sig_stab_ATTN;
for i = 1:size(fig3_stabmat_SA_ds1,1)
    fig3_stabmat_SA_ds1(i,i) = 0;
end
round_ds1corr = round(r_stab_ds1*1000)/1000;
round_ds1p = round(p_stab_ds1*1000)/1000;
fig5_ds1_stab_attn = r_stab_ds1;
fig5_ds1_typ_attn = r_typ_ds1;

figSup1_ds1_attn(1) = r_typ_ds1;
figSup1_ds1_attn(2) = corr_stab_ATTN(4,5);
figSup1_ds1_attn_sig(1) = p_typ_ds1;
figSup1_ds1_attn_sig(2) = sig_stab_ATTN(4,5);

fig2_typ_SA_ds1 = Avg_Typ;
round_ds1corr_typ = round(r_typ_ds1*1000)/1000;
round_ds1p_typ = round(p_typ_ds1*1000)/1000;

fig2_opt_SA_ds1 = Avg_Opt;
round_ds1corr_opt = round(r_opt_ds1*1000)/1000;
round_ds1p_opt = round(p_opt_ds1*1000)/1000;

fig2_disc_SA_ds1 = Avg_Disc;
round_ds1corr_disc = round(r_disc_ds1*1000)/1000;
round_ds1p_disc = round(p_disc_ds1*1000)/1000;

bf_ds1 = 0.05/(sum(sum(tril(ones([size(corr_stab_ATTN,1)]),-1))));

sup2_ds1_attn = quantile(behav_ATTN,[0 0.25 0.5 0.75 1]);
for i = 1:length(sup2_ds1_attn) - 1
    log = behav_ATTN > sup2_ds1_attn(i) & behav_ATTN <  sup2_ds1_attn(i+1);
    ds1_attn_quants(i) = tanh(nanmean(atanh(Typ{1,2}(log))));

    ds1_attn_err(i) = tanh(std(atanh(Typ{1,2}(log))))/sqrt(sum(log));
end




clear behav_ATTN

% DS2
load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS2_StabTypDisc_OptCutoff']);
fig2_stab_SA_ds2 = Stab{1,2};
fig2_behav_SA_ds2 = behav_ATTN;
fig3_stabmat_SA_ds2 = corr_stab_ATTN;
    fig3_sig_SA_ds2 = sig_stab_ATTN;
for i = 1:size(fig3_stabmat_SA_ds2,1)
    fig3_stabmat_SA_ds2(i,i) = 0;
end
round_ds2corr_attn = round(corr_stab_ATTN(1,2)*1000)/1000;
round_ds2p_attn = round(sig_stab_ATTN(1,2)*1000)/1000;


figSup1_ds2_attn(1) = corr_typ_ATTN(1,2);
figSup1_ds2_attn(2) = corr_typ_ATTN(3,4);
figSup1_ds2_attn(3) = corr_typ_ATTN(5,6);
figSup1_ds2_attn(4) = corr_typ_ATTN(7,8);
figSup1_ds2_attn(5) = corr_typ_ATTN(9,10);
figSup1_ds2_wm(1) = corr_typ_WM(1,2);
figSup1_ds2_wm(2) = corr_typ_WM(3,4);
figSup1_ds2_wm(3) = corr_typ_WM(5,6);
figSup1_ds2_wm(4) = corr_typ_WM(7,8);
figSup1_ds2_wm(5) = corr_typ_WM(9,10);

figSup1_ds2_attn_sig(1) = sig_typ_ATTN(1,2);
figSup1_ds2_attn_sig(2) = sig_typ_ATTN(3,4);
figSup1_ds2_attn_sig(3) = sig_typ_ATTN(5,6);
figSup1_ds2_attn_sig(4) = sig_typ_ATTN(7,8);
figSup1_ds2_attn_sig(5) = sig_typ_ATTN(9,10);
figSup1_ds2_wm_sig(1) = sig_typ_WM(1,2);
figSup1_ds2_wm_sig(2) = sig_typ_WM(3,4);
figSup1_ds2_wm_sig(3) = sig_typ_WM(5,6);
figSup1_ds2_wm_sig(4) = sig_typ_WM(7,8);
figSup1_ds2_wm_sig(5) = sig_typ_WM(9,10);


fig2_stab_WM_ds2 = Stab{3,4};
fig2_behav_WM_ds2 = behav_WM;
fig3_stabmat_WM_ds2 = corr_stab_WM;
    fig3_sig_WM_ds2 = sig_stab_WM;
for i = 1:size(fig3_stabmat_WM_ds2,1)
    fig3_stabmat_WM_ds2(i,i) = 0;
end
round_ds2corr_wm = round(corr_stab_WM(3,4)*1000)/1000;
round_ds2p_wm = round(sig_stab_WM(3,4)*1000)/1000;

fig5_ds2_stab_attn = corr_stab_ATTN(1,2);
fig5_ds2_typ_attn = corr_typ_ATTN(1,2);
fig5_ds2_stab_wm = corr_stab_WM(3,4);
fig5_ds2_typ_wm = corr_typ_WM(3,4);


fig2_typ_SA_ds2 = Typ{1,2};
fig2_typ_WM_ds2 = Typ{3,4};
round_ds2corr_attn_typ = round(corr_typ_ATTN(1,2)*1000)/1000;
round_ds2p_attn_typ = round(sig_typ_ATTN(1,2)*1000)/1000;
round_ds2corr_wm_typ = round(corr_typ_WM(3,4)*1000)/1000;
round_ds2p_wm_typ = round(sig_typ_WM(3,4)*1000)/1000;

fig2_opt_SA_ds2 = Opt_SA{1,2};
fig2_opt_WM_ds2 = Opt_WM{3,4};
round_ds2corr_attn_opt = round(corr_opt_ATTN(1,2)*1000)/1000;
round_ds2p_attn_opt = round(sig_opt_ATTN(1,2)*1000)/1000;
round_ds2corr_wm_opt = round(corr_opt_WM(3,4)*1000)/1000;
round_ds2p_wm_opt = round(sig_opt_WM(3,4)*1000)/1000;

fig2_disc_SA_ds2 = Disc{1,2};
fig2_disc_WM_ds2 = Disc{3,4};
round_ds2corr_attn_disc = round(corr_disc_ATTN(1,2)*1000)/1000;
round_ds2p_attn_disc = round(sig_disc_ATTN(1,2)*1000)/1000;
round_ds2corr_wm_disc = round(corr_disc_WM(3,4)*1000)/1000;
round_ds2p_wm_disc = round(sig_disc_WM(3,4)*1000)/1000;

bf_ds2 = 0.05/(sum(sum(tril(ones([size(corr_stab_ATTN,1)]),-1))));

sup2_ds2_attn = quantile(behav_ATTN,[0 0.25 0.5 0.75 1]);
for i = 1:length(sup2_ds2_attn) - 1
    log = behav_ATTN > sup2_ds2_attn(i) & behav_ATTN <  sup2_ds2_attn(i+1);
    ds2_attn_quants(i) = tanh(nanmean(atanh(Typ{1,2}(log))));
    ds2_attn_err(i) = tanh(std(atanh(Typ{1,2}(log))))/sqrt(sum(log));
end

sup2_ds2_wm = quantile(behav_WM,[0 0.25 0.5 0.75 1]);
for i = 1:length(sup2_ds2_wm) - 1
    log = behav_WM > sup2_ds2_wm(i) & behav_WM <  sup2_ds2_wm(i+1);
    ds2_wm_quants(i) = tanh(nanmean(atanh(Typ{3,4}(log))));
    ds2_wm_err(i) = tanh(std(atanh(Typ{3,4}(log))))/sqrt(sum(log));
end




clear Stab behav_ATTN behav_WM Typ Opt_SA Opt_WM Disc

% DS3
load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_StabTypDisc_OptCutoff']);
fig2_stab_SA_ds3 = Stab{1,2};
fig2_behav_SA_ds3 = behav_ATTN;
fig3_stabmat_SA_ds3 = corr_stab_ATTN;
    fig3_sig_SA_ds3 = sig_stab_ATTN;
for i = 1:size(fig3_stabmat_SA_ds3,1)
    fig3_stabmat_SA_ds3(i,i) = 0;
end
round_ds3corr_attn = round(corr_stab_ATTN(1,2)*1000)/1000;
round_ds3p_attn = round(sig_stab_ATTN(1,2)*1000)/1000;

fig2_stab_WM_ds3 = Stab{3,4};
fig2_behav_WM_ds3 = behav_WM;
fig3_stabmat_WM_ds3 = corr_stab_WM;
    fig3_sig_WM_ds3 = sig_stab_WM;
for i = 1:size(fig3_stabmat_WM_ds3,1)
    fig3_stabmat_WM_ds3(i,i) = 0;
end
round_ds3corr_wm = round(corr_stab_WM(3,4)*1000)/1000;
round_ds3p_wm = round(sig_stab_WM(3,4)*1000)/1000;

fig5_ds3_stab_attn = corr_stab_ATTN(1,2);
fig5_ds3_typ_attn = corr_typ_ATTN(1,2);
fig5_ds3_stab_wm = corr_stab_WM(3,4);
fig5_ds3_typ_wm = corr_typ_WM(3,4);

fig2_typ_SA_ds3 = Typ{1,2};
fig2_typ_WM_ds3 = Typ{3,4};
round_ds3corr_attn_typ = round(corr_typ_ATTN(1,2)*1000)/1000;
round_ds3p_attn_typ = round(sig_typ_ATTN(1,2)*1000)/1000;
round_ds3corr_wm_typ = round(corr_typ_WM(3,4)*1000)/1000;
round_ds3p_wm_typ = round(sig_typ_WM(3,4)*1000)/1000;

fig2_opt_SA_ds3 = Opt_SA{1,2};
fig2_opt_WM_ds3 = Opt_WM{3,4};
round_ds3corr_attn_opt = round(corr_opt_ATTN(1,2)*1000)/1000;
round_ds3p_attn_opt = round(sig_opt_ATTN(1,2)*1000)/1000;
round_ds3corr_wm_opt = round(corr_opt_WM(3,4)*1000)/1000;
round_ds3p_wm_opt = round(sig_opt_WM(3,4)*1000)/1000;

fig2_disc_SA_ds3 = Disc{1,2};
fig2_disc_WM_ds3 = Disc{3,4};
round_ds3corr_attn_disc = round(corr_disc_ATTN(1,2)*1000)/1000;
round_ds3p_attn_disc = round(sig_disc_ATTN(1,2)*1000)/1000;
round_ds3corr_wm_disc = round(corr_disc_WM(3,4)*1000)/1000;
round_ds3p_wm_disc = round(sig_disc_WM(3,4)*1000)/1000;

bf_ds3 = 0.05/(sum(sum(tril(ones([size(corr_stab_ATTN,1)]),-1))));

figSup1_ds3_attn(1) = corr_typ_ATTN(1,2);
figSup1_ds3_attn(2) = corr_typ_ATTN(3,4);
figSup1_ds3_attn(3) = corr_typ_ATTN(5,6);
figSup1_ds3_attn(4) = corr_typ_ATTN(7,8);
figSup1_ds3_attn(5) = corr_typ_ATTN(9,10);
figSup1_ds3_attn(6) = corr_typ_ATTN(1,2);
figSup1_ds3_attn(7) = corr_typ_ATTN(3,4);
figSup1_ds3_attn(8) = corr_typ_ATTN(5,6);
load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_RestTyps.mat']);
figSup1_ds3_attn(9) = OverRest_typ_ATTN;

figSup1_ds3_wm(1) = corr_typ_WM(1,2);
figSup1_ds3_wm(2) = corr_typ_WM(3,4);
figSup1_ds3_wm(3) = corr_typ_WM(5,6);
figSup1_ds3_wm(4) = corr_typ_WM(7,8);
figSup1_ds3_wm(5) = corr_typ_WM(9,10);
figSup1_ds3_wm(6) = corr_typ_WM(11,12);
figSup1_ds3_wm(7) = corr_typ_WM(13,14);
figSup1_ds3_wm(8) = corr_typ_WM(15,16);
figSup1_ds3_wm(9) = OverRest_typ_WM;

figSup1_ds3_attn_sig(1) = sig_typ_ATTN(1,2);
figSup1_ds3_attn_sig(2) = sig_typ_ATTN(3,4);
figSup1_ds3_attn_sig(3) = sig_typ_ATTN(5,6);
figSup1_ds3_attn_sig(4) = sig_typ_ATTN(7,8);
figSup1_ds3_attn_sig(5) = sig_typ_ATTN(9,10);
figSup1_ds3_attn_sig(6) = sig_typ_ATTN(1,2);
figSup1_ds3_attn_sig(7) = sig_typ_ATTN(3,4);
figSup1_ds3_attn_sig(8) = sig_typ_ATTN(5,6);
figSup1_ds3_attn_sig(9) = sig_OverRest_typ_ATTN;

figSup1_ds3_wm_sig(1) = sig_typ_WM(1,2);
figSup1_ds3_wm_sig(2) = sig_typ_WM(3,4);
figSup1_ds3_wm_sig(3) = sig_typ_WM(5,6);
figSup1_ds3_wm_sig(4) = sig_typ_WM(7,8);
figSup1_ds3_wm_sig(5) = sig_typ_WM(9,10);
figSup1_ds3_wm_sig(6) = sig_typ_WM(1,2);
figSup1_ds3_wm_sig(7) = sig_typ_WM(3,4);
figSup1_ds3_wm_sig(8) = sig_typ_WM(5,6);
figSup1_ds3_wm_sig(9) = sig_OverRest_typ_WM;

sup2_ds3_attn = quantile(behav_ATTN,[0 0.25 0.5 0.75 1]);
for i = 1:length(sup2_ds3_attn) - 1
    log = behav_ATTN > sup2_ds3_attn(i) & behav_ATTN <  sup2_ds3_attn(i+1);
    ds3_attn_quants(i) = tanh(nanmean(atanh(Typ{1,2}(log))));
    ds3_attn_err(i) = tanh(std(atanh(Typ{1,2}(log))))/sqrt(sum(log));
end

sup2_ds3_wm = quantile(behav_WM,[0 0.25 0.5 0.75 1]);
for i = 1:length(sup2_ds3_wm) - 1
    log = behav_WM > sup2_ds3_wm(i) & behav_WM <  sup2_ds3_wm(i+1);
    ds3_wm_quants(i) = tanh(nanmean(atanh(Typ{3,4}(log))));
    ds3_wm_err(i) = tanh(std(atanh(Typ{3,4}(log))))/sqrt(sum(log));
end

% ----------------------------------------------------------------------------------------------
% ------------------    PLOT FIGURE 2     ------------------------------------------------------
% ----------------------------------------------------------------------------------------------

xpos1 = [0.19, 0.3];
xpos2  = [0.64, 0.3];
ypos1 = [0.96, 0.05];
ypos2 = [0.65, 0.05];
ypos3 = [0.34, 0.05];


% PLOT
figure
% DS1 ATTN
subplot(3,2,1)
scatter(fig2_stab_SA_ds1, fig2_behav_SA_ds1,100,[0.4 0 0],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_stab_SA_ds1) | isnan(fig2_behav_SA_ds1);
coefficients = polyfit(fig2_stab_SA_ds1(~ind), fig2_behav_SA_ds1(~ind), 1);
xFit = linspace(min(fig2_stab_SA_ds1), max(fig2_stab_SA_ds1), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([0 4.5]);
% xlabel('Stability')
ylabel("gradCPT d'")
set(gca,'FontSize',20)
if round_ds1p >= .001
annotation('textbox', [xpos1(1), ypos1(1), xpos1(2), ypos1(2)], 'string', {['partial r_s = ' num2str(round_ds1corr)];['p = ' num2str(round_ds1p)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos1(1), ypos1(1), xpos1(2), ypos1(2)], 'string', {['partial r_s = ' num2str(round_ds1corr)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end


% DS2 ATTN STABILITY
subplot(3,2,3)
scatter(fig2_stab_SA_ds2, fig2_behav_SA_ds2,100,[0.4 0 0],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_stab_SA_ds2) | isnan(fig2_behav_SA_ds2);
coefficients = polyfit(fig2_stab_SA_ds2(~ind), fig2_behav_SA_ds2(~ind), 1);
xFit = linspace(min(fig2_stab_SA_ds2), max(fig2_stab_SA_ds2), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([0 4.5]);
% xlabel('Stability')
ylabel("gradCPT d'")
set(gca,'FontSize',20)
if round_ds2p_attn >= .001
annotation('textbox', [xpos1(1), ypos2(1), xpos1(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_attn)];['p = ' num2str(round_ds2p_attn)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos1(1), ypos2(1), xpos1(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_attn)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end


% DS2 WM STABILITY
subplot(3,2,4)
scatter(fig2_stab_WM_ds2, fig2_behav_WM_ds2,100,[0.4 0 0],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_stab_WM_ds2) | isnan(fig2_behav_WM_ds2);
coefficients = polyfit(fig2_stab_WM_ds2(~ind), fig2_behav_WM_ds2(~ind), 1);
xFit = linspace(min(fig2_stab_WM_ds2), max(fig2_stab_WM_ds2), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([0 4.5]);
% xlabel('Stability')
ylabel("Pashler's K")
set(gca,'FontSize',20)
if round_ds2p_wm >= .001
annotation('textbox', [xpos2(1), ypos2(1), xpos2(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_wm)];['p = ' num2str(round_ds2p_wm)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos2(1), ypos2(1), xpos2(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_wm)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end


% DS3 ATTN STABILITY
subplot(3,2,5)
scatter(fig2_stab_SA_ds3, fig2_behav_SA_ds3,100,[0.4 0 0],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_stab_SA_ds3) | isnan(fig2_behav_SA_ds3);
coefficients = polyfit(fig2_stab_SA_ds3(~ind), fig2_behav_SA_ds3(~ind), 1);
xFit = linspace(min(fig2_stab_SA_ds3), max(fig2_stab_SA_ds3), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([50 100]);
xlabel('Stability')
ylabel('0-back Accuracy')
set(gca,'FontSize',20)
if round_ds3p_attn >= .001
    annotation('textbox', [xpos1(1), ypos3(1)  xpos1(2), ypos3(2)], 'string', {['partial r_s = ' num2str(round_ds3corr_attn)];['p = ' num2str(round_ds3p_attn)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos1(1), ypos3(1)  xpos1(2), ypos3(2)],  'string', {['partial r_s = ' num2str(round_ds3corr_attn)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end
% hold off

% DS3 WM STABILITY
subplot(3,2,6)
scatter(fig2_stab_WM_ds3, fig2_behav_WM_ds3,100,[0.4 0 0],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_stab_WM_ds3) | isnan(fig2_behav_WM_ds3);
coefficients = polyfit(fig2_stab_WM_ds3(~ind), fig2_behav_WM_ds3(~ind), 1);
xFit = linspace(min(fig2_stab_WM_ds3), max(fig2_stab_WM_ds3), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([50 100]);
xlabel('Stability')
ylabel('2-back Accuracy')
set(gca,'FontSize',20)
if round_ds3p_wm >= .001
annotation('textbox', [xpos2(1), ypos3(1)  xpos2(2), ypos3(2)],  'string', {['partial r_s = ' num2str(round_ds3corr_wm)];['p = ' num2str(round_ds3p_wm)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos2(1), ypos3(1)  xpos2(2), ypos3(2)],  'string', {['partial r_s = ' num2str(round_ds3corr_wm)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end



% Label title of figure
annotation('textbox', [0.20, 0.95, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize', 20,'EdgeColor','none')
annotation('textbox', [0.67, 0.95, 0.45, 0.05], 'string', 'B. Working Memory','FontSize', 20,'EdgeColor','none')

set(gcf, 'Position',  [10, 10, 600, 750])

% Label datasets
annotation('textarrow', [0.04, 0.3],[0.89, 0.3], 'String', {['Dataset 1']},'FontSize',25, 'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)
annotation('textarrow', [0.04, 0.3], [0.57, 0.3], 'String', {['Dataset 2']},'FontSize',25,'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)
annotation('textarrow', [0.03, 0.3], [0.25, 0.3], 'String', {['Dataset 3']},'FontSize',25,'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)


% Figure 2.2 Typicality scatters 
figure
%nexttile
% DS1 attn typicality
subplot(3,2,1)
scatter(fig2_typ_SA_ds1, fig2_behav_SA_ds1,100,[0 0.1 0.5],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_typ_SA_ds1) | isnan(fig2_behav_SA_ds1);
coefficients = polyfit(fig2_typ_SA_ds1(~ind), fig2_behav_SA_ds1(~ind), 1);
xFit = linspace(min(fig2_typ_SA_ds1), max(fig2_typ_SA_ds1), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([0 4.5]);
% xlabel('Typicality')
ylabel("gradCPT d'")
set(gca,'FontSize',20)
set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds1p_typ >= .001
annotation('textbox', [xpos1(1), ypos1(1), xpos1(2), ypos1(2)], 'string', {['partial r_s = ' num2str(round_ds1corr_typ)];['p = ' num2str(round_ds1p_typ)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos1(1), ypos1(1), xpos1(2), ypos1(2)], 'string', {['partial r_s = ' num2str(round_ds1corr_typ)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end

% DS2 ATTN TYPICALITY
subplot(3,2,3)
scatter(fig2_typ_SA_ds2, fig2_behav_SA_ds2,100,[0 0.1 0.5],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_typ_SA_ds2) | isnan(fig2_behav_SA_ds2);
coefficients = polyfit(fig2_typ_SA_ds2(~ind), fig2_behav_SA_ds2(~ind), 1);
xFit = linspace(min(fig2_typ_SA_ds2), max(fig2_typ_SA_ds2), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([0 4.5]);
% xlabel('Typicality')
ylabel("gradCPT d'")
set(gca,'FontSize',20)
% set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds2p_attn_typ >= .001
annotation('textbox', [xpos1(1), ypos2(1), xpos1(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_attn_typ)];['p = ' num2str(round_ds2p_attn_typ)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos1(1), ypos2(1), xpos1(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_attn_typ)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end
  
 % DS2 WM TYPICALITY
subplot(3,2,4)
scatter(fig2_typ_WM_ds2, fig2_behav_WM_ds2,100,[0 0.1 0.5],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_typ_WM_ds2) | isnan(fig2_behav_WM_ds2);
coefficients = polyfit(fig2_typ_WM_ds2(~ind), fig2_behav_WM_ds2(~ind), 1);
xFit = linspace(min(fig2_typ_WM_ds2), max(fig2_typ_WM_ds2), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([0 4.5]);
% xlabel('Typicality')
 ylabel("Pashler's K")
set(gca,'FontSize',20)
% set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds2p_wm_typ >= .001
annotation('textbox', [xpos2(1), ypos2(1), xpos2(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_wm_typ)];['p = ' num2str(round_ds2p_wm_typ)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos2(1), ypos2(1), xpos2(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_wm_typ)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end

% DS3 ATTN TYPICALITY
subplot(3,2,5)
scatter(fig2_typ_SA_ds3, fig2_behav_SA_ds3,100,[0 0.1 0.5],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_typ_SA_ds3) | isnan(fig2_behav_SA_ds3);
coefficients = polyfit(fig2_typ_SA_ds3(~ind), fig2_behav_SA_ds3(~ind), 1);
xFit = linspace(min(fig2_typ_SA_ds3), max(fig2_typ_SA_ds3), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([50 100]);
xlabel('Typicality')
 ylabel('0-back Accuracy')
set(gca,'FontSize',20)
% set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds3p_attn_typ >= .001
    annotation('textbox', [xpos1(1), ypos3(1), xpos1(2), ypos3(2)], 'string', {['partial r_s = ' num2str(round_ds3corr_attn_typ)];['p = ' num2str(round_ds3p_attn_typ)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos1(1), ypos3(1), xpos1(2), ypos3(2)], 'string', {['partial r_s = ' num2str(round_ds3corr_attn_typ)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end

% DS3 WM TYPICALITY
subplot(3,2,6)
scatter(fig2_typ_WM_ds3, fig2_behav_WM_ds3,100,[0 0.1 0.5],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_typ_WM_ds3) | isnan(fig2_behav_WM_ds3);
coefficients = polyfit(fig2_typ_WM_ds3(~ind), fig2_behav_WM_ds3(~ind), 1);
xFit = linspace(min(fig2_typ_WM_ds3), max(fig2_typ_WM_ds3), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([50 100]);
xlabel('Typicality')
ylabel('2-back Accuracy')
set(gca,'FontSize',20)
% set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds3p_wm_typ >= .001
annotation('textbox', [xpos2(1), ypos3(1), xpos2(2), ypos3(2)], 'string', {['partial r_s = ' num2str(round_ds3corr_wm_typ)];['p = ' num2str(round_ds3p_wm_typ)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos2(1), ypos3(1), xpos2(2), ypos3(2)], 'string', {['partial r_s = ' num2str(round_ds3corr_wm_typ)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end

% Label title of figure
annotation('textbox', [0.20, 0.95, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize', 20,'EdgeColor','none')
annotation('textbox', [0.67, 0.95, 0.45, 0.05], 'string', 'B. Working Memory','FontSize', 20,'EdgeColor','none')

set(gcf, 'Position',  [10, 10, 600, 750])


% Optimality scatters
% Figure 2.3 Optimality scatters 
figure

% DS1 attn Optimality
subplot(3,2,1)
scatter(fig2_opt_SA_ds1, fig2_behav_SA_ds1,100,[0.3 0.6 0.8],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_opt_SA_ds1) | isnan(fig2_behav_SA_ds1);
coefficients = polyfit(fig2_opt_SA_ds1(~ind), fig2_behav_SA_ds1(~ind), 1);
xFit = linspace(min(fig2_opt_SA_ds1), max(fig2_opt_SA_ds1), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([0 4.5]);
% xlabel('Optimality')
ylabel("gradCPT d'")
set(gca,'FontSize',20)
set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds1p_opt >= .001
annotation('textbox', [xpos1(1), ypos1(1), xpos1(2), ypos1(2)], 'string', {['partial r_s = ' num2str(round_ds1corr_opt)];['p = ' num2str(round_ds1p_opt)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos1(1), ypos1(1), xpos1(2), ypos1(2)], 'string', {['partial r_s = ' num2str(round_ds1corr_opt)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end

% DS2 ATTN optimality
subplot(3,2,3)
scatter(fig2_opt_SA_ds2, fig2_behav_SA_ds2,100,[0.3 0.6 0.8],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_opt_SA_ds2) | isnan(fig2_behav_SA_ds2);
coefficients = polyfit(fig2_opt_SA_ds2(~ind), fig2_behav_SA_ds2(~ind), 1);
xFit = linspace(min(fig2_opt_SA_ds2), max(fig2_opt_SA_ds2), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([0 4.5]);
% xlabel('Optimality')
 ylabel("gradCPT d'")
set(gca,'FontSize',20)
% set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds2p_attn_opt >= .001
annotation('textbox', [xpos1(1), ypos2(1), xpos1(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_attn_opt)];['p = ' num2str(round_ds2p_attn_opt)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos1(1), ypos2(1), xpos1(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_attn_opt)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end
  
 % DS2 WM Optimality
subplot(3,2,4)
scatter(fig2_opt_WM_ds2, fig2_behav_WM_ds2,100,[0.3 0.6 0.8],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_opt_WM_ds2) | isnan(fig2_behav_WM_ds2);
coefficients = polyfit(fig2_opt_WM_ds2(~ind), fig2_behav_WM_ds2(~ind), 1);
xFit = linspace(min(fig2_opt_WM_ds2), max(fig2_opt_WM_ds2), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([0 4.5]);
% xlabel('Optimality')
ylabel("Pashler's K")
set(gca,'FontSize',20)
% set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds2p_wm_opt >= .001
annotation('textbox', [xpos2(1), ypos2(1), xpos2(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_wm_opt)];['p = ' num2str(round_ds2p_wm_opt)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos2(1), ypos2(1), xpos2(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_wm_opt)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end

% DS3 ATTN Optimality
subplot(3,2,5)
scatter(fig2_opt_SA_ds3, fig2_behav_SA_ds3,100,[0.3 0.6 0.8],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_opt_SA_ds3) | isnan(fig2_behav_SA_ds3);
coefficients = polyfit(fig2_opt_SA_ds3(~ind), fig2_behav_SA_ds3(~ind), 1);
xFit = linspace(min(fig2_opt_SA_ds3), max(fig2_opt_SA_ds3), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([50 100]);
xlabel('Optimality')
ylabel('0-back Accuracy')
set(gca,'FontSize',20)
% set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds3p_attn_opt >= .001
    annotation('textbox', [xpos1(1), ypos3(1), xpos1(2), ypos3(2)], 'string', {['partial r_s = ' num2str(round_ds3corr_attn_opt)];['p = ' num2str(round_ds3p_attn_opt)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos1(1), ypos3(1), xpos1(2), ypos3(2)], 'string', {['partial r_s = ' num2str(round_ds3corr_attn_opt)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end




% DS3 WM Optimality
subplot(3,2,6)
scatter(fig2_opt_WM_ds3, fig2_behav_WM_ds3,100,[0.3 0.6 0.8],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_opt_WM_ds3) | isnan(fig2_behav_WM_ds3);
coefficients = polyfit(fig2_opt_WM_ds3(~ind), fig2_behav_WM_ds3(~ind), 1);
xFit = linspace(min(fig2_opt_WM_ds3), max(fig2_opt_WM_ds3), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 1]);
ylim([50 100]);
xlabel('Optimality')
 ylabel('2-back Accuracy')
set(gca,'FontSize',20)
% set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds3p_wm_opt >= .001
annotation('textbox', [xpos2(1), ypos3(1), xpos2(2), ypos3(2)], 'string', {['partial r_s = ' num2str(round_ds3corr_wm_opt)];['p = ' num2str(round_ds3p_wm_opt)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos2(1), ypos3(1), xpos2(2), ypos3(2)], 'string', {['partial r_s = ' num2str(round_ds3corr_wm_opt)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end

% Label title of figure
annotation('textbox', [0.20, 0.95, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize', 20,'EdgeColor','none')
annotation('textbox', [0.67, 0.95, 0.45, 0.05], 'string', 'B. Working Memory','FontSize', 20,'EdgeColor','none')

set(gcf, 'Position',  [10, 10, 600, 750])





% Discriminabilty

% Figure 2.3 Discriminability scatters 
figure

% DS1 attn Discriminability
subplot(3,2,1)
scatter(fig2_disc_SA_ds1, fig2_behav_SA_ds1,100,[0 0.4 0.0],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_disc_SA_ds1) | isnan(fig2_behav_SA_ds1);
coefficients = polyfit(fig2_disc_SA_ds1(~ind), fig2_behav_SA_ds1(~ind), 1);
xFit = linspace(min(fig2_disc_SA_ds1), max(fig2_disc_SA_ds1), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 4]);
ylim([0 4.5]);
% xlabel('Discriminability')
 ylabel("gradCPT d'")
set(gca,'FontSize',20)
set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds1p_disc >= .001
annotation('textbox', [xpos1(1), ypos1(1), xpos1(2), ypos1(2)], 'string', {['partial r_s = ' num2str(round_ds1corr_disc)];['p = ' num2str(round_ds1p_disc)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos1(1), ypos1(1), xpos1(2), ypos1(2)], 'string', {['partial r_s = ' num2str(round_ds1corr_disc)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end

% DS2 ATTN Discriminability
subplot(3,2,3)
scatter(fig2_disc_SA_ds2, fig2_behav_SA_ds2,100,[0 0.4 0.0],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_disc_SA_ds2) | isnan(fig2_behav_SA_ds2);
coefficients = polyfit(fig2_disc_SA_ds2(~ind), fig2_behav_SA_ds2(~ind), 1);
xFit = linspace(min(fig2_disc_SA_ds2), max(fig2_disc_SA_ds2), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 4]);
ylim([0 4.5]);
% xlabel('Discriminability')
ylabel("gradCPT d'")
set(gca,'FontSize',20)
% set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds2p_attn_disc >= .001
annotation('textbox', [xpos1(1), ypos2(1), xpos1(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_attn_disc)];['p = ' num2str(round_ds2p_attn_disc)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos1(1), ypos2(1), xpos1(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_attn_disc)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end
  
 % DS2 WM Discriminability
subplot(3,2,4)
scatter(fig2_disc_WM_ds2, fig2_behav_WM_ds2,100,[0 0.4 0.0],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_disc_WM_ds2) | isnan(fig2_behav_WM_ds2);
coefficients = polyfit(fig2_disc_WM_ds2(~ind), fig2_behav_WM_ds2(~ind), 1);
xFit = linspace(min(fig2_disc_WM_ds2), max(fig2_disc_WM_ds2), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 4]);
ylim([0 4.5]);
% xlabel('Discriminability')
ylabel("Pashler's K")
set(gca,'FontSize',20)
% set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds2p_wm_disc >= .001
annotation('textbox', [xpos2(1), ypos2(1), xpos2(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_wm_disc)];['p = ' num2str(round_ds2p_wm_disc)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos2(1), ypos2(1), xpos2(2), ypos2(2)], 'string', {['partial r_s = ' num2str(round_ds2corr_wm_disc)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end

% DS3 ATTN Discriminability
subplot(3,2,5)
scatter(fig2_disc_SA_ds3, fig2_behav_SA_ds3,100,[0 0.4 0.0],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_disc_SA_ds3) | isnan(fig2_behav_SA_ds3);
coefficients = polyfit(fig2_disc_SA_ds3(~ind), fig2_behav_SA_ds3(~ind), 1);
xFit = linspace(min(fig2_disc_SA_ds3), max(fig2_disc_SA_ds3), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 4]);
ylim([50 100]);
xlabel('Discriminability')
ylabel('0-back Accuracy')
set(gca,'FontSize',20)
% set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds3p_attn_disc >= .001
    annotation('textbox', [xpos1(1), ypos3(1), xpos1(2), ypos3(2)], 'string', {['partial r_s = ' num2str(round_ds3corr_attn_disc)];['p = ' num2str(round_ds3p_attn_disc)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos1(1), ypos3(1), xpos1(2), ypos3(2)], 'string', {['partial r_s = ' num2str(round_ds3corr_attn_disc)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end

% DS3 WM Discriminability
subplot(3,2,6)
scatter(fig2_disc_WM_ds3, fig2_behav_WM_ds3,100,[0 0.4 0.0],'filled','MarkerFaceAlpha',.5);hold on
ind = isnan(fig2_disc_WM_ds3) | isnan(fig2_behav_WM_ds3);
coefficients = polyfit(fig2_disc_WM_ds3(~ind), fig2_behav_WM_ds3(~ind), 1);
xFit = linspace(min(fig2_disc_WM_ds3), max(fig2_disc_WM_ds3), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2);
xlim([0 4]);
ylim([50 100]);
xlabel('Discriminability')
 ylabel('2-back Accuracy')
set(gca,'FontSize',20)
% set(gcf, 'Position',  [100, 100, 500, 400])
% axis square
if round_ds3p_wm_disc >= .001
annotation('textbox', [xpos2(1), ypos3(1), xpos2(2), ypos3(2)], 'string', {['partial r_s = ' num2str(round_ds3corr_wm_disc)];['p = ' num2str(round_ds3p_wm_disc)]},'FontSize', 20,'EdgeColor','none')
else
    annotation('textbox', [xpos2(1), ypos3(1), xpos2(2), ypos3(2)], 'string', {['partial r_s = ' num2str(round_ds3corr_wm_disc)];['p < .001']},'FontSize', 20,'EdgeColor','none')
end

% Label title of figure
annotation('textbox', [0.20, 0.95, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize', 20,'EdgeColor','none')
annotation('textbox', [0.67, 0.95, 0.45, 0.05], 'string', 'B. Working Memory','FontSize', 20,'EdgeColor','none')

set(gcf, 'Position',  [10, 10, 600, 750])


%% Create figure 3 -- Stability Matrices
% ----------------------------------------------------------------------------------------------
% ------------------    SETUP FIGURE 3    -----------------------------------------------------
% ----------------------------------------------------------------------------------------------

gradCPT_run_names = {'gradCPT','gradCPT','gradCPT','rest','rest'};
Chun_run_names = {'gradCPT','VSTM','MOT','movie','rest'};
HCP_run_names = {'0-back','2-back','relational','motor','language','social','emotion','gambling','rest'};


% find max/min for color bar
mx(1) = max(max(fig3_stabmat_SA_ds1));
mx(2) = max(max(fig3_stabmat_SA_ds2));
mx(3) = max(max(fig3_stabmat_WM_ds2));
mx(4) = max(max(fig3_stabmat_SA_ds3));
mx(5) = max(max(fig3_stabmat_WM_ds3));

mn(1) = min(min(fig3_stabmat_SA_ds1));
mn(2) = min(min(fig3_stabmat_SA_ds2));
mn(3) = min(min(fig3_stabmat_WM_ds2));
mn(4) = min(min(fig3_stabmat_SA_ds2));
mn(5) = min(min(fig3_stabmat_WM_ds3));

c_hi = abs(max(mx));
c_lo = abs(min(mn));
c_define = max([c_hi c_lo]);
if c_define == c_hi
    clim = [-1*c_hi c_hi];
elseif c_define == c_lo
    clim = [-1*c_lo c_lo];
end


% Find which values are significant 
% off-diagonal bonferroni corrected
% on-diagonal non-corrected
runs1 = [1 1 1 2 2];
for i = 1:size(fig3_stabmat_SA_ds1,1)
for j = 1:size(fig3_stabmat_SA_ds1,1)
if fig3_sig_SA_ds1(i,j) < bf_ds1;
    stars1(i,j) = 1;
else
    stars1(i,j) = NaN;
end % if < bf_ds1
if runs1(i) == runs1(j)
    if fig3_sig_SA_ds1(i,j) < 0.05
        stars1(i,j) = 1;
    end
end % if i == j
end %for j
end % for i
[SS1x_ATTN, SS1y_ATTN] = find(triu(stars1 == 1));


runs2 = [1 1 2 2 3 3 4 4 5 5];
for i = 1:size(fig3_stabmat_SA_ds2,1)
for j = 1:size(fig3_stabmat_SA_ds2,1)
if fig3_sig_SA_ds2(i,j) < bf_ds2;
    stars2_SA(i,j) = 1;
else
    stars2_SA(i,j) = NaN;
end % if < bf_ds1
if fig3_sig_WM_ds2(i,j) < bf_ds2;
    stars2_WM(i,j) = 1;
else
    stars2_SA(i,j) = NaN;
end % if < bf_ds1
if runs2(i) == runs2(j)
    if fig3_sig_SA_ds2(i,j) < 0.05
        stars2_SA(i,j) = 1;
    end
    if fig3_sig_WM_ds2(i,j) < 0.05
        stars2_WM(i,j) = 1;
    end
end % if i == j
end %for j
end % for i
[SS2x_ATTN, SS2y_ATTN] = find(triu(stars2_SA == 1));
[SS2x_WM, SS2y_WM] = find(triu(stars2_WM == 1));


runs3 = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 9 9];
for i = 1:size(fig3_stabmat_SA_ds3,1)
for j = 1:size(fig3_stabmat_SA_ds3,1)
if fig3_sig_SA_ds3(i,j) < bf_ds3;
    stars3_SA(i,j) = 1;
else
    stars3_SA(i,j) = NaN;
end % if < bf_ds3
if fig3_sig_WM_ds3(i,j) < bf_ds3;
    stars3_WM(i,j) = 1;
else
    stars3_SA(i,j) = NaN;
end % if < bf_ds3
if runs3(i) == runs3(j)
    if fig3_sig_SA_ds3(i,j) < 0.05
        stars3_SA(i,j) = 1;
    end
    if fig3_sig_WM_ds3(i,j) < 0.05
        stars3_WM(i,j) = 1;
    end
end % if i == j
end %for j
end % for i
[SS3x_ATTN, SS3y_ATTN] = find(triu(stars3_SA == 1));
[SS3x_WM, SS3y_WM] = find(triu(stars3_WM == 1));







figure;
% DS1 attention stability matrix
subplot(3,2,1);
imagesc(tril(fig3_stabmat_SA_ds1));hold on
plot(SS1x_ATTN, SS1y_ATTN,'w*','MarkerSize',15,'LineWidth',1.5);
colormap(redwhiteblack);
caxis(clim)
hh = colorbar;
ylabel(hh,"Spearman's rho",'FontSize',16,'Rotation',270);
xticks([1:5]);
xticklabels(gradCPT_run_names);
yticks([1:5]);
yticklabels(gradCPT_run_names);
xtickangle(45);
axis square
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',15,'FontWeight','bold','box','off')
% Heavy black lines
plot([0.0 1.5],[0.5 0.5], 'k', 'LineWidth', 8);
plot([1.5 2.5],[1.5 1.5], 'k', 'LineWidth', 8);
plot([2.5 3.5],[2.5 2.5], 'k', 'LineWidth', 8);
plot([0.5 0.5], [0.0 3.5],'k', 'LineWidth', 8); %keep 
plot([0.5 0.5], [0.0 1.5],'k', 'LineWidth', 8);
plot([1.5 1.5], [0.5 1.5],'k', 'LineWidth', 8);
plot([2.5 2.5], [1.5 2.5],'k', 'LineWidth', 8);
plot([0 3.5],[3.5 3.5], 'k', 'LineWidth', 8);% keep
plot([3.5 3.5], [2.5 3.5],'k', 'LineWidth', 8);
% Light gray lines
plot([3.5 4.5], [3.5 3.5], 'Color',[.35 .35 .35], 'LineWidth', 4);
plot([3.5 3.5], [3.5 5.5],'Color',[.35 .35 .35], 'LineWidth', 4);% keep
plot([4.5 5.5], [4.5 4.5], 'Color',[.35 .35 .35], 'LineWidth', 4);
plot([5.5 5.5], [4.5 5.5],'Color',[.35 .35 .35], 'LineWidth', 4);
plot([4.5 4.5], [3.5 4.5],'Color',[.35 .35 .35], 'LineWidth', 4);
plot([3.5 5.5], [5.5 5.5], 'Color',[.35 .35 .35], 'LineWidth', 4);% keep

% DS2 attention stability matrix
subplot(3,2,3);
imagesc(tril(fig3_stabmat_SA_ds2)); hold on
plot(SS2x_ATTN, SS2y_ATTN,'w*','MarkerSize',10,'LineWidth',1.5);
colormap(redwhiteblack);
caxis(clim)
hh = colorbar;
ylabel(hh,"Spearman's rho",'FontSize',16,'Rotation',270);
xticks([1.5:2:9.5]);
xticklabels(Chun_run_names);
yticks([1.5:2:9.5]);
yticklabels(Chun_run_names);
xtickangle(45);
axis square
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',15,'FontWeight','bold','box','off')
% Heavy black lines
plot([0 1.5],[0.5 0.5], 'k', 'LineWidth', 7);
plot([1.5 2.5],[1.5 1.5], 'k', 'LineWidth', 7);
plot([0.5 0.5], [0 2.5],'k', 'LineWidth', 7);
plot([0 2.5],[2.5 2.5], 'k', 'LineWidth', 7);%keep
plot([2.5 2.5], [1.5 2.5],'k', 'LineWidth', 7);
plot([1.5 1.5], [0.5 1.5],'k', 'LineWidth', 7);
% Light gray lines
plot([2.5 3.5],[2.5 2.5], 'Color',[.35 .35 .35], 'LineWidth', 3);
plot([3.5 4.5],[3.5 3.5], 'Color',[.35 .35 .35], 'LineWidth', 3);% keep
plot([2.5 2.5], [2.5 4.5],'Color',[.35 .35 .35], 'LineWidth', 3);% keep
plot([2.5 5.5],[4.5 4.5], 'Color',[.35 .35 .35], 'LineWidth', 3);%keep
plot([3.5 3.5], [2.5 3.5],'Color',[.35 .35 .35], 'LineWidth', 3);%keep
plot([4.5 4.5], [3.5 6.5],'Color',[.35 .35 .35], 'LineWidth', 3);%keep

plot([4.5 7.5],[6.5 6.5], 'Color',[.35 .35 .35], 'LineWidth', 3);%keep
plot([6.5 6.5], [5.5 8.5],'Color',[.35 .35 .35], 'LineWidth', 3);%keep
plot([6.5 9.5],[8.5 8.5], 'Color',[.35 .35 .35], 'LineWidth', 3);%keep
plot([8.5 8.5], [7.5 10.5],'Color',[.35 .35 .35], 'LineWidth', 3);
plot([8.5 10.5],[10.5 10.5], 'Color',[.35 .35 .35], 'LineWidth', 3);
plot([10.5 10.5], [9.5 10.5],'Color',[.35 .35 .35], 'LineWidth', 3);

plot([9.5 9.5], [8.5 9.5],'Color',[.35 .35 .35], 'LineWidth', 3);
plot([7.5 7.5], [6.5 7.5],'Color',[.35 .35 .35], 'LineWidth', 3);
plot([5.5 5.5], [4.5 5.5],'Color',[.35 .35 .35], 'LineWidth', 3);
plot([5.5 6.5], [5.5 5.5],'Color',[.35 .35 .35], 'LineWidth', 3);
plot([7.5 8.5], [7.5 7.5],'Color',[.35 .35 .35], 'LineWidth', 3);
plot([9.5 10.5], [9.5 9.5],'Color',[.35 .35 .35], 'LineWidth', 3);

% DS2 working memory stability matrix
subplot(3,2,4);
imagesc(tril(fig3_stabmat_WM_ds2)); hold on
plot(SS2x_WM, SS2y_WM,'w*','MarkerSize',10,'LineWidth',1.5);
colormap(redwhiteblack);
caxis(clim)
hh = colorbar;
ylabel(hh,"Spearman's rho",'FontSize',16,'Rotation',270);
xticks([1.5:2:9.5]);
xticklabels(Chun_run_names);
yticks([1.5:2:9.5]);
yticklabels(Chun_run_names);
xtickangle(45);
axis square
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',15,'FontWeight','bold','box','off')
% Light gray lines
plot([0 1.5],[0.5 0.5], 'Color',[.35 .35 .35], 'LineWidth', 3);
plot([1.5 2.5],[1.5 1.5], 'Color',[.35 .35 .35], 'LineWidth', 3);
plot([0.5 0.5], [0 2.5], 'Color',[.35 .35 .35], 'LineWidth', 3);
plot([0 2.5],[2.5 2.5], 'Color',[.35 .35 .35], 'LineWidth', 3);%keep
plot([2.5 2.5], [1.5 2.5], 'Color',[.35 .35 .35], 'LineWidth', 3);
plot([1.5 1.5], [0.5 1.5], 'Color',[.35 .35 .35], 'LineWidth', 3);

plot([2.5 3.5],[2.5 2.5], 'Color',[.35 .35 .35], 'LineWidth', 3);
plot([3.5 4.5],[3.5 3.5], 'Color',[.35 .35 .35], 'LineWidth', 3);% keep
plot([2.5 2.5], [2.5 4.5],'Color',[.35 .35 .35], 'LineWidth', 3);% keep
plot([2.5 5.5],[4.5 4.5], 'Color',[.35 .35 .35], 'LineWidth', 3);%keep
plot([3.5 3.5], [2.5 3.5],'Color',[.35 .35 .35], 'LineWidth', 3);%keep
plot([4.5 4.5], [3.5 6.5],'Color',[.35 .35 .35], 'LineWidth', 3);%keep
% Heavy black lines
plot([2.5 3.5],[2.5 2.5], 'k', 'LineWidth', 7);
plot([3.5 4.5],[3.5 3.5], 'k', 'LineWidth', 7);
plot([2.5 2.5], [2.5 4.5], 'k', 'LineWidth', 7);
plot([2.5 4.5],[4.5 4.5], 'k', 'LineWidth', 7);
plot([3.5 3.5], [2.5 3.5], 'k', 'LineWidth', 7);
plot([4.5 4.5], [3.5 4.5], 'k', 'LineWidth', 7);
% Back to light gray lines
plot([4.5 4.5], [4.5 6.5],'Color',[.35 .35 .35], 'LineWidth', 3);%keep
plot([4.5 5.5],[4.5 4.5],'Color',[.35 .35 .35], 'LineWidth', 3);%keep
plot([4.5 7.5],[6.5 6.5], 'Color',[.35 .35 .35], 'LineWidth', 3);%keep
plot([6.5 6.5], [5.5 8.5],'Color',[.35 .35 .35], 'LineWidth', 3);%keep
plot([6.5 9.5],[8.5 8.5], 'Color',[.35 .35 .35], 'LineWidth', 3);%keep
plot([8.5 8.5], [7.5 10.5],'Color',[.35 .35 .35], 'LineWidth', 3);
plot([8.5 10.5],[10.5 10.5], 'Color',[.35 .35 .35], 'LineWidth', 3);
plot([10.5 10.5], [9.5 10.5],'Color',[.35 .35 .35], 'LineWidth', 3);

plot([9.5 9.5], [8.5 9.5],'Color',[.35 .35 .35], 'LineWidth', 3);
plot([7.5 7.5], [6.5 7.5],'Color',[.35 .35 .35], 'LineWidth', 3);
plot([5.5 5.5], [4.5 5.5],'Color',[.35 .35 .35], 'LineWidth', 3);
plot([5.5 6.5], [5.5 5.5],'Color',[.35 .35 .35], 'LineWidth', 3);
plot([7.5 8.5], [7.5 7.5],'Color',[.35 .35 .35], 'LineWidth', 3);
plot([9.5 10.5], [9.5 9.5],'Color',[.35 .35 .35], 'LineWidth', 3);


% DS3 attention stability matrix
subplot(3,2,5);
imagesc(tril(fig3_stabmat_SA_ds3)); hold on
plot(SS3x_ATTN, SS3y_ATTN,'w*','MarkerSize',5,'LineWidth',1.5);
colormap(redwhiteblack);
caxis(clim)
hh = colorbar;
ylabel(hh,"Spearman's rho",'FontSize',16,'Rotation',270);
xticks([1.5:2:15.5,18.5]);
xticklabels(HCP_run_names);
yticks([1.5:2:15.5,18.5]);
yticklabels(HCP_run_names);
xtickangle(45);
axis square
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',15,'FontWeight','bold','box','off')
% Heavy black lines
plot([0 1.5],[0.5 0.5], 'k', 'LineWidth', 4);
plot([1.5 2.5],[1.5 1.5], 'k', 'LineWidth', 4);
plot([0.5 0.5], [0 2.5],'k', 'LineWidth', 4);
plot([0 2.5],[2.5 2.5], 'k', 'LineWidth', 4);%keep
plot([2.5 2.5], [1.5 2.5],'k', 'LineWidth', 4);
plot([1.5 1.5], [0.5 1.5],'k', 'LineWidth', 4);
% Light gray lines
plot([2.5 3.5],[2.5 2.5], 'Color',[.35 .35 .35], 'LineWidth', 2);
plot([3.5 4.5],[3.5 3.5], 'Color',[.35 .35 .35], 'LineWidth', 2);% keep
plot([2.5 2.5], [2.5 4.5],'Color',[.35 .35 .35], 'LineWidth', 2);% keep
plot([2.5 5.5],[4.5 4.5], 'Color',[.35 .35 .35], 'LineWidth', 2);%keep
plot([3.5 3.5], [2.5 3.5],'Color',[.35 .35 .35], 'LineWidth', 2);%keep
plot([4.5 4.5], [3.5 6.5],'Color',[.35 .35 .35], 'LineWidth', 2);%keep

plot([4.5 7.5],[6.5 6.5], 'Color',[.35 .35 .35], 'LineWidth', 2);%keep
plot([6.5 6.5], [5.5 8.5],'Color',[.35 .35 .35], 'LineWidth', 2);%keep
plot([6.5 9.5],[8.5 8.5], 'Color',[.35 .35 .35], 'LineWidth', 2);%keep
plot([8.5 8.5], [7.5 10.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([8.5 11.5],[10.5 10.5], 'Color',[.35 .35 .35], 'LineWidth', 2);
plot([10.5 10.5], [9.5 12.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([12.5 12.5], [11.5 14.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([14.5 14.5], [13.5 16.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([16.5 16.5], [15.5 18.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([18.5 18.5], [17.5 18.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([11.5 11.5], [10.5 11.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([13.5 13.5], [12.5 13.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([15.5 15.5], [14.5 15.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([17.5 17.5], [16.5 17.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([17.5 17.5], [16.5 17.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([18.5 18.5], [17.5 18.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([19.5 19.5], [18.5 19.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([20.5 20.5], [19.5 20.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([16.5 16.5], [16.5 20.5],'Color',[.35 .35 .35], 'LineWidth', 2);% Long line for left side of rest
plot([19.5 19.5], [18.5 19.5],'Color',[.35 .35 .35], 'LineWidth', 2);

plot([9.5 9.5], [8.5 9.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([7.5 7.5], [6.5 7.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([5.5 5.5], [4.5 5.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([5.5 6.5], [5.5 5.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([7.5 8.5], [7.5 7.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([9.5 10.5], [9.5 9.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([11.5 12.5], [11.5 11.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([13.5 14.5], [13.5 13.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([15.5 16.5], [15.5 15.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([17.5 18.5], [17.5 17.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([17.5 18.5], [17.5 17.5],'Color',[.35 .35 .35], 'LineWidth', 2);% SHort lines on top of jagged line down diagonal
plot([18.5 19.5], [18.5 18.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([19.5 20.5], [19.5 19.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([10.5 13.5],[12.5 12.5], 'Color',[.35 .35 .35], 'LineWidth', 2);
plot([12.5 15.5],[14.5 14.5], 'Color',[.35 .35 .35], 'LineWidth', 2);
plot([14.5 17.5],[16.5 16.5], 'Color',[.35 .35 .35], 'LineWidth', 2);
plot([16.5 20.5],[20.5 20.5], 'Color',[.35 .35 .35], 'LineWidth', 2); %Bottom line


% DS3 working memory stability matrix
subplot(3,2,6);
imagesc(tril(fig3_stabmat_WM_ds3)); hold on
plot(SS3x_WM, SS3y_WM,'w*','MarkerSize',5,'LineWidth',1.5);
colormap(redwhiteblack);
caxis(clim)
hh = colorbar;
ylabel(hh,"Spearman's rho",'FontSize',16,'Rotation',270);
xticks([1.5:2:15.5,18.5]);
xticklabels(HCP_run_names);
yticks([1.5:2:15.5,18.5]);
yticklabels(HCP_run_names);
xtickangle(45);
axis square
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',15,'Fontweight','bold','box','off')
% light gray lines
plot([0 1.5],[0.5 0.5], 'Color',[.35 .35 .35], 'LineWidth', 2);
plot([1.5 2.5],[1.5 1.5], 'Color',[.35 .35 .35], 'LineWidth',2);
plot([0.5 0.5], [0 2.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([0 2.5],[2.5 2.5], 'Color',[.35 .35 .35], 'LineWidth', 2);%keep
plot([2.5 2.5], [1.5 2.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([1.5 1.5], [0.5 1.5],'Color',[.35 .35 .35], 'LineWidth', 2);

plot([2.5 3.5],[2.5 2.5], 'Color',[.35 .35 .35], 'LineWidth', 2);
plot([3.5 4.5],[3.5 3.5], 'Color',[.35 .35 .35], 'LineWidth', 2);% keep
plot([2.5 2.5], [2.5 4.5],'Color',[.35 .35 .35], 'LineWidth', 2);% keep
plot([2.5 5.5],[4.5 4.5], 'Color',[.35 .35 .35], 'LineWidth', 2);%keep
plot([3.5 3.5], [2.5 3.5],'Color',[.35 .35 .35], 'LineWidth', 2);%keep
plot([4.5 4.5], [3.5 6.5],'Color',[.35 .35 .35], 'LineWidth', 2);%keep
% Heavy black lines
plot([2.5 3.5],[2.5 2.5], 'k', 'LineWidth', 4);
plot([3.5 4.5],[3.5 3.5], 'k', 'LineWidth', 4);% keep
plot([2.5 2.5], [2.5 4.5], 'k', 'LineWidth', 4);% keep
plot([2.5 4.5],[4.5 4.5], 'k', 'LineWidth', 4);%keep
plot([3.5 3.5], [2.5 3.5], 'k', 'LineWidth', 4);%keep
plot([4.5 4.5], [3.5 4.5], 'k', 'LineWidth', 4);%keep
% Back to light gray lines
plot([4.5 7.5],[6.5 6.5], 'Color',[.35 .35 .35], 'LineWidth', 2);%keep
plot([6.5 6.5], [5.5 8.5],'Color',[.35 .35 .35], 'LineWidth', 2);%keep
plot([6.5 9.5],[8.5 8.5], 'Color',[.35 .35 .35], 'LineWidth', 2);%keep
plot([8.5 8.5], [7.5 10.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([8.5 11.5],[10.5 10.5], 'Color',[.35 .35 .35], 'LineWidth', 2);
plot([10.5 10.5], [9.5 12.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([12.5 12.5], [11.5 14.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([14.5 14.5], [13.5 16.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([16.5 16.5], [15.5 18.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([18.5 18.5], [17.5 18.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([11.5 11.5], [10.5 11.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([13.5 13.5], [12.5 13.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([15.5 15.5], [14.5 15.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([17.5 17.5], [16.5 17.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([17.5 17.5], [16.5 17.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([18.5 18.5], [17.5 18.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([19.5 19.5], [18.5 19.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([20.5 20.5], [19.5 20.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([16.5 16.5], [16.5 20.5],'Color',[.35 .35 .35], 'LineWidth', 2);% Long line for left side of rest
plot([19.5 19.5], [18.5 19.5],'Color',[.35 .35 .35], 'LineWidth', 2);

plot([9.5 9.5], [8.5 9.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([7.5 7.5], [6.5 7.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([5.5 5.5], [4.5 5.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([5.5 6.5], [5.5 5.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([7.5 8.5], [7.5 7.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([9.5 10.5], [9.5 9.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([11.5 12.5], [11.5 11.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([13.5 14.5], [13.5 13.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([15.5 16.5], [15.5 15.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([17.5 18.5], [17.5 17.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([17.5 18.5], [17.5 17.5],'Color',[.35 .35 .35], 'LineWidth', 2);% SHort lines on top of jagged line down diagonal
plot([18.5 19.5], [18.5 18.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([19.5 20.5], [19.5 19.5],'Color',[.35 .35 .35], 'LineWidth', 2);
plot([10.5 13.5],[12.5 12.5], 'Color',[.35 .35 .35], 'LineWidth', 2);
plot([12.5 15.5],[14.5 14.5], 'Color',[.35 .35 .35], 'LineWidth', 2);
plot([14.5 17.5],[16.5 16.5], 'Color',[.35 .35 .35], 'LineWidth', 2);
plot([16.5 20.5],[20.5 20.5], 'Color',[.35 .35 .35], 'LineWidth', 2); %Bottom line


annotation('textbox', [0.27, 0.935, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize',40,'EdgeColor','none')
annotation('textbox', [0.67, 0.935, 0.45, 0.05], 'string', 'B. Working Memory','FontSize',40,'EdgeColor','none')

set(gcf, 'Position',  [100, 0, 1700, 1000])

% Label datasets
annotation('textarrow', [0.13, 0.15],[0.9, 0.15], 'String', {['Dataset 1']},'FontSize',35, 'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)
annotation('textarrow', [0.13, 0.15], [0.6, 0.15], 'String', {['Dataset 2']},'FontSize',35,'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)
annotation('textarrow', [0.13, 0.15], [0.3, 0.15], 'String', {['Dataset 3']},'FontSize',35,'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)
% Label sample size in each box
annotation('textarrow', [0.19, 0.3], [0.9, 0.25], 'String', {["gradCPT (d')"];['N = ' num2str(min([length(fig2_stab_SA_ds1(~isnan(fig2_stab_SA_ds1))), length(fig2_behav_SA_ds1(~isnan(fig2_behav_SA_ds1)))]))]},'FontSize',30,'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)
annotation('textarrow', [0.19, 0.3], [0.6, 0.25], 'String', {["gradCPT (d')"];['N = ' num2str(min([length(fig2_stab_SA_ds2(~isnan(fig2_stab_SA_ds2))), length(fig2_behav_SA_ds2(~isnan(fig2_behav_SA_ds2)))]))]},'FontSize',30,'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)
annotation('textarrow', [0.218, 0.3], [0.3, 0.25], 'String', {["0-back accuracy"];['N = ' num2str(min([length(fig2_stab_SA_ds3(~isnan(fig2_stab_SA_ds3))), length(fig2_behav_SA_ds3(~isnan(fig2_behav_SA_ds3)))]))]},'FontSize',30,'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)
annotation('textarrow', [0.53, 0.3], [0.6, 0.25], 'String', {["Pashler's K"];['N = ' num2str(min([length(fig2_stab_WM_ds2(~isnan(fig2_stab_WM_ds2))), length(fig2_behav_WM_ds2(~isnan(fig2_behav_WM_ds2)))]))]},'FontSize',30,'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)
annotation('textarrow', [0.53, 0.3], [0.3, 0.25], 'String', {["2-back accuracy"];['N = ' num2str(min([length(fig2_stab_WM_ds3(~isnan(fig2_stab_WM_ds3))), length(fig2_behav_WM_ds3(~isnan(fig2_behav_WM_ds3)))]))]},'FontSize',30,'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)



% ----------------------------------------------------------------------------------------------
% ------------------    PLOT FIGURE 3 Stability and Typicality Matrices     --------------------
% ------------------    ALSO Supplementary Figure Opt and Disc Matrices
% ----------------------------------------------------------------------------------------------

%%


load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS1_NetworkStabTypOptDisc_forReview.mat')

NetStab_SA_DS1 = Net_stab_SA;
NetTyp_SA_DS1 = Net_typ_SA;
NetOpt_SA_DS1 = Net_opt_SA;
NetDisc_SA_DS1 = Net_disc_SA;

NetSA_Stab(:,:,1) = NetStab_SA_DS1;
NetSA_Typ(:,:,1) = NetTyp_SA_DS1;
NetSA_Opt(:,:,1) = NetOpt_SA_DS1;
NetSA_Disc(:,:,1) = NetDisc_SA_DS1;

load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS2_NetworkStabTypOptDisc_forReview.mat')

NetStab_SA_DS2 = Net_stab_SA;
NetTyp_SA_DS2 = Net_typ_SA;
NetOpt_SA_DS2 = Net_opt_SA;
NetDisc_SA_DS2 = Net_disc_SA;
NetStab_WM_DS2 = Net_stab_WM;
NetTyp_WM_DS2 = Net_typ_WM;
NetOpt_WM_DS2 = Net_opt_WM;
NetDisc_WM_DS2 = Net_disc_WM;

NetSA_Stab(:,:,2) = NetStab_SA_DS2;
NetWM_Stab(:,:,1) = NetStab_WM_DS2;
NetSA_Typ(:,:,2) = NetTyp_SA_DS2;
NetWM_Typ(:,:,1) = NetTyp_WM_DS2;
NetSA_Opt(:,:,2) = NetOpt_SA_DS2;
NetWM_Opt(:,:,1) = NetOpt_WM_DS2;
NetSA_Disc(:,:,2) = NetDisc_SA_DS2;
NetWM_Disc(:,:,1) = NetDisc_WM_DS2;

load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_NetworkStabTypOptDisc_forReview.mat')

NetStab_SA_DS3 = Net_stab_SA;
NetTyp_SA_DS3 = Net_typ_SA;
NetOpt_SA_DS3 = Net_opt_SA;
NetDisc_SA_DS3 = Net_disc_SA;
NetStab_WM_DS3 = Net_stab_WM;
NetTyp_WM_DS3 = Net_typ_WM;
NetOpt_WM_DS3 = Net_opt_WM;
NetDisc_WM_DS3 = Net_disc_WM;

NetSA_Stab(:,:,3) = NetStab_SA_DS3;
NetWM_Stab(:,:,2) = NetStab_WM_DS3;
NetSA_Typ(:,:,3) = NetTyp_SA_DS3;
NetWM_Typ(:,:,2) = NetTyp_WM_DS3;
NetSA_Opt(:,:,3) = NetOpt_SA_DS3;
NetWM_Opt(:,:,2) = NetOpt_WM_DS3;
NetSA_Disc(:,:,3) = NetDisc_SA_DS3;
NetWM_Disc(:,:,2) = NetDisc_WM_DS3;

meanNetSA_Stab = tanh(nanmean(atanh(NetSA_Stab),3));
meanNetWM_Stab = tanh(nanmean(atanh(NetWM_Stab),3));

meanNetSA_Typ = tanh(nanmean(atanh(NetSA_Typ),3));
meanNetWM_Typ = tanh(nanmean(atanh(NetWM_Typ),3));

meanNetSA_Opt = tanh(nanmean(atanh(NetSA_Opt),3));
meanNetWM_Opt = tanh(nanmean(atanh(NetWM_Opt),3));

meanNetSA_Disc = tanh(nanmean(atanh(NetSA_Disc),3));
meanNetWM_Disc = tanh(nanmean(atanh(NetWM_Disc),3));

% Normalize matrices
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
    
for net1 = 1:max(Shen_network_labels)
    for net2 = 1:max(Shen_network_labels)

%     num = 1:length(Shen_network_labels);
    which1 = Shen_network_labels == net1; % which indices are actually the network(want to make sure we don't include these in MC stats)
   which2 = Shen_network_labels == net2;
   
   
    tmpmat = blankmat;
    tmpmat(which1,which2) = 1;
    tmpmat(which2,which1) = 1;
%    
    net_log = tmpmat(trlind) == 1;
%

net_size = sum(net_log);

NetStab_SA_DS1_Norm(net1,net2) = NetStab_SA_DS1(net1,net2)/net_size;
NetStab_SA_DS2_Norm(net1,net2) = NetStab_SA_DS2(net1,net2)/net_size;
NetStab_SA_DS3_Norm(net1,net2) = NetStab_SA_DS3(net1,net2)/net_size;
NetTyp_SA_DS1_Norm(net1,net2) = NetTyp_SA_DS1(net1,net2)/net_size;
NetTyp_SA_DS2_Norm(net1,net2) = NetTyp_SA_DS2(net1,net2)/net_size;
NetTyp_SA_DS3_Norm(net1,net2) = NetTyp_SA_DS3(net1,net2)/net_size;

NetStab_WM_DS2_Norm(net1,net2) = NetStab_WM_DS2(net1,net2)/net_size;
NetStab_WM_DS3_Norm(net1,net2) = NetStab_WM_DS3(net1,net2)/net_size;
NetTyp_WM_DS2_Norm(net1,net2) = NetTyp_WM_DS2(net1,net2)/net_size;
NetTyp_WM_DS3_Norm(net1,net2) = NetTyp_WM_DS3(net1,net2)/net_size;

    end % for net1
end %for net2
NetSA_Stab_Norm(:,:,1) = NetStab_SA_DS1_Norm;
NetSA_Typ_Norm(:,:,1) = NetTyp_SA_DS1_Norm;
NetSA_Stab_Norm(:,:,2) = NetStab_SA_DS2_Norm;
NetWM_Stab_Norm(:,:,1) = NetStab_WM_DS2_Norm;
NetSA_Typ_Norm(:,:,2) = NetTyp_SA_DS2_Norm;
NetWM_Typ_Norm(:,:,1) = NetTyp_WM_DS2_Norm;
NetSA_Stab_Norm(:,:,3) = NetStab_SA_DS3_Norm;
NetWM_Stab_Norm(:,:,2) = NetStab_WM_DS3_Norm;
NetSA_Typ_Norm(:,:,3) = NetTyp_SA_DS3_Norm;
NetWM_Typ_Norm(:,:,2) = NetTyp_WM_DS3_Norm;

meanNetSA_Stab_Norm = tanh(nanmean(atanh(NetSA_Stab_Norm),3));
meanNetWM_Stab_Norm = tanh(nanmean(atanh(NetWM_Stab_Norm),3));

meanNetSA_Typ_Norm = tanh(nanmean(atanh(NetSA_Typ_Norm),3));
meanNetWM_Typ_Norm = tanh(nanmean(atanh(NetWM_Typ_Norm),3));

% Finished normalizing matrices :)

cc = max(abs([max(max([meanNetSA_Stab_Norm; meanNetWM_Stab_Norm; meanNetSA_Typ_Norm; meanNetWM_Typ_Norm])), min(min([meanNetSA_Stab_Norm; meanNetWM_Stab_Norm; meanNetSA_Typ_Norm; meanNetWM_Typ_Norm]))]));

clims = [-1 1];
clims_Norm = [-1*cc cc];

clims_disc = [-4 4];

cc2 = max(abs([max(max([NetStab_SA_DS1_Norm; NetTyp_SA_DS1_Norm; NetStab_SA_DS2_Norm; NetStab_WM_DS2_Norm; NetTyp_SA_DS2_Norm; NetTyp_WM_DS2_Norm; NetStab_SA_DS3_Norm; NetStab_WM_DS3_Norm; NetTyp_SA_DS3_Norm; NetTyp_WM_DS3_Norm])), min(min([NetStab_SA_DS1_Norm; NetTyp_SA_DS1_Norm; NetStab_SA_DS2_Norm; NetStab_WM_DS2_Norm; NetTyp_SA_DS2_Norm; NetTyp_WM_DS2_Norm; NetStab_SA_DS3_Norm; NetStab_WM_DS3_Norm; NetTyp_SA_DS3_Norm; NetTyp_WM_DS3_Norm]))]));
clims2 = [-1 1];
clims2_Norm = [-1*cc2 cc2];

% Plot 
figure
% Average stability across datasets for sustained attention
ax(1) = subplot(2,4,1);
imagesc(meanNetSA_Stab);
colormap(ax(1), redwhiteblack);
caxis(clims);
axis square;
yticks([1:10])
yticklabels(network_names);
xticks([]);
set(gca, 'FontSize',20)
hh = colorbar;
 ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
hh.Label.Position = [4.5 0];

% Average Normalized stability for sustained attention
ax(2) = subplot(2,4,2);
imagesc(meanNetSA_Stab_Norm);
colormap(ax(2), redwhiteblack);
caxis(clims_Norm);
axis square;
yticks([]);
xticks([]);
set(gca, 'FontSize',20)
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
hh.Label.Position = [5.5 0];

% Average stability across datasets for working memory 
ax(3) = subplot(2,4,3);
imagesc(meanNetWM_Stab);
colormap(ax(3), redwhiteblack);
caxis(clims);
axis square;
yticks([]);
xticks([]);
set(gca, 'FontSize',20)
hh = colorbar;
 ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
hh.Label.Position = [4.5 0];

% Average Normalized stability for working memory 
ax(4) = subplot(2,4,4);
imagesc(meanNetWM_Stab_Norm);
colormap(ax(4), redwhiteblack);
caxis(clims_Norm);
axis square;
yticks([]);
xticks([]);
set(gca, 'FontSize',20)
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
hh.Label.Position = [5.5 0];

% Average typicality across datasets for sustained attention
ax(5) = subplot(2,4,5);
imagesc(real(meanNetSA_Typ));
colormap(ax(5), bluewhiteblack);
caxis(clims);
axis square;
yticks([1:10])
yticklabels(network_names);
xticks([1:10]);
xticklabels(network_names);
xtickangle(45);
set(gca, 'FontSize',20)
hh = colorbar;
 ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
hh.Label.Position = [4.5 0];

% Average Normalized typicality for sustained attention
ax(6) = subplot(2,4,6);
imagesc(real(meanNetSA_Typ_Norm));
colormap(ax(6), bluewhiteblack);
caxis(clims_Norm);
axis square;
yticks([]);
xticks([1:10]);
xticklabels(network_names);
xtickangle(45);
set(gca, 'FontSize',20)
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
hh.Label.Position = [5.5 0];

% Average typicality across datasets for working memory 
ax(7) = subplot(2,4,7);
imagesc(real(meanNetWM_Typ));
colormap(ax(7), bluewhiteblack);
caxis(clims);
axis square;
yticks([]);
xticks([1:10]);
xticklabels(network_names);
xtickangle(45);
set(gca, 'FontSize',20)
hh = colorbar;
 ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
hh.Label.Position = [4.5 0];

% Average Normalized typicality for working memory 
ax(8) = subplot(2,4,8);
imagesc(real(meanNetWM_Typ_Norm));
colormap(ax(8), bluewhiteblack);
caxis(clims_Norm);
axis square;
yticks([]);
xticks([1:10]);
xticklabels(network_names);
xtickangle(45);
set(gca, 'FontSize',20)
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
hh.Label.Position = [5.5 0];
set(gcf,'Position',[0 0 1800 1000])

sz = 0.28; % how big each plot is (standardize to make it easy to change all subplot sizes at once)

annotation('textbox', [0.2, 0.95, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize',40,'EdgeColor','none')
annotation('textbox', [0.68, 0.95, 0.45, 0.05], 'string', 'B. Working Memory','FontSize',40,'EdgeColor','none')


annotation('textbox', [0.16, 0.87, 0.45, 0.05], 'string', "Pearson's r",'FontSize',25,'EdgeColor','none')
annotation('textbox', [0.61, 0.87, 0.45, 0.05], 'string', "Pearson's r",'FontSize',25,'EdgeColor','none')

annotation('textbox', [0.33, 0.87, 0.45, 0.05], 'string', {['Normalized by network size']}','FontSize',25,'EdgeColor','none')
annotation('textbox', [0.78, 0.87, 0.45, 0.05], 'string', {['Normalized by network size']},'FontSize',25,'EdgeColor','none')

% stability and typ
annotation('textarrow', [0.01, 0.8],[.78 .9], 'String', {['Stability']},'FontSize',35, 'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)
annotation('textarrow', [0.01, 0.8],[.42 .9], 'String', {['Typicality']},'FontSize',35,'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)

ha=get(gcf,'children');
 set(ha(2),'position',[.71 .2 sz sz]) % bottom right corner
 set(ha(4),'position',[.5 .2 sz sz])
 set(ha(6),'position',[.26 .2 sz sz])
 set(ha(8),'position',[.05 .2 sz sz])
% 
  set(ha(10),'position',[.71 .55 sz sz]) % top right
  set(ha(12),'position',[.5 .55 sz sz])
  set(ha(14),'position',[.26 .55 sz sz])
  set(ha(16),'position',[.05 .55 sz sz]) % top left



%% Repeat for Disc and optimality for review

for net1 = 1:max(Shen_network_labels)
    for net2 = 1:max(Shen_network_labels)

%     num = 1:length(Shen_network_labels);
    which1 = Shen_network_labels == net1; % which indices are actually the network(want to make sure we don't include these in MC stats)
   which2 = Shen_network_labels == net2;
   
   
    tmpmat = blankmat;
    tmpmat(which1,which2) = 1;
    tmpmat(which2,which1) = 1;
%    
    net_log = tmpmat(trlind) == 1;
%

net_size = sum(net_log);

NetOpt_SA_DS1_Norm(net1,net2) = NetOpt_SA_DS1(net1,net2)/net_size;
NetOpt_SA_DS2_Norm(net1,net2) = NetOpt_SA_DS2(net1,net2)/net_size;
NetOpt_SA_DS3_Norm(net1,net2) = NetOpt_SA_DS3(net1,net2)/net_size;
NetDisc_SA_DS1_Norm(net1,net2) = NetDisc_SA_DS1(net1,net2)/net_size;
NetDisc_SA_DS2_Norm(net1,net2) = NetDisc_SA_DS2(net1,net2)/net_size;
NetDisc_SA_DS3_Norm(net1,net2) = NetDisc_SA_DS3(net1,net2)/net_size;

NetOpt_WM_DS2_Norm(net1,net2) = NetOpt_WM_DS2(net1,net2)/net_size;
NetOpt_WM_DS3_Norm(net1,net2) = NetOpt_WM_DS3(net1,net2)/net_size;
NetDisc_WM_DS2_Norm(net1,net2) = NetDisc_WM_DS2(net1,net2)/net_size;
NetDisc_WM_DS3_Norm(net1,net2) = NetDisc_WM_DS3(net1,net2)/net_size;

    end % for net1
end %for net2
NetSA_Opt_Norm(:,:,1) = NetOpt_SA_DS1_Norm;
NetSA_Disc_Norm(:,:,1) = NetDisc_SA_DS1_Norm;
NetSA_Opt_Norm(:,:,2) = NetOpt_SA_DS2_Norm;
NetWM_Opt_Norm(:,:,1) = NetOpt_WM_DS2_Norm;
NetSA_Disc_Norm(:,:,2) = NetDisc_SA_DS2_Norm;
NetWM_Disc_Norm(:,:,1) = NetDisc_WM_DS2_Norm;
NetSA_Opt_Norm(:,:,3) = NetOpt_SA_DS3_Norm;
NetWM_Opt_Norm(:,:,2) = NetOpt_WM_DS3_Norm;
NetSA_Disc_Norm(:,:,3) = NetDisc_SA_DS3_Norm;
NetWM_Disc_Norm(:,:,2) = NetDisc_WM_DS3_Norm;

meanNetSA_Opt_Norm = tanh(nanmean(atanh(NetSA_Opt_Norm),3));
meanNetWM_Opt_Norm = tanh(nanmean(atanh(NetWM_Opt_Norm),3));

meanNetSA_Disc_Norm = tanh(nanmean(atanh(NetSA_Disc_Norm),3));
meanNetWM_Disc_Norm = tanh(nanmean(atanh(NetWM_Disc_Norm),3));

% Finished normalizing matrices :)

cc = max(abs([max(max([meanNetSA_Opt_Norm; meanNetWM_Opt_Norm; meanNetSA_Disc_Norm; meanNetWM_Disc_Norm])), min(min([meanNetSA_Opt_Norm; meanNetWM_Opt_Norm; meanNetSA_Disc_Norm; meanNetWM_Disc_Norm]))]));

clims = [-1 1];
clims_Norm = [-1*cc cc];

cc2 = max(abs([max(max([NetOpt_SA_DS1_Norm; NetDisc_SA_DS1_Norm; NetOpt_SA_DS2_Norm; NetOpt_WM_DS2_Norm; NetDisc_SA_DS2_Norm; NetDisc_WM_DS2_Norm; NetOpt_SA_DS3_Norm; NetOpt_WM_DS3_Norm; NetDisc_SA_DS3_Norm; NetDisc_WM_DS3_Norm])), min(min([NetOpt_SA_DS1_Norm; NetDisc_SA_DS1_Norm; NetOpt_SA_DS2_Norm; NetOpt_WM_DS2_Norm; NetDisc_SA_DS2_Norm; NetDisc_WM_DS2_Norm; NetOpt_SA_DS3_Norm; NetOpt_WM_DS3_Norm; NetDisc_SA_DS3_Norm; NetDisc_WM_DS3_Norm]))]));
clims2 = [-1 1];
clims2_Norm = [-1*cc2 cc2];

% Plot 
figure
% Average Optility across datasets for sustained attention
ax(1) = subplot(2,4,1);
imagesc(meanNetSA_Opt);
colormap(ax(1), lightbluewhiteblack);
caxis(clims);
axis square;
yticks([1:10])
yticklabels(network_names);
xticks([]);
set(gca, 'FontSize',20)
hh = colorbar;
 ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
hh.Label.Position = [4.5 0];

% Average Normalized Optility for sustained attention
ax(2) = subplot(2,4,2);
imagesc(meanNetSA_Opt_Norm);
colormap(ax(2), lightbluewhiteblack);
caxis(clims_Norm);
axis square;
yticks([]);
xticks([]);
set(gca, 'FontSize',20)
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
hh.Label.Position = [5.5 0];

% Average Optility across datasets for working memory 
ax(3) = subplot(2,4,3);
imagesc(meanNetWM_Opt);
colormap(ax(3), lightbluewhiteblack);
caxis(clims);
axis square;
yticks([]);
xticks([]);
set(gca, 'FontSize',20)
hh = colorbar;
 ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
hh.Label.Position = [4.5 0];

% Average Normalized Optility for working memory 
ax(4) = subplot(2,4,4);
imagesc(meanNetWM_Opt_Norm);
colormap(ax(4), lightbluewhiteblack);
caxis(clims_Norm);
axis square;
yticks([]);
xticks([]);
set(gca, 'FontSize',20)
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
hh.Label.Position = [5.5 0];

% Average Discicality across datasets for sustained attention
ax(5) = subplot(2,4,5);
imagesc(real(meanNetSA_Disc));
colormap(ax(5), greenwhiteblack);
caxis([-4 4]);
axis square;
yticks([1:10])
yticklabels(network_names);
xticks([1:10]);
xticklabels(network_names);
xtickangle(45);
set(gca, 'FontSize',20)
hh = colorbar;
 ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
hh.Label.Position = [4.5 0];

% Average Normalized Discicality for sustained attention
ax(6) = subplot(2,4,6);
imagesc(real(meanNetSA_Disc_Norm));
colormap(ax(6), greenwhiteblack);
caxis(clims_Norm);
axis square;
yticks([]);
xticks([1:10]);
xticklabels(network_names);
xtickangle(45);
set(gca, 'FontSize',20)
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
hh.Label.Position = [5.5 0];

% Average Discicality across datasets for working memory 
ax(7) = subplot(2,4,7);
imagesc(real(meanNetWM_Disc));
colormap(ax(7), greenwhiteblack);
caxis([-4 4]);
axis square;
yticks([]);
xticks([1:10]);
xticklabels(network_names);
xtickangle(45);
set(gca, 'FontSize',20)
hh = colorbar;
 ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
hh.Label.Position = [4.5 0];

% Average Normalized Discicality for working memory 
ax(8) = subplot(2,4,8);
imagesc(real(meanNetWM_Disc_Norm));
colormap(ax(8), greenwhiteblack);
caxis(clims_Norm);
axis square;
yticks([]);
xticks([1:10]);
xticklabels(network_names);
xtickangle(45);
set(gca, 'FontSize',20)
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
hh.Label.Position = [5.5 0];
set(gcf,'Position',[0 0 1800 1000])

sz = 0.28; % how big each plot is (standardize to make it easy to change all subplot sizes at once)

annotation('textbox', [0.2, 0.95, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize',40,'EdgeColor','none')
annotation('textbox', [0.68, 0.95, 0.45, 0.05], 'string', 'B. Working Memory','FontSize',40,'EdgeColor','none')


annotation('textbox', [0.16, 0.87, 0.45, 0.05], 'string', "Pearson's r",'FontSize',25,'EdgeColor','none')
annotation('textbox', [0.61, 0.87, 0.45, 0.05], 'string', "Pearson's r",'FontSize',25,'EdgeColor','none')

annotation('textbox', [0.33, 0.87, 0.45, 0.05], 'string', {['Normalized by network size']}','FontSize',25,'EdgeColor','none')
annotation('textbox', [0.78, 0.87, 0.45, 0.05], 'string', {['Normalized by network size']},'FontSize',25,'EdgeColor','none')

% Optility and Disc
annotation('textarrow', [0.01, 0.8],[.78 .9], 'String', {['Optimality']},'FontSize',35, 'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)
annotation('textarrow', [0.01, 0.8],[.42 .9], 'String', {['Discriminability']},'FontSize',35,'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)

ha=get(gcf,'children');
 set(ha(2),'position',[.71 .2 sz sz]) % bottom right corner
 set(ha(4),'position',[.5 .2 sz sz])
 set(ha(6),'position',[.26 .2 sz sz])
 set(ha(8),'position',[.05 .2 sz sz])
% 
  set(ha(10),'position',[.71 .55 sz sz]) % top right
  set(ha(12),'position',[.5 .55 sz sz])
  set(ha(14),'position',[.26 .55 sz sz])
  set(ha(16),'position',[.05 .55 sz sz]) % top left





%%

% Plot supplementary figure S3 - Network stability and typicality by
% dataset

figure
% Stability ================================
% DS1 Stability 
subplot(3,2,1);
imagesc(NetStab_SA_DS1);
colormap(redwhiteblack);
caxis(clims);
axis square;
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 Stability ATTN
subplot(3,2,3);
imagesc(NetStab_SA_DS2);
colormap(redwhiteblack);
caxis(clims);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
axis square;
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 Stability WM
subplot(3,2,4);
imagesc(NetStab_WM_DS2);
colormap(redwhiteblack);
caxis(clims);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
axis square;
xticklabels([]);
yticklabels([]);
% DS3 Stability ATTN
subplot(3,2,5);
imagesc(NetStab_SA_DS3);
colormap(redwhiteblack);
caxis(clims);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
axis square;
xticklabels(network_names);
xticks([1:10]);
xtickangle(45)
yticklabels(network_names);
yticks([1:10]);
% DS3 Stability WM
subplot(3,2,6);
imagesc(NetStab_WM_DS3);
colormap(redwhiteblack);
caxis(clims);
axis square;
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
xticklabels(network_names);
xticks([1:10]);
xtickangle(45)
yticklabels([]);
%set(gca, 'FontSize',25);
annotation('textbox', [0.14, 0.94, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize',25,'EdgeColor','none')
annotation('textbox', [0.58, 0.94, 0.45, 0.05], 'string', 'B. Working Memory','FontSize',25,'EdgeColor','none')
set(gcf, 'Position',  [00, 100, 850, 1000]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)

figure
% Stability Normalized
subplot(3,2,1);
imagesc(NetStab_SA_DS1_Norm);
colormap(redwhiteblack);
caxis(clims2_Norm);
axis square;
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 Stability ATTN normalized
subplot(3,2,3);
imagesc(NetStab_SA_DS2_Norm);
colormap(redwhiteblack);
caxis(clims2_Norm);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 Stability WM normalized
subplot(3,2,4);
imagesc(NetStab_WM_DS2_Norm);
colormap(redwhiteblack);
caxis(clims2_Norm);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
xticklabels([]);
yticklabels([]);
% DS3 Stability ATTN normalized
subplot(3,2,5);
imagesc(NetStab_SA_DS3_Norm);
colormap(redwhiteblack);
caxis(clims2_Norm);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
xticklabels(network_names);
xticks([1:10]);
yticklabels(network_names);
yticks([1:10]);
% DS3 Stability WM normalized
subplot(3,2,6);
imagesc(NetStab_WM_DS3_Norm);
colormap(redwhiteblack);
caxis(clims2_Norm);
axis square;
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
xticklabels(network_names);
xticks([1:10]);
xtickangle(45);
yticklabels([]);
%set(gca, 'FontSize',25);
%set(gca, 'FontSize',25);
annotation('textbox', [0.14, 0.94, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize',25,'EdgeColor','none')
annotation('textbox', [0.58, 0.94, 0.45, 0.05], 'string', 'B. Working Memory','FontSize',25,'EdgeColor','none')
set(gcf, 'Position',  [00, 100, 850, 1000]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)

figure
% Typicality 
% DS1 Typicality
subplot(3,2,1);
imagesc(real(NetTyp_SA_DS1));
colormap(bluewhiteblack);
caxis(clims);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
axis square;
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 TYp ATTN
subplot(3,2,3);
imagesc(real(NetTyp_SA_DS2));
colormap(bluewhiteblack);
caxis(clims);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
axis square;
% DS2 TYp WM
subplot(3,2,4);
imagesc(real(NetTyp_WM_DS2));
colormap(bluewhiteblack);
caxis(clims);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
axis square;
xticklabels([]);
yticklabels([]);
% DS3 Typicality ATTN
subplot(3,2,5);
imagesc(real(NetTyp_SA_DS3));
colormap(bluewhiteblack);
caxis(clims);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
axis square;
xticklabels(network_names);
xticks([1:10]);
yticklabels(network_names);
yticks([1:10]);
% DS3 Typicality WM
subplot(3,2,6);
imagesc(real(NetTyp_WM_DS3));
colormap(bluewhiteblack);
caxis(clims);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
axis square;
xticklabels(network_names);
xticks([1:10]);
xtickangle(45);
yticklabels([]);
%set(gca, 'FontSize',25);
annotation('textbox', [0.14, 0.94, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize',25,'EdgeColor','none')
annotation('textbox', [0.58, 0.94, 0.45, 0.05], 'string', 'B. Working Memory','FontSize',25,'EdgeColor','none')
set(gcf, 'Position',  [00, 100, 850, 1000]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)

figure
% Typicality Normalized
% DS1 Typicality norm
subplot(3,2,1);
imagesc(real(NetTyp_SA_DS1_Norm));
colormap(bluewhiteblack);
caxis(clims2_Norm);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 TYp ATTN normalized
subplot(3,2,3);
imagesc(real(NetTyp_SA_DS2_Norm));
colormap(bluewhiteblack);
caxis(clims2_Norm);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 TYp WM normalized
subplot(3,2,4);
imagesc(real(NetTyp_WM_DS2_Norm));
colormap(bluewhiteblack);
caxis(clims2_Norm);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
xticklabels([]);
yticklabels([]);
% DS3 Typicality ATTN normalized
subplot(3,2,5);
imagesc(real(NetTyp_SA_DS3_Norm));
colormap(bluewhiteblack);
caxis(clims2_Norm);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
xticks([1:10]);
xticklabels(network_names);
xticks([1:10]);
yticklabels(network_names);
yticks([1:10]);
% DS3 Typicality WM normalized
subplot(3,2,6);
imagesc(real(NetTyp_WM_DS3_Norm));
colormap(bluewhiteblack);
caxis(clims2_Norm);
set(gca, 'FontSize',25);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
xticklabels(network_names);
xticks([1:10]);
xtickangle(45)
yticklabels([]);
%set(gca, 'FontSize',25);
annotation('textbox', [0.14, 0.94, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize',25,'EdgeColor','none')
annotation('textbox', [0.58, 0.94, 0.45, 0.05], 'string', 'B. Working Memory','FontSize',25,'EdgeColor','none')
set(gcf, 'Position',  [00, 100, 850, 1000]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)

%%
figure
% Optimality ================================
% DS1 Optimality
subplot(3,2,1);
imagesc(NetOpt_SA_DS1);
colormap(lightbluewhiteblack);
caxis(clims);
axis square;
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 Optimality ATTN
subplot(3,2,3);
imagesc(NetOpt_SA_DS2);
colormap(lightbluewhiteblack);
caxis(clims);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
axis square;
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 Optimality WM
subplot(3,2,4);
imagesc(NetOpt_WM_DS2);
colormap(lightbluewhiteblack);
caxis(clims);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
axis square;
xticklabels([]);
yticklabels([]);
% DS3 Optimality ATTN
subplot(3,2,5);
imagesc(NetOpt_SA_DS3);
colormap(lightbluewhiteblack);
caxis(clims);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
axis square;
xticklabels(network_names);
xticks([1:10]);
xtickangle(45)
yticklabels(network_names);
yticks([1:10]);
% DS3 Optimality WM
subplot(3,2,6);
imagesc(NetOpt_WM_DS3);
colormap(lightbluewhiteblack);
caxis(clims);
axis square;
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
xticklabels(network_names);
xticks([1:10]);
xtickangle(45)
yticklabels([]);
%set(gca, 'FontSize',25);
annotation('textbox', [0.14, 0.94, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize',25,'EdgeColor','none')
annotation('textbox', [0.58, 0.94, 0.45, 0.05], 'string', 'B. Working Memory','FontSize',25,'EdgeColor','none')
set(gcf, 'Position',  [00, 100, 850, 1000]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)

figure
% Optimality Normalized
subplot(3,2,1);
imagesc(NetOpt_SA_DS1_Norm);
colormap(lightbluewhiteblack);
caxis(clims2_Norm);
axis square;
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 Optimality ATTN normalized
subplot(3,2,3);
imagesc(NetOpt_SA_DS2_Norm);
colormap(lightbluewhiteblack);
caxis(clims2_Norm);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 Optimality WM normalized
subplot(3,2,4);
imagesc(NetOpt_WM_DS2_Norm);
colormap(lightbluewhiteblack);
caxis(clims2_Norm);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
xticklabels([]);
yticklabels([]);
% DS3 Optimality ATTN normalized
subplot(3,2,5);
imagesc(NetOpt_SA_DS3_Norm);
colormap(lightbluewhiteblack);
caxis(clims2_Norm);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
xticklabels(network_names);
xticks([1:10]);
yticklabels(network_names);
yticks([1:10]);
% DS3 Optimality WM normalized
subplot(3,2,6);
imagesc(NetOpt_WM_DS3_Norm);
colormap(lightbluewhiteblack);
caxis(clims2_Norm);
axis square;
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
xticklabels(network_names);
xticks([1:10]);
xtickangle(45);
yticklabels([]);
%set(gca, 'FontSize',25);
%set(gca, 'FontSize',25);
annotation('textbox', [0.14, 0.94, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize',25,'EdgeColor','none')
annotation('textbox', [0.58, 0.94, 0.45, 0.05], 'string', 'B. Working Memory','FontSize',25,'EdgeColor','none')
set(gcf, 'Position',  [00, 100, 850, 1000]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)



figure
% Discriminability ================================
% DS1 Discriminability
subplot(3,2,1);
imagesc(NetDisc_SA_DS1);
colormap(greenwhiteblack);
caxis(clims_disc);
axis square;
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 Discriminability ATTN
subplot(3,2,3);
imagesc(NetDisc_SA_DS2);
colormap(greenwhiteblack);
caxis(clims_disc);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
axis square;
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 Discriminability WM
subplot(3,2,4);
imagesc(NetDisc_WM_DS2);
colormap(greenwhiteblack);
caxis(clims_disc);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
axis square;
xticklabels([]);
yticklabels([]);
% DS3 Discriminability ATTN
subplot(3,2,5);
imagesc(NetDisc_SA_DS3);
colormap(greenwhiteblack);
caxis(clims_disc);
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
axis square;
xticklabels(network_names);
xticks([1:10]);
xtickangle(45)
yticklabels(network_names);
yticks([1:10]);
% DS3 Discriminability WM
subplot(3,2,6);
imagesc(NetDisc_WM_DS3);
colormap(greenwhiteblack);
caxis(clims_disc);
axis square;
hh = colorbar;
ylabel(hh,"Pearson's r",'FontSize',16,'Rotation',270);
xticklabels(network_names);
xticks([1:10]);
xtickangle(45)
yticklabels([]);
%set(gca, 'FontSize',25);
annotation('textbox', [0.14, 0.94, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize',25,'EdgeColor','none')
annotation('textbox', [0.58, 0.94, 0.45, 0.05], 'string', 'B. Working Memory','FontSize',25,'EdgeColor','none')
set(gcf, 'Position',  [00, 100, 850, 1000]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)

figure
% Discriminability Normalized
subplot(3,2,1);
imagesc(NetDisc_SA_DS1_Norm);
colormap(greenwhiteblack);
caxis(clims2_Norm);
axis square;
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 Discriminability ATTN normalized
subplot(3,2,3);
imagesc(NetDisc_SA_DS2_Norm);
colormap(greenwhiteblack);
caxis(clims2_Norm);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
yticklabels(network_names);
yticks([1:10]);
xticklabels([]);
% DS2 Discriminability WM normalized
subplot(3,2,4);
imagesc(NetDisc_WM_DS2_Norm);
colormap(greenwhiteblack);
caxis(clims2_Norm);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
xticklabels([]);
yticklabels([]);
% DS3 Discriminability ATTN normalized
subplot(3,2,5);
imagesc(NetDisc_SA_DS3_Norm);
colormap(greenwhiteblack);
caxis(clims2_Norm);
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
axis square;
xticklabels(network_names);
xticks([1:10]);
yticklabels(network_names);
yticks([1:10]);
% DS3 Discriminability WM normalized
subplot(3,2,6);
imagesc(NetDisc_WM_DS3_Norm);
colormap(greenwhiteblack);
caxis(clims2_Norm);
axis square;
hh = colorbar;
ylabel(hh,"Normalized r",'FontSize',16,'Rotation',270);
xticklabels(network_names);
xticks([1:10]);
xtickangle(45);
yticklabels([]);
%set(gca, 'FontSize',25);
%set(gca, 'FontSize',25);
annotation('textbox', [0.14, 0.94, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize',25,'EdgeColor','none')
annotation('textbox', [0.58, 0.94, 0.45, 0.05], 'string', 'B. Working Memory','FontSize',25,'EdgeColor','none')
set(gcf, 'Position',  [00, 100, 850, 1000]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)


% ----------------------------------------------------------------------------------------------
% ------------------    PLOT FIGURE 5 Stability Network Anatomy     --------------------
% ----------------------------------------------------------------------------------------------
%%
addpath('/Users/annacorriveau/Documents/GitHub/Violinplot-Matlab');

% Load lesioned values
load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS1_StabAnat.mat');
% for n = 1:length(r_stab_ds1); %loop through networks to save relevant values
lesioned_stab_ds1_attn = cell2mat(r_stab_ds1);
lesioned_typ_ds1_attn = cell2mat(r_typ_ds1);
% end

load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS2_StabAnat.mat');
%load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS2_StabAnat_resavetyp.mat');
for n = 1:size(corr_stab_ATTN,2); %loop through networks to save relevant values
lesioned_stab_ds2_attn(n) = corr_stab_ATTN{n}(1,2);
lesioned_stab_ds2_wm(n) = corr_stab_WM{n}(3,4);

lesioned_typ_ds2_attn(n) = corr_typ_ATTN{n}(1,2);
lesioned_typ_ds2_wm(n) = corr_typ_WM{n}(3,4);
end


load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_StabAnat.mat');
%load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_StabAnat_resavetyp.mat');
for n = 1:size(corr_stab_ATTN,2); %loop through networks to save relevant values
 lesioned_stab_ds3_attn(n) = corr_stab_ATTN{n}(1,2);
 lesioned_stab_ds3_wm(n) = corr_stab_WM{n}(3,4);

    lesioned_typ_ds3_attn(n) = corr_typ_ATTN{n}(1,2);
    lesioned_typ_ds3_wm(n) = corr_typ_WM{n}(3,4);
 end

% Load MC shuffled null distributions
% load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS1_NullStabConcat.mat']);
% null_stab_ds1_attn = cell2mat(ds1_stab_concat);
load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS1_NullStabAnat_clean3SA.mat')
null_stab_ds1_attn = cell2mat(r_stab_ds1);

load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS1_NullTypAnat.mat']);
null_typ_ds1_attn = cell2mat(r_typ_ds1);
null_typ_ds1_attn = null_typ_ds1_attn(:,1:250);



load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS2_NullStabAnat1SA.mat');
for i = 1:size(corr_stab_ATTN,2)
    for p = 1:size(corr_stab_ATTN,1)
        null_stab_ds2_attn(p,i) = corr_stab_ATTN{p,i}(1,2);
    end
end
%null_stab_ds2_attn = ds2_stab_ATTN_concat;
load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS2_NullStabAnat1WM.mat');
for i = 1:size(corr_stab_WM,2)
    for p = 1:size(corr_stab_WM,1)
        null_stab_ds2_wm(p,i) = corr_stab_WM{p,i}(3,4);
    end
end
% null_stab_ds2_wm = ds2_stab_WM_concat;

load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS2_NullTypAnat1SA.mat');
for i = 1:size(corr_typ_ATTN,2)
    for p = 1:size(corr_typ_ATTN,1)
        null_typ_ds2_attn(p,i) = corr_typ_ATTN{p,i}(1,2);
    end
end
%null_typ_ds2_attn = null_typ_ds2_attn(:,1:250);
%null_typ_ds2_attn = ds2_typ_ATTN_concat;
load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS2_NullTypAnat1WM.mat');
for i = 1:size(corr_typ_WM,2)
    for p = 1:size(corr_typ_WM,1)
        null_typ_ds2_wm(p,i) = corr_typ_WM{p,i}(3,4);
    end
end
%null_typ_ds2_wm = null_typ_ds2_wm(:,1:250);
%null_typ_ds2_wm = ds2_typ_WM_concat;


 % DS3 stability
load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_NullStabAnat1SA.mat');

for i = 1:size(corr_stab_ATTN,2)
    for p = 1:size(corr_stab_ATTN,1)
        null_stab_ds3_attn(p,i) = corr_stab_ATTN{p,i}(1,2);
    end
end

load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_NullStabAnat1WM.mat');

for i = 1:size(corr_stab_WM,2)-1
    for p = 1:size(corr_stab_WM,1)
        null_stab_ds3_wm(p,i) = corr_stab_WM{p,i}(3,4);
    end
end
ct = size(null_stab_ds3_wm,2);
load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_NullStabAnat3WM.mat');
for i = 1:size(corr_stab_WM,2)-1
    for p = 1:size(corr_stab_WM,1)
        null_stab_ds3_wm(p,i+ct) = corr_stab_WM{p,i}(3,4);
    end
end
null_stab_ds3_wm = null_stab_ds3_wm(:,1:1000);

% Ds3 typicality
load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_NullTypAnat1SA.mat');
for i = 1:size(corr_typ_ATTN,2)-1
    for p = 1:size(corr_typ_ATTN,1)
        null_typ_ds3_attn(p,i) = corr_typ_ATTN{p,i}(1,2);
    end
end
ct = size(null_typ_ds3_attn,2);
load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_NullTypAnat__3SA.mat');
for i = 1:size(corr_typ_ATTN,2)-1
    for p = 1:size(corr_typ_ATTN,1)
        null_typ_ds3_attn(p, i + ct) = corr_typ_ATTN{p,i}(1,2);
    end
end
ct = size(null_typ_ds3_attn,2);
load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_NullTypAnat__SA4.mat');
for i = 1:size(corr_typ_ATTN,2)-1
    for p = 1:size(corr_typ_ATTN,1)
        null_typ_ds3_attn(p, i + ct) = corr_typ_ATTN{p,i}(1,2);
    end
end
null_typ_ds3_attn = null_typ_ds3_attn(:,1:250);

load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_NullTypAnat1WM.mat');
for i = 1:size(corr_typ_WM,2)-1
    for p = 1:size(corr_typ_WM,1)
        null_typ_ds3_wm(p,i) = corr_typ_WM{p,i}(3,4);
    end
end
ct = size(null_typ_ds3_wm,2);
load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_NullTypAnat__2WM.mat');
for i = 1:size(corr_typ_WM,2)-1
    for p = 1:size(corr_typ_WM,1)
        null_typ_ds3_wm(p, i + ct) = corr_typ_WM{p,i}(3,4);
    end
end
ct = size(null_typ_ds3_wm,2);
load('/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/outputs/DS3_NullTypAnat__WM4.mat');
for i = 1:size(corr_typ_WM,2)-1
    for p = 1:size(corr_typ_WM,1)
        null_typ_ds3_wm(p, i + ct) = corr_typ_WM{p,i}(3,4);
    end
end

null_typ_ds3_wm = null_typ_ds3_wm(:,1:250);

% find p-values (how many values are more extreme than our observed
% lesioned values?)

% Find if network significantly contributed to stability/typicality's
% correlation with performance
%     diff_ds1_stab_attn = fig5_ds1_stab_attn - lesioned_stab_ds1_attn;
%     diff_ds2_stab_attn = fig5_ds2_stab_attn - lesioned_stab_ds2_attn;
%     diff_ds3_stab_attn = fig5_ds3_stab_attn - lesioned_stab_ds3_attn;
% 
%     diff_ds2_stab_wm = fig5_ds2_stab_wm - lesioned_stab_ds2_wm;
%     diff_ds3_stab_wm = fig5_ds3_stab_wm - lesioned_stab_ds3_wm;

for n = 1:length(lesioned_stab_ds1_attn)

    % DS1 attention stability 

    if (min([sum(null_stab_ds1_attn(n,:) < lesioned_stab_ds1_attn(n)),   sum(null_stab_ds1_attn(n,:) > lesioned_stab_ds1_attn(n))]) + 1)/(size(null_stab_ds1_attn,2) + 1) < .025;
        sig_ds1_stabanat_attn(n) = 1; % significant
    else
        sig_ds1_stabanat_attn(n) = NaN; % not significant
    end


 % DS2 attention stability

    if (min([sum(null_stab_ds2_attn(n,:) < lesioned_stab_ds2_attn(n)),   sum(null_stab_ds2_attn(n,:) > lesioned_stab_ds2_attn(n))]) + 1)/(size(null_stab_ds2_attn,2) + 1) < .025;
        sig_ds2_stabanat_attn(n) = 1; % significant
    else
        sig_ds2_stabanat_attn(n) = NaN; % not significant
    end

 % DS2 working mem stability

    if (min([sum(null_stab_ds2_wm(n,:) < lesioned_stab_ds2_wm(n)),   sum(null_stab_ds2_wm(n,:) > lesioned_stab_ds2_wm(n))]) + 1)/(size(null_stab_ds2_wm,2) + 1) < .025;
        sig_ds2_stabanat_wm(n) = 1; % significant
    else
        sig_ds2_stabanat_wm(n) = NaN; % not significant
    end


     % DS3 attention stability

    if (min([sum(null_stab_ds3_attn(n,:) < lesioned_stab_ds3_attn(n)),   sum(null_stab_ds3_attn(n,:) > lesioned_stab_ds3_attn(n))]) + 1)/(size(null_stab_ds3_attn,2) + 1) < .025;
        sig_ds3_stabanat_attn(n) = 1; % significant
    else
        sig_ds3_stabanat_attn(n) = NaN; % not significant
    end

 % DS3 working mem stability

    if (min([sum(null_stab_ds3_wm(n,:) < lesioned_stab_ds3_wm(n)),   sum(null_stab_ds3_wm(n,:) > lesioned_stab_ds3_wm(n))]) + 1)/(size(null_stab_ds3_wm,2) + 1) < .025;
        sig_ds3_stabanat_wm(n) = 1; % significant
    else
        sig_ds3_stabanat_wm(n) = NaN; % not significant
    end


end % for n = 


%%
% Create violin plots for null values
xlims = [0 11];

ds1_SA_ylim = [0.2 0.5];
figure
subplot(3,2,1);
yline([fig5_ds1_stab_attn]);
violinplot(null_stab_ds1_attn',{},'ViolinColor',[0.4 0 0],'ViolinAlpha',0.3,'ShowData',false);hold on  
xlim(xlims);
xticks([]);
ylim([.2 .5]);
for m = 1:length(lesioned_stab_ds1_attn)
    plot([m-.4 m+.4],[lesioned_stab_ds1_attn(m) lesioned_stab_ds1_attn(m)],'LineWidth',4,'Color',[0.4 0 0])

end
scatter([1:length(sig_ds1_stabanat_attn)], [sig_ds1_stabanat_attn * .5],200,'k*','LineWidth',2);
set(gca,'YLim',ds1_SA_ylim,'YTick',ds1_SA_ylim(1):0.1:ds1_SA_ylim(2),'FontSize',25)


ds2_SA_ylim = [.35 .5];
ds2_WM_ylim = [.2 .35];

% Stability anatomy DS2 sustained attenion
subplot(3,2,3);
yline([fig5_ds2_stab_attn]);hold on
violinplot(null_stab_ds2_attn',{},'ViolinColor',[0.4 0 0],'ViolinAlpha',0.3,'ShowData',false);hold on  
xlim(xlims);
xticks([]);
for m = 1:length(lesioned_stab_ds2_attn)
    plot([m-.4 m+.4],[lesioned_stab_ds2_attn(m) lesioned_stab_ds2_attn(m)],'LineWidth',4,'Color',[0.4 0 0])

end
scatter([1:length(sig_ds2_stabanat_attn)], [sig_ds2_stabanat_attn * ds2_SA_ylim(2)],200,'k*','LineWidth',2);
set(gca,'YLim',ds2_SA_ylim,'YTick',ds2_SA_ylim(1):0.05:ds2_SA_ylim(2),'FontSize',25)


% Stability anatomy DS2 working memory
subplot(3,2,4);
yline([fig5_ds2_stab_wm]);hold on
violinplot(null_stab_ds2_wm',{},'ViolinColor',[0.4 0 0],'ViolinAlpha',0.3,'ShowData',false);hold on  
xlim(xlims);
xticks([]);
for m = 1:length(lesioned_stab_ds2_wm)
    plot([m-.4 m+.4],[lesioned_stab_ds2_wm(m) lesioned_stab_ds2_wm(m)],'LineWidth',4,'Color',[0.4 0 0])

end
scatter([1:length(sig_ds2_stabanat_wm)], [sig_ds2_stabanat_wm * ds2_WM_ylim(2)],200,'k*','LineWidth',2);
set(gca,'YLim',ds2_WM_ylim,'YTick',ds2_WM_ylim(1):0.05:ds2_WM_ylim(2),'FontSize',25)


ds3_SA_ylim = [0.12 0.22];
ds3_WM_ylim = [0.15 0.25];

% Stability anatomy DS3 sustained attenion
subplot(3,2,5);
yline([fig5_ds3_stab_attn]);hold on
violinplot(null_stab_ds3_attn',{},'ViolinColor',[0.4 0 0],'ViolinAlpha',0.3,'ShowData',false);hold on  
xlim(xlims);
xticklabels(network_names);
xtickangle(45);
yticks([1:3]);
%ylim(ds3_SA_ylim);
set(gca,'YLim',ds3_SA_ylim,'YTick',ds3_SA_ylim(1):0.05:ds3_SA_ylim(2),'FontSize',25)
for m = 1:length(lesioned_stab_ds3_attn)
    plot([m-.4 m+.4],[lesioned_stab_ds3_attn(m) lesioned_stab_ds3_attn(m)],'LineWidth',4,'Color',[0.4 0 0])

end
scatter([1:length(sig_ds3_stabanat_attn)], [sig_ds3_stabanat_attn * ds3_SA_ylim(2)],200,'k*','LineWidth',2);

% Stability anatomy DS3 working memory
subplot(3,2,6);
yline([fig5_ds3_stab_wm]);hold on
violinplot(null_stab_ds3_wm',{},'ViolinColor',[0.4 0 0],'ViolinAlpha',0.3,'ShowData',false);hold on  
xlim(xlims);
xticklabels(network_names);
xtickangle(45);
%ylim(ds3_WM_ylim);
for m = 1:length(lesioned_stab_ds3_wm)
    plot([m-.4 m+.4],[lesioned_stab_ds3_wm(m) lesioned_stab_ds3_wm(m)],'LineWidth',4,'Color',[0.4 0 0])

end
set(gca,'YLim',ds3_WM_ylim,'YTick',ds3_WM_ylim(1):0.05:ds3_WM_ylim(2),'FontSize',25)
scatter([1:length(sig_ds3_stabanat_wm)], [sig_ds3_stabanat_wm * ds3_WM_ylim(2)],200,'k*','LineWidth',2);

set(gcf, 'Position', [0 0 1800 1000]);
annotation('textbox', [0.2, 0.96, 0.6, 0.05], 'string', {['Sustained Attention']},'FontSize',30,'FontWeight','bold','EdgeColor','none')
annotation('textbox', [0.7, 0.96, 0.6, 0.05], 'string', {['Working Memory']},'FontSize',30,'FontWeight','bold','EdgeColor','none')

% ----------------------------------------------------------------------------------------------
% ------------------    PLOT FIGURE 6 Typicality Network Anatomy     --------------------
% ----------------------------------------------------------------------------------------------
%%


for n = 1:length(lesioned_typ_ds1_attn)

    % DS1 attention typicality

    if (min([sum(null_typ_ds1_attn(n,:) < lesioned_typ_ds1_attn(n)),   sum(null_typ_ds1_attn(n,:) > lesioned_typ_ds1_attn(n))]) + 1)/(size(null_typ_ds1_attn,2) + 1) < .025;
        sig_ds1_typanat_attn(n) = 1; % significant
    else
        sig_ds1_typanat_attn(n) = NaN; % not significant
    end


    % DS2 attention typicality
    if (min([sum(null_typ_ds2_attn(n,:) < lesioned_typ_ds2_attn(n)),   sum(null_typ_ds2_attn(n,:) > lesioned_typ_ds2_attn(n))]) + 1)/(size(null_typ_ds2_attn,2) + 1) < .025;
        sig_ds2_typanat_attn(n) = 1; % significant
    else
        sig_ds2_typanat_attn(n) = NaN; % not significant
    end

    % DS2 working mem typicality
    if (min([sum(null_typ_ds2_wm(n,:) < lesioned_typ_ds2_wm(n)),   sum(null_typ_ds2_wm(n,:) > lesioned_typ_ds2_wm(n))]) + 1)/(size(null_typ_ds2_wm,2) + 1) < .025;
        sig_ds2_typanat_wm(n) = 1; % significant
    else
        sig_ds2_typanat_wm(n) = NaN; % not significant
    end


    % DS3 attention typicality
    if (min([sum(null_typ_ds3_attn(n,:) < lesioned_typ_ds3_attn(n)),   sum(null_typ_ds3_attn(n,:) > lesioned_typ_ds3_attn(n))]) + 1)/(size(null_typ_ds3_attn,2) + 1) < .025;
        sig_ds3_typanat_attn(n) = 1; % significant
    else
        sig_ds3_typanat_attn(n) = NaN; % not significant
    end

    % DS3 working mem typicality
    if (min([sum(null_typ_ds3_wm(n,:) < lesioned_typ_ds3_wm(n)),   sum(null_typ_ds3_wm(n,:) > lesioned_typ_ds3_wm(n))]) + 1)/(size(null_typ_ds3_wm,2) + 1) < .025;
        sig_ds3_typanat_wm(n) = 1; % significant
    else
        sig_ds3_typanat_wm(n) = NaN; % not significant
    end



end % for n = 

%%

% Create violin plots for null values
xlims = [0 11];
ds1_typSA_ylim = [.3 .6];

figure
subplot(3,2,1);
yline([fig5_ds1_typ_attn]);
violinplot(null_typ_ds1_attn',{},'ViolinColor',[0 0.1 0.4],'ViolinAlpha',0.3,'ShowData',false);hold on  
xlim(xlims);
xticklabels([]);
for m = 1:length(lesioned_typ_ds1_attn)
    plot([m-.4 m+.4],[lesioned_typ_ds1_attn(m) lesioned_typ_ds1_attn(m)],'LineWidth',4,'Color',[0 0.1 0.4])

end
set(gca,'YLim',ds1_typSA_ylim,'YTick',ds1_typSA_ylim(1):0.1:ds1_typSA_ylim(2),'FontSize',25)
scatter([1:length(sig_ds1_typanat_attn)], [sig_ds1_typanat_attn * ds1_typSA_ylim(2)],200,'k*','LineWidth',2);


ds2_typSA_ylim = [-0.05 .21];
ds2_typWM_ylim = [0.05 .31];

% Typicality anatomy DS2 sustained attenion
subplot(3,2,3);
yline([fig5_ds2_typ_attn]);hold on
violinplot(null_typ_ds2_attn',{},'ViolinColor',[0 0.1 0.4],'ViolinAlpha',0.3,'ShowData',false);hold on  
xlim(xlims);
xticklabels([]);
for m = 1:length(lesioned_typ_ds2_attn)
    plot([m-.4 m+.4],[lesioned_typ_ds2_attn(m) lesioned_typ_ds2_attn(m)],'LineWidth',4,'Color',[0 0.1 0.4]);

end
set(gca,'YLim',ds2_typSA_ylim,'YTick',ds2_typSA_ylim(1):0.1:ds2_typSA_ylim(2),'FontSize',25)
scatter([1:length(sig_ds2_typanat_attn)], [sig_ds2_typanat_attn * ds2_typSA_ylim(2)],200,'k*','LineWidth',2);


% Typicality anatomy DS2 working memory
subplot(3,2,4);
yline([fig5_ds2_typ_wm]);hold on
violinplot(null_typ_ds2_wm',{},'ViolinColor',[0 0.1 0.4],'ViolinAlpha',0.3,'ShowData',false);hold on  
xlim(xlims);
xticklabels([]);
for m = 1:length(lesioned_typ_ds2_wm)
    plot([m-.4 m+.4],[lesioned_typ_ds2_wm(m) lesioned_typ_ds2_wm(m)],'LineWidth',4,'Color',[0 0.1 0.4]);

end
set(gca,'YLim',ds2_typWM_ylim,'YTick',ds2_typWM_ylim(1):0.1:ds2_typWM_ylim(2),'FontSize',25);
scatter([1:length(sig_ds2_typanat_wm)], [sig_ds2_typanat_wm * ds2_typWM_ylim(2)],200,'k*','LineWidth',2);

ds3_typSA_ylim = [.1 .24];
ds3_typWM_ylim = [.15 .29];

% Typicality anatomy DS3 sustained attenion
subplot(3,2,5);
yline([fig5_ds3_typ_attn]);hold on
violinplot(null_typ_ds3_attn',{},'ViolinColor',[0 0.1 0.4],'ViolinAlpha',0.3,'ShowData',false);hold on  
xlim(xlims);
xticklabels(network_names);
xtickangle(45);
for m = 1:length(lesioned_typ_ds3_attn)
    plot([m-.4 m+.4],[lesioned_typ_ds3_attn(m) lesioned_typ_ds3_attn(m)],'LineWidth',4,'Color',[0 0.1 0.4]);

end
set(gca,'YLim',ds3_typSA_ylim,'YTick',ds3_typSA_ylim(1):0.05:ds3_typSA_ylim(2),'FontSize',25)
scatter([1:length(sig_ds3_typanat_attn)], [sig_ds3_typanat_attn * ds3_typSA_ylim(2)],200,'k*','LineWidth',2);

% Typicality anatomy DS3 working memory
subplot(3,2,6);
yline([fig5_ds3_typ_wm]);hold on
violinplot(null_typ_ds3_wm',{},'ViolinColor',[0 0.1 0.4],'ViolinAlpha',0.3,'ShowData',false);hold on  
xlim(xlims);
xticklabels(network_names);
xtickangle(45);
for m = 1:length(lesioned_typ_ds3_wm)
    plot([m-.4 m+.4],[lesioned_typ_ds3_wm(m) lesioned_typ_ds3_wm(m)],'LineWidth',4,'Color',[0 0.1 0.4]);

end
set(gca,'YLim',ds3_typWM_ylim,'YTick',ds3_typWM_ylim(1):0.05:ds3_typWM_ylim(2),'FontSize',25)
scatter([1:length(sig_ds3_typanat_wm)], [sig_ds3_typanat_wm * ds3_typWM_ylim(2)],200,'k*','LineWidth',2);

set(gcf, 'Position', [0 0 1800 1000]);
annotation('textbox', [0.2, 0.96, 0.6, 0.05], 'string', {['Sustained Attention']},'FontSize',30,'FontWeight','bold','EdgeColor','none')
annotation('textbox', [0.7, 0.96, 0.6, 0.05], 'string', {['Working Memory']},'FontSize',30,'FontWeight','bold','EdgeColor','none')


% ----------------------------------------------------------------------------------------------
% -----------    PLOT FIGURE S1 Within-run Typicality relationship with performance  --------------------
% ---------------------------------------------------------------------------------------------
%% Figure S1

Chun_run_names = {'gradCPT','VSTM','MOT','movie','rest'};
HCP_run_names = {'0-back','2-back','relational','motor','language','social','emotion','gambling','rest'};


% Are values significant?
% DS1 Typicality
for s = 1:length(figSup1_ds1_attn_sig)
    if figSup1_ds1_attn_sig(s) < 0.05
        sup1_sig_ds1_attn(s) = 1;
    else
        sup1_sig_ds1_attn(s) = NaN;
    end
end

% DS2 Attention Typicality
for s = 1:length(figSup1_ds2_attn_sig)
    if figSup1_ds2_attn_sig(s) < 0.05
        sup1_sig_ds2_attn(s) = 1;
    else
        sup1_sig_ds2_attn(s) = NaN;
    end
end

% DS2 WM Typicality
for s = 1:length(figSup1_ds2_wm_sig)
    if figSup1_ds2_wm_sig(s) < 0.05
        sup1_sig_ds2_wm(s) = 1;
    else
        sup1_sig_ds2_wm(s) = NaN;
    end
end


% DS3 Attention Typicality
for s = 1:length(figSup1_ds3_attn_sig)
    if figSup1_ds3_attn_sig(s) < 0.05
        sup1_sig_ds3_attn(s) = 1;
    else
        sup1_sig_ds3_attn(s) = NaN;
    end
end

% DS3 WM Typicality
for s = 1:length(figSup1_ds3_wm_sig)
    if figSup1_ds3_wm_sig(s) < 0.05
        sup1_sig_ds3_wm(s) = 1;
    else
        sup1_sig_ds3_wm(s) = NaN;
    end
end


figure
% DS1 within-run typ
subplot(3,2,1);
bar(figSup1_ds1_attn,'FaceColor',[0 0.3 .7]);hold on
scatter([1:length(sup1_sig_ds1_attn)], [sup1_sig_ds1_attn * .6],100,'k*','LineWidth',2);
xticks([1:length(figSup1_ds1_attn)]);
xticklabels({'gradCPT','rest'});
ylabel("Spearman's rho");
xtickangle(45);
ylim([-.5 0.7])
set(gca, 'FontSize',25);

% DS2 within-run typ with attention
subplot(3,2,3);
bar(figSup1_ds2_attn, 'FaceColor',[0 0.3 .7]);hold on
scatter([1:length(sup1_sig_ds2_attn)], [sup1_sig_ds2_attn * .25],100,'k*','LineWidth',2);
xticks([1:length(figSup1_ds2_attn)]);
xticklabels(Chun_run_names);
ylabel("Spearman's rho");
xtickangle(45);
ylim([-.15 0.35])
set(gca, 'FontSize',25);

% dS2 within-run typ with WM
subplot(3,2,4);
bar(figSup1_ds2_wm, 'FaceColor',[0 0.3 .7]);hold on
scatter([1:length(sup1_sig_ds2_wm)], [sup1_sig_ds2_wm * .25],100,'k*','LineWidth',2);
% title('Dataset 2 typicality working mem');
xticks([1:length(figSup1_ds2_wm)]);
xticklabels(Chun_run_names);
xtickangle(45);
ylim([-.15 0.35])
ylabel("Spearman's rho");
set(gca, 'FontSize',25);

% dS3 within-run typ with WM
subplot(3,2,5);
bar(figSup1_ds3_attn, 'FaceColor',[0 0.3 .7]);hold on
scatter([1:length(sup1_sig_ds3_attn)], [sup1_sig_ds3_attn * .4],100,'k*','LineWidth',2);
xticks([1:length(figSup1_ds3_attn)]);
xticklabels(HCP_run_names);
ylabel("Spearman's rho");
ylim([-.15 0.5])
set(gca, 'FontSize',25);
xtickangle(45);

subplot(3,2,6);
bar(figSup1_ds3_wm, 'FaceColor',[0 0.3 .7]);hold on
scatter([1:length(sup1_sig_ds3_wm)], [sup1_sig_ds3_wm * .4],100,'k*','LineWidth',2);
% title('Dataset 3 typicality working mem');
xticks([1:length(figSup1_ds3_wm)]);
xticklabels(HCP_run_names);
ylabel("Spearman's rho");
xtickangle(45);
ylim([-.15 0.5])
set(gca, 'FontSize',25);
set(gcf, 'Position',[100 100 1200 1000]);

annotation('textbox', [0.15, 0.935, 0.46, 0.05], 'string', 'A. Sustained Attention','FontSize',40,'EdgeColor','none')
annotation('textbox', [0.58, 0.935, 0.45, 0.05], 'string', 'B. Working Memory','FontSize',40,'EdgeColor','none')

% Label datasets
annotation('textarrow', [0.026, 0.3],[0.9, 0.3], 'String', {['Dataset 1']},'FontSize',35, 'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)
annotation('textarrow', [0.03, 0.3], [0.6, 0.3], 'String', {['Dataset 2']},'FontSize',35,'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)
annotation('textarrow', [0.02, 0.3], [0.3, 0.3], 'String', {['Dataset 3']},'FontSize',35,'HeadStyle', 'none', 'LineStyle', 'none','TextRotation',90)


%% Figure S2
% ----------------------------------------------------------------------------------------------
% -----------    PLOT FIGURE S2 Typicality performance quartiles  --------------------
% ---------------------------------------------------------------------------------------------

sa_typ_mat(1,:) = ds1_attn_quants;
sa_typ_mat(2,:) = ds2_attn_quants;
sa_typ_mat(3,:) = ds3_attn_quants;
err_sa(1,:) = ds1_attn_err;
err_sa(2,:) = ds2_attn_err;
err_sa(3,:) = ds3_attn_err;

wm_typ_mat(1,:) = NaN([1,length(ds2_wm_quants)]);
wm_typ_mat(2,:) = ds2_wm_quants;
wm_typ_mat(3,:) = ds3_wm_quants;

err_wm(1,:) = NaN([1,length(ds2_wm_quants)]);
err_wm(2,:) = ds2_wm_err;
err_wm(3,:) = ds3_wm_err;

%x = 1:length(ds1_attn_quants);

figure;

subplot(1,2,1);
for r = 1:size(sa_typ_mat,1)
plot([1:size(sa_typ_mat,2)], sa_typ_mat(r,:),'LineWidth',6,'DisplayName',['Dataset ' num2str(r)]);hold on
errorbar(sa_typ_mat(r,:),err_sa(r,:),'Color','k','DisplayName','');
end
xlim([0.5 4.5]);
ylim([.55 .77]);
xlabel('Performance quartile');
xticklabels({'low','','','high'})
ylabel('Mean Typicality')
set(gca,'FontSize',25);
legend
title('Sustained Attention')

subplot(1,2,2);
for r2 = 1:size(wm_typ_mat,1)
plot([1:size(wm_typ_mat,2)], wm_typ_mat(r2,:),'LineWidth',6,'DisplayName',['Dataset ' num2str(r2)]);hold on
errorbar(wm_typ_mat(r2,:),err_wm(r2,:),'Color','k','DisplayName','');
end
xlim([0.5 4.5]);
ylim([.55 .77]);
xticklabels({'low','','','high'})
 xlabel('Performance quartile');
 ylabel('Mean Typicality')
 set(gca,'FontSize',25)
 set(gcf,'Position',[100 100 1000 400]);
 legend
title('Working Memory')




%% Figure S3
% ----------------------------------------------------------------------------------------------
% -----------    PLOT FIGURE S4 Within-connectome standard deviation   --------------------
% ---------------------------------------------------------------------------------------------


    load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS1_data/gradCPT_flatmats_motion_behav_data']);
for r = 1:size(fc_flat,2)
    for p = 1:size(fc_flat{1},1)
        if sum(isnan(fc_flat{r}(p,:))) == length(fc_flat{r}(p,:))
            figSup3_ds1(p,r) = NaN;
        else
        figSup3_ds1(p,r) = nanstd(fc_flat{r}(p,:));     
        end
    end
end

runs = [1 1 1 2 2];
for i = 1:max(runs)
log = find(runs == i);
    for l = 1:length(log)
        ds1_FC_pre(l,:) = real(nanmean(fc_flat{log(l)}));
        
    end
ds1_FC_pre(isinf(ds1_FC_pre)) = NaN;
ds1_FC_std(i) = nanstd(tanh(nanmean(ds1_FC_pre)));
%
clear ds1_FC_pre
end



    load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS2_data/Chun_flatmats_motion_behav']);

for r = 1:size(fc_flat,2)
    for p = 1:size(fc_flat{1},1)
                if sum(isnan(fc_flat{r}(p,:))) == length(fc_flat{r}(p,:))
            figSup3_ds2(p,r) = NaN;
        else
        figSup3_ds2(p,r) = nanstd(fc_flat{r}(p,:));       
                end
    end
end

runs = [1 1 2 2 3 3 4 4 5 5];
for i = 1:max(runs)
log = find(runs == i);
    for l = 1:length(log);
        ds2_FC_pre(l,:) = nanmean(fc_flat{log(l)});
        
    end
    ds2_FC_pre(isinf(ds2_FC_pre)) = NaN;
ds2_FC_std(i) = nanstd(tanh(nanmean(ds2_FC_pre)));
clear ds2_FC_pre
end

    load(['/Users/annacorriveau/Documents/GitHub/FCStability_Typicality/data/DS3_data/HCP_flatmats_motion_behav_data']);
for r = 1:size(fc_flat,2)
    for p = 1:size(fc_flat{1},1)
                if sum(isnan(fc_flat{r}(p,:))) == length(fc_flat{r}(p,:))
            figSup3_ds3(p,r) = NaN;
        else
        figSup3_ds3(p,r) = nanstd(fc_flat{r}(p,:));  
                end
    end
end

runs = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 9 9];
for i = 1:max(runs)
log = find(runs == i);
    for l = 1:length(log);
        ds3_FC_pre(l,:) = nanmean(fc_flat{log(l)});
        
    end
    ds3_FC_pre(isinf(ds3_FC_pre)) = NaN;
ds3_FC_std(i) = nanstd(tanh(nanmean(ds3_FC_pre)));
clear ds3_FC_pre
end

v1 = figure;
vs1 = subplot(1,3,1);hold on
vs1.Position = [0.04 0.05 0.13 0.9]; 
boxplot(figSup3_ds1,'PlotStyle','compact','Color','k');
hv = findobj(gca,'Tag','Box','Type','line');
hv(1).Color = [0.9 0.6 0];
hv(2).Color = [0.9 0.6 0];
hv(3).Color = [0.4 0 0.5];
hv(4).Color = [0.4 0 0.5];
hv(5).Color = [0.4 0 0.5];

hh = findobj(gca,'Tag','Line');
% boxplot(DS1_SD,'PlotStyle','compact','ColorGroup',ds1_group);
xticks([1:5]);
xtickangle(45);
xticklabels({'gradCPT','gradCPT','gradCPT','rest','rest'});
ylabel('Connectome standard deviation')
title('Dataset 1');
set(gca,'FontSize',20);
ylim([ .1 .4])


vs2 = subplot(1,3,2);
vs2.Position = [.195 .05 .275 .9];
boxplot(figSup3_ds2,'PlotStyle','compact','Color','k');hold on
hv = findobj(gca,'Tag','Box');
hv(1).Color = [0.9 0.6 0];
hv(2).Color = [0.9 0.6 0];
hv(3).Color = [0.2 0.6 0.9];
hv(4).Color = [0.2 0.6 0.9];
hv(5).Color = [0.8 0.2 0.4];
hv(6).Color = [0.8 0.2 0.4];
hv(7).Color = [0.8 0.7 0];
hv(8).Color = [0.8 0.7 0];
hv(9).Color = [0.4 0 0.5];
hv(10).Color = [0.4 0 0.5];

xticks([1:10]);
xtickangle(45);
xticklabels({'gradCPT','gradCPT','VSTM','VSTM','MOT','MOT','movie','movie','rest','rest'});
title('Dataset 2');
set(gca,'FontSize',20);

vs3 = subplot(1,3,3);
vs3.Position = [.46 .05 .52 .9];
boxplot(figSup3_ds3,'PlotStyle','compact','Color','k');hold on
hv = findobj(gca,'Tag','Box');
hv(1).Color = [0.9 0.6 0];
hv(2).Color = [0.9 0.6 0];
hv(3).Color = [0.9 0.6 0];
hv(4).Color = [0.9 0.6 0];
hv(5).Color = [0.2 0 0.6];
hv(6).Color = [0.2 0 0.6];
hv(7).Color = [0.8 0.4 0.7];
hv(8).Color = [0.8 0.4 0.7];
hv(9).Color = [0.0 0.3 0.3];
hv(10).Color = [0.0 0.3 0.3];
hv(11).Color = [0.9 0.4 0.2];
hv(12).Color = [0.9 0.4 0.2];
hv(13).Color = [0.1 0.1 0.8];
hv(14).Color = [0.1 0.1 0.8];
hv(15).Color = [ 0 0.7 0.5];
hv(16).Color = [ 0 0.7 0.5];
hv(17).Color = [0.8 0.7 0];
hv(18).Color = [0.8 0.7 0];
hv(19).Color = [0.4 0 0.5];
hv(20).Color = [0.4 0 0.5];

xticks([1:20]);
xtickangle(45);
xticklabels({'0-back','0-back','2-back','2-back','relational','relational','motor','motor','language','language','social','social','emotional','emotional','gambling','gambling','rest','rest','rest','rest'});
set(gca,'FontSize',20);
title('Dataset 3');
set(v1,'units','points','position',[0 50 1800 600])



% AND S5
%%
h1col = [0.4 0 0.5;  0.9 0.6 0];
v2 = figure
vs1 = subplot(1,3,1);hold on
for i = 1:length(ds1_FC_std)
gg = bar(i,ds1_FC_std(i),'FaceColor','flat');hold on
gg.CData = h1col(i,:);
end
% gg.CData(2,:) = h1col(1,:);
% gg.CData(1,:) = h1col(2,:);
xticks([1 2]);
xticklabels({'gradCPT','rest'});
xtickangle(45);
ylim([0 .18]);
ylabel({['Average task connectome standard deviation']});
title('Dataset 1');
set(gca,'FontSize',20);hold off
vs1.Position = [0.05 0.2 0.14 0.7];

h2col = [0.4 0 0.5; 0.8 0.7 0; 0.8 0.2 0.4; 0.2 0.6 0.9;  0.9 0.6 0];

vs2 = subplot(1,3,2);
vs2.Position = [.23 .2 .28 .7];
for ii = 1:length(ds2_FC_std)
h2 = bar(ii,ds2_FC_std(ii),'FaceColor','flat');hold on
h2.CData = h2col(ii,:);
end
xticks([1:5]);
xticklabels({'gradCPT','VSTM','MOT','movie','rest'});
xtickangle(45);
ylim([0 .18]);
title('Dataset 2');
set(gca,'FontSize',20);


h3col = [0.4 0 0.5; 0.8 0.7 0; 0 0.7 0.5; 0.1 0.1 0.8 ; 0.9 0.4 0.2; 0.0 0.3 0.3; 0.8 0.4 0.7; 0.2 0 0.6;  0.9 0.6 0];


vs3 = subplot(1,3,3);
vs3.Position = [.545 .2 .45 .7];
for iii = 1:length(ds3_FC_std)
h3 = bar(iii,ds3_FC_std(iii),'FaceColor','flat');hold on
h3.CData = h3col(iii,:);
end
xticks([1:9]);
xticklabels({'0-back','2-back','relational','motor','language','social','emotional','gambling','rest'});
xtickangle(45);
ylim([0 .18]);
title('Dataset 3');
set(gca,'FontSize',20);


set(v2,'units','points','position',[0 50 1800 600])





