
%% Retrieve info for each ROI pair and concatenate all mice
persist_S0_T = [];
persist_S1_T = [];
persist_oris_T = [];
trans_T=[];
trans_oris_T=[];
newadd_T=[];
newadd_oris_T=[];

%repeat for each mouse
persist_oris (1,:) = archive_long {4,1};
persist_oris (2,:) = archive_long {4,2};

persist_oris_T = horzcat (persist_oris_T, persist_oris);

trans_oris = archive_nolong{4,1};
trans_oris_T = horzcat (trans_oris_T, trans_oris);

newadd_oris = archive_nolong{4,2};
newadd_oris_T = horzcat (newadd_oris_T, newadd_oris);

persist_S0 = archive_long {3,1};
persist_S0_T = horzcat (persist_S0_T, persist_S0);

persist_S1 = archive_long {3,2};
persist_S1_T = horzcat (persist_S1_T, persist_S1);

trans = archive_nolong {3,1};
trans_T = horzcat (trans_T, trans);

newadd = archive_nolong {3,2};
newadd_T = horzcat (newadd_T, newadd);

clear persist_oris trans_oris newadd_oris persist_S0 persist_S1 trans newadd

%% Dynamics
dyn_summary = [];

perc_pers = size(persist_S0_T,2)/(size(persist_S0_T,2) + size(trans_T,2));
perc_trans = size(trans_T,2)/(size(persist_S0_T,2) + size(trans_T,2));
perc_newadd = size(newadd_T,2)/(size(persist_S0_T,2) + size(trans_T,2));

dyn_summary (1,1)= perc_trans;
dyn_summary (1,2)= perc_newadd;
dyn_summary (1,3)= perc_pers;

clear perc_trans perc_newadd perc_pers
%% Study difference in between populations

%plot
[s14_pers_S0,s14_pers_S0_T]  = plot_SEM (persist_S0_T);
[s14_pers_S1,s14_pers_S1_T]  = plot_SEM (persist_S1_T);
[s14_trans,s14_trans_T]  = plot_SEM (trans_T);
[s14_newadd,s14_newadd_T]  = plot_SEM (newadd_T);

s14_summary (1) = s14_trans_T;
s14_summary (2) = s14_pers_S0_T;
s14_summary (3) = s14_pers_S1_T;
s14_summary (4) = s14_newadd_T;

s14_trans (size (s14_trans,1)+1:size (s14_pers_S0,1))=NaN;
s14_pers_S0 (size (s14_pers_S0,1)+1:size (s14_trans,1))=NaN;
S0_percell = horzcat (s14_trans, s14_pers_S0);
figure()
v = violinplot(S0_percell);
ylim ([0 8])
% 
% s14_add_S1 (size (s14_add_S1,1)+1:size (s14_surv_S1,1))=NaN;
% S1_percell = horzcat (s14_add_S1, s14_surv_S1);
% figure()
% v = violinplot(S1_percell);
% ylim ([0 8])
% 
% Pers_percell = horzcat (s14_pers_S0, s14_pers_S1);
% figure()
% v = violinplot(Pers_percell);
% ylim ([0 8])

%% synchrony of persistent neurons
map = [];

% S0
[r2,lags2] = xcorr(persist_S0_T);

max_r2 = max(r2);
norm_r2 = r2./max_r2;

syn2_S0= norm_r2(3360,:);
figure()
histogram (syn2_S0, 'Normalization','probability');
ylim ([0 0.3]);

figure()
cdfplot(syn2_S0)

syn2_S0_m = reshape(syn2_S0,size(persist_S0_T,2),[]);

figure()
h_S0=heatmap (syn2_S0_m, 'Colormap',map);

norm_r2_subset = norm_r2(3210:3510,:);
mean_subset = mean (norm_r2_subset,2);
x = 1:size (mean_subset,1);
figure(); plot (x, mean_subset);
ylim ([0.5 0.8])

clear max_r2 norm_r2
clear r2 lags2
clear norm_r2_subset mean_subset x

med_S0 = median (syn2_S0);
Q_S0 = quantile (syn2_S0,[0.25 0.75]);

% S1
[r2,lags2] = xcorr(persist_S1_T);

max_r2 = max(r2);
norm_r2 = r2./max_r2;

syn2_S1= norm_r2(3360,:);
figure()
histogram (syn2_S1, 'Normalization','probability');
ylim ([0 0.3]);

figure()
cdfplot(syn2_S1)

syn2_S1_m = reshape(syn2_S1,size(persist_S1_T,2),[]);

figure()
h_S1=heatmap (syn2_S1_m, 'Colormap',map);

norm_r2_subset = norm_r2(3210:3510,:);
mean_subset = mean (norm_r2_subset,2);
x = 1:size (mean_subset,1);
figure(); plot (x, mean_subset);
ylim ([0.5 0.8])

clear max_r2 norm_r2
clear r2 lags2
clear norm_r2_subset mean_subset x

med_S1 = median (syn2_S1);
Q_S1 = quantile (syn2_S1,[0.25 0.75]);

%Statistic test

[h,p] = kstest2(syn2_S0,syn2_S1);

figure ()
cdfplot (syn2_S0)
hold on
cdfplot (syn2_S1)

%% synchrony in initial response (30s after gratings exposure)

% S0
[r2,lags2] = xcorr(persist_S0_T(304:1203,:));

max_r2 = max(r2);
norm_r2 = r2./max_r2;

syn2_S0= norm_r2(900,:);
figure()
histogram (syn2_S0, 'Normalization','probability');
ylim ([0 0.3]);

syn2_S0_m = reshape(syn2_S0,size(persist_S0_T,2),[]);

figure()
h_S0=heatmap (syn2_S0_m, 'Colormap',map);

norm_r2_subset = norm_r2(750:1050,:);
mean_subset = mean (norm_r2_subset,2);
x = 1:size (mean_subset,1);
figure(); plot (x, mean_subset);

clear max_r2 norm_r2
clear r2 lags2
clear norm_r2_subset mean_subset x

% S1
[r2,lags2] = xcorr(persist_S1_T(304:1203,:));

max_r2 = max(r2);
norm_r2 = r2./max_r2;

syn2_S1= norm_r2(900,:);
figure()
histogram (syn2_S1, 'Normalization','probability');
ylim ([0 0.3]);

syn2_S1_m = reshape(syn2_S1,size(persist_S1_T,2),[]);

figure()
h_S1=heatmap (syn2_S1_m, 'Colormap',map);

norm_r2_subset = norm_r2(750:1050,:);
mean_subset = mean (norm_r2_subset,2);
x = 1:size (mean_subset,1);
figure(); plot (x, mean_subset);

clear max_r2 norm_r2
clear r2 lags2
clear norm_r2_subset mean_subset x

[h,p] = kstest2(syn2_S0,syn2_S1);

%% calcium transients

[~, ~, ~] =  select_Ca_transients (persist_S0_T);