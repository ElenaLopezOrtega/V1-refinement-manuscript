%% Segments

% get values
VS = myPoolSRPCN.values.SegmentsSpines;
VS_cell_number = [];
BS = myPoolCNTCN.values.SegmentsSpines;
BS_cell_number = [];

%% long term changes (baseline 1, baseline 2, day 5 and day 7

sessions = [1 2 4 6];

[fVS] = long_term_changes (VS, sessions);
[fBS] = long_term_changes (BS, sessions);

% histogram
binRange = 0:0.05:2.5;
hcx = histcounts(fVS(:,3),[binRange Inf]);
hcy = histcounts(fBS(:,3),[binRange Inf]);
figure
hb = bar(binRange,[hcx;hcy]');
hb(1).FaceColor = '0.85 0.85 0.85';
hb(2).FaceColor = '0 0 0';


%% long term changes - segments per neuron

[VS_neurons, VS_neuron_order, VS_segm_neuron] = segments_neuron (VS, VS_cell_number, sessions);
[BS_neurons, BS_neuron_order, BS_segm_neuron] = segments_neuron (BS, BS_cell_number, sessions);

% histogram neurons
clear hcx
clear hcy
clear hb

binRange = 0:0.05:2.5;
hcx = histcounts(VS_neurons(:,3),[binRange Inf]);
hcy = histcounts(BS_neurons(:,3),[binRange Inf]);
figure
hb = bar(binRange,[hcx;hcy]');
hb(1).FaceColor = '0.75 0.75 0.75';
hb(2).FaceColor = '0 0 0';

%% Fraction of total segments that increase and decrease, and cells they belong to

BS_SD = std (fBS(:,3), 'omitnan');
BS_mean = mean (fBS (:,3), 'omitnan');
threshold_inc = BS_mean+BS_SD;
threshold_dec = BS_mean-BS_SD;

%visual stimulation
[VS_segm_prc, VS_neuron_inc, VS_neuron_dec] = increase_decrease (threshold_inc, threshold_dec, VS, sessions, VS_cell_number);
[BS_segm_prc, BS_neuron_inc, BS_neuron_dec] = increase_decrease (threshold_inc, threshold_dec, BS, sessions, BS_cell_number);

%% Spines

% get values
VS_spines = myPoolSRPCN.values.Spines;
BS_spines = myPoolCNTCN.values.Spines;

%% long term changes (baseline 1, baseline 2, day 5 and day 7

 sessions = [1 2 4 6];
% sessions = [1 2 4 5 7];

[fVS_spines] = long_term_changes (VS_spines, sessions);
[fBS_spines] = long_term_changes (BS_spines, sessions);

% Kolmogorov_Smirnov test for spine GluA1 distriburion on exp. day 5
fVS_spines_d5 = fVS_spines(:,3);
fBS_spines_d5 = fBS_spines(:,3);

[h,p] = kstest2(fVS_spines_d5, fBS_spines_d5)

% histogram
clear hcx
clear hcy
clear hb

binRange = 0:0.05:3;
hcx = histcounts(fVS_spines(:,3),[binRange Inf]);
hcy = histcounts(fBS_spines(:,3),[binRange Inf]);
figure
hb = bar(binRange,[hcx;hcy]');
hb(1).FaceColor = '0.85 0.85 0.85';
hb(2).FaceColor = '0 0 0';

%cumulative histogram
VS_higher = find (fVS_spines(:,3)>3);
fVS_spines(VS_higher,3)=3;

BS_higher = find (fBS_spines(:,3)>3);
fBS_spines(BS_higher,3)=3;

figure;
[h1,stats1] = cdfplot(fVS_spines (:,3))
hold on
[h2,stats2] = cdfplot(fBS_spines (:,3))


%fraction of total spines
BS_SD_spines = std (fBS_spines(:,3), 'omitnan');
BS_mean_spines = mean (fBS_spines (:,3), 'omitnan');
threshold_spines_inc = BS_mean_spines + BS_SD_spines;
threshold_spines_dec = BS_mean_spines - BS_SD_spines;

[VS_prc_spines,VS_spines_n] = increase_decrease_spines(threshold_spines_inc, threshold_spines_dec, fVS_spines, sessions);
[BS_prc_spines,BS_spines_n] = increase_decrease_spines(threshold_spines_inc, threshold_spines_dec, fBS_spines, sessions);

% fVS_spines_inc = find (fVS_spines (:, 3)>threshold_spines_inc);
% VS_prc_spines (1) = size (fVS_spines_inc, 1)/ size (fVS_spines,1);
% 
% fVS_spines_dec = find (fVS_spines (:, 3)<threshold_spines_dec);
% VS_prc_spines (2) = size (fVS_spines_dec, 1)/ size (fVS_spines,1);
% 
% fBS_spines_inc = find (fBS_spines (:, 3)>threshold_spines_inc);
% BS_prc_spines (1) = size (fBS_spines_inc, 1)/ size (fBS_spines,1);
% 
% fBS_spines_dec = find (fBS_spines (:, 3)<threshold_spines_dec);
% BS_prc_spines (2) = size (fBS_spines_dec, 1)/ size (fBS_spines,1);

%% spines per segment

for i = 1:size (VS,1)
prompt = 'Which segment do you want to analyzed? ';
VS_spseg = input(prompt);

[VS_prc_spseg,VS_spseg_n] = increase_decrease_spines(threshold_spines_inc, threshold_spines_dec, VS_spseg, sessions);

VS_prc_spseg_total(i,:)= VS_prc_spseg;
VS_spseg_n_total (i) = VS_spseg_n;

end

VS_prc_spseg_total(:,3)= 1 - VS_prc_spseg_total (:,1) - VS_prc_spseg_total (:,2) ;
VS_sum_spines = sum(VS_spseg_n_total);


for i = 1:size (BS,1)
prompt = 'Which segment do you want to analyzed? ';
BS_spseg = input(prompt);

[BS_prc_spseg,BS_spseg_n] = increase_decrease_spines(threshold_spines_inc, threshold_spines_dec, BS_spseg, sessions);

BS_prc_spseg_total(i,:)= BS_prc_spseg;
BS_spseg_n_total (i) = BS_spseg_n;

end

BS_prc_spseg_total(:,3)= 1 - BS_prc_spseg_total (:,1) - BS_prc_spseg_total (:,2) ;
BS_sum_spines = sum(BS_spseg_n_total);

%spine selection
% VS_spseg = VS_spseg(:, [1 2 4 6]);
% VS_spseg(isnan(VS_spseg)) = 0;
% VS_spseg_mean = mean (VS_spseg, 2);
% del = VS_spseg_mean ==0;
% VS_spseg(del, :) = [];
% 
% VS_spseg_mean_exp = mean(VS_spseg(:,[3 4]),2);
% del2 = VS_spseg_mean_exp ==0;
% VS_spseg(del2, :) = [];
% 
% %spine distribution per segments
% VS_spseg_inc = find (VS_spseg (:, 3)>threshold_spines_inc);
% VS_prc_spseg (1) = size (VS_spseg_inc, 1)/ size (VS_spseg,1);
% 
% VS_spseg_dec = find (VS_spseg (:, 3)<threshold_spines_dec);
% VS_prc_spseg (2) = size (VS_spseg_dec, 1)/ size (VS_spseg,1);
% 
% BS_spseg_inc = find (BS_spseg (:, 3)>threshold_spines_inc);
% BS_prc_spseg (1) = size (BS_spseg_inc, 1)/ size (VS_spseg,1);
% 
% BS_spseg_dec = find (BS_spseg (:, 3)<threshold_spines_dec);
% BS_prc_spseg (2) = size (BS_spseg_dec, 1)/ size (BS_spseg,1);


