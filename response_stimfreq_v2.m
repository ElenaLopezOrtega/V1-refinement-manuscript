function [mean_PYY_f, mean_PYY_f_resp, dff_cells_oris_resp,idx, mean_PYY_cnt, mean_PYY_cnt_resp] = response_stimfreq_v2 (PYY,PYY_Gr, f, f_cnt, dff_cells_oris_pref_ori)

%Calculate averaged responses at a particular frequency
%Select responsive ROIS: ROIS which response to stimulation is higher
%than max response that happens by chance when exposed to gray screen
%Calculate responsive ROIs' averaged response to a particular frequency

PYY_cnt = PYY(f_cnt,:); % signal at control frequency (usually 1Hz, cells respond to this frequency during normal stimulation)
mean_PYY_cnt = mean (PYY_cnt);

PYY_f = PYY(f,:); %sigal at frequency of interest
mean_PYY_f = mean (PYY_f);

PYY_G_f = PYY_Gr(f,:);
PYY_G_f_mean = mean (PYY_G_f);
PYY_G_f_stdv = std (PYY_G_f);
threshold = PYY_G_f_mean + 3*PYY_G_f_stdv; %Z-score of 3
%  max_PYY_G_f = max (PYY_G_f); % thershold for signal at frequency of interest

idx = find (PYY_f > threshold); %select responsive cells at frequence of interest
% idx = find (PYY_f > max_PYY_G_f); %select responsive cells at frequence of interest
PYY_f_resp = PYY_f (idx);
mean_PYY_f_resp = mean (PYY_f_resp); %response of responsive cells at the frequence of interest

PYY_cnt_resp = PYY_cnt (idx);
mean_PYY_cnt_resp = mean (PYY_cnt_resp); % response of responsive cells at the control frequency

dff_cells_oris_resp = dff_cells_oris_pref_ori(:,idx);
