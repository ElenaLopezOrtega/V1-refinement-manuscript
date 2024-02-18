% Code for analysis of calcium timeseries data during SRP stimulation. 
% 10s baseline, 100s phase-reversing sinusoidal? gratings, reversal every
% 500ms 
%
% 

%% Studying invidual ROIs

dff_cells = [];
dff_cells = mfRegionTraces_denoised_df_f;


% dff = mfRegionTraces_denoised_df_f;
% dff_cells = vertcat(dff_cells, dff);

% dff_cells = dff_cells (:, 1:23520);

% plot all ROIs
% for i = 1:size (dff_cells,1)
%     figure()
%     plot (dff_cells(i,:));
% end 

% if removing cells is necessary
%idx = any(dff_cells' >= 100);
idx = any(abs (dff_cells') >= 15); % percentile 99 of max responses from all ROIs include in the analysis in 14.27
col = find(idx);

if isempty(col)==0
    dff_cells (col,:) = NaN;
end   

clear idx col

%% Separate gray screen form orientations

dff_cells_G =[];
dff_cells_oris = [];

x='S0';
 
if x == 'S0'
   dff_cells_G = dff_cells(:, 1:3360);
   dff_cells_oris = dff_cells (:, 3361:size(dff_cells,2));
   e=1;

elseif x == 'S1' 
    dff_cells_oris = dff_cells (:, 1:20160);
    dff_cells_G = dff_cells(:,20161:size(dff_cells,2));
    e=2;
else
    dff_cells_oris = dff_cells (:, 1:20160);
    dff_cells_G = dff_cells(:,20161:size(dff_cells,2));
    e=3;
end
   
dff_cells_oris = dff_cells_oris';
dff_cells_G = dff_cells_G';

archive_dff {1,e} = dff_cells_oris;
archive_dff {2,e} = dff_cells_G;

clear x 

%% Study activity - Active ROIs with responses over 3ZScores

%% Select Ca transients (define as local peaks above 3SD of the population)  
clear summary_frames
clear summary_cells
clear Is_active

[summary_frames, summary_cells, Is_active] = select_Ca_transients (dff_cells_oris);

% Variables explanation:
% first colum: all ROIs
% second column: ROIs with more than 500 frames above threshold
% third column: ROIs with more than 650 frames above threshold
% forth column: ROIs with mopre than 800 frames above threshold
% summary frames: 1st raw: averaged number of frames, 2nd raw: average response
% amplitude, 3rd raw: averaged max response amplitude

archive_activity {1,e} = Is_active;
archive_activity {2,e} = vertcat(summary_cells, summary_frames);
%% Study responsiveness: Responsive cells should have a prefered orientation and present a 1Hz component on power analysis

%% Study orientation preference 
clear Is_ori_select
clear ori_pref ori_orth
clear dff_select_oris
clear percent_cells_select

[Is_ori_select, ori_pref, ori_orth, dff_select_oris] = pref_orientation (dff_cells_oris, Is_active);

%Plot oientation seletive ROIs
% for j = 1: size (dff_select_oris,2)
%     figure()
%     plot (dff_select_oris (:,j));
% end  

% Calculate percentage of cells that have a prefered orientation
percent_cells_select = sum (Is_ori_select)/sum (Is_active);

%% Select preferred and orthogonal orientation for orientation selective cells

clear oris orth dff_select_oris_2 dff_select_prefori dff_select_orthori

oris = ori_pref;
oris(isnan(oris))=0;
oris = nonzeros (oris);

orth= ori_orth;
orth(isnan(orth))=0;
orth = nonzeros (orth);

dff_select_oris_2 = reshape (dff_select_oris, 3360,6,[]);
% dff_oris_select_2 = dff_oris_select_2(304:3308,:,:);

for i = 1:size (dff_select_oris_2,3)
    dff_select_prefori(:,i) = dff_select_oris_2 (:,oris(i),i);
end

for i = 1:size (dff_select_oris_2,3)
    dff_select_orthori(:,i) = dff_select_oris_2 (:,orth(i),i);
end

archive_select {1,e} = Is_ori_select;
archive_select {2,e} = dff_select_oris;
archive_select {3,e} = ori_pref;
archive_select {4,e} = dff_select_prefori;
archive_select {5,e} = ori_orth;
archive_select {6,e} = dff_select_orthori;
archive_select {7,e} = percent_cells_select;

clear dff_select_oris_2 i oris orth
%% Power analysis of ROIs with a prefered orientation
clear PYY_pref PYY_orth PYY_G f
clear dff_select_prefori_grat dff_select_orthori_grat dff_select_G

%power analysis during grating stimulation

dff_select_prefori_grat = dff_select_prefori(304:3308,:); %select response during grating stimulation
dff_select_orthori_grat = dff_select_orthori(304:3308,:);

[PYY_pref,f] = power_analysis (dff_select_prefori_grat);
[PYY_orth] = power_analysis (dff_select_orthori_grat);

% PYY_mean = mean (PYY,2);
% figure; plot(f(51:1504),PYY_mean(51:1504))
% axis ([0.5 15 -0.5 10])
% ylabel('Power spectral density')
% xlabel('Frequency (Hz)')

% power to gray screen
dff_select_G = dff_cells_G (304:3308, Is_ori_select==1);

[PYY_G] = power_analysis (dff_select_G);

% PYY_norm = PYY./PYY_G;
% figure; h = heatmap (t,f(80:205),PYY_norm(80:205,:),'Colormap',jet);

% PYY_1Hz = PYY(98,:);
% PYY_1Hz_mean = mean (PYY_1Hz, 2);
% 
% PYY_XHz = PYY(121,:);
% PYY_G_1Hz = PYY_G(98,:);


clear dff_select_prefori_grat dff_select_orthori_grat dff_select_G

%% Power analysis of ROIs with a prefered orientation during ONSET
clear PYY_pref PYY_orth PYY_G f
clear dff_select_prefori_grat dff_select_orthori_grat dff_select_G

%If you are using mouse archives
e=1; %sessions 
clear dff_select_prefori dff_select_orthori Is_ori_select dff_cells_G Is_active ori_pref
dff_select_prefori = archive_select{4,e};
dff_select_orthori = archive_select{6,e};
Is_ori_select = archive_select{1,e};
dff_cells_G = archive_dff{2,e};
Is_active = archive_activity {1,e};
ori_pref = archive_select {3,e};
dff_cells_oris = archive_dff {1,e};

%power analysis during grating stimulation

% dff_select_prefori_grat = dff_select_prefori(304:3308,:); %select response during grating stimulation
% dff_select_orthori_grat = dff_select_orthori(304:3308,:);
dff_select_prefori_grat = dff_select_prefori(304:754,:); %select response during grating stimulation
dff_select_orthori_grat = dff_select_orthori(304:754,:);

[PYY_pref,f] = power_analysis (dff_select_prefori_grat);
[PYY_orth] = power_analysis (dff_select_orthori_grat);

% PYY_mean = mean (PYY,2);
% figure; plot(f(51:1504),PYY_mean(51:1504))
% axis ([0.5 15 -0.5 10])
% ylabel('Power spectral density')
% xlabel('Frequency (Hz)')

% power to gray screen
% dff_select_G = dff_cells_G (304:3308, Is_ori_select==1);
dff_select_G = dff_cells_G (304:754, Is_ori_select==1);

[PYY_G] = power_analysis (dff_select_G);

% PYY_norm = PYY./PYY_G;
% figure; h = heatmap (t,f(80:205),PYY_norm(80:205,:),'Colormap',jet);

% PYY_1Hz = PYY(98,:);
% PYY_1Hz_mean = mean (PYY_1Hz, 2);
% 
% PYY_XHz = PYY(121,:);
% PYY_G_1Hz = PYY_G(98,:);

clear dff_select_prefori_grat dff_select_orthori_grat dff_select_G


 %% Select ROIs that respond at 1Hz
clear mean_PYY
clear mean_PYY_resp
clear dff_resp_G
clear idx_resp

f_int = 98; % position of ~1Hz in frequency vector
f1Hz = 98; % position of ~1Hz in frequency vector
 
[mean_PYY, mean_PYY_resp, dff_resp_prefori, dff_resp_orthori, idx_resp, ~, ~] = response_stimfreq_v3 (PYY_pref,PYY_G, f_int, f1Hz, dff_select_prefori, dff_select_orthori);

idx = find(Is_ori_select ==1);
idx2 = idx (idx_resp);
Is_resp = zeros (1, size (Is_ori_select,2));
Is_resp (idx2)=1;

percent_cells_resp = sum (Is_resp)/sum (Is_active);

for i = 1: size (dff_resp_prefori,2)
    figure()
    plot(dff_resp_prefori(:,i));
end    

% for i = 1: size (dff_resp_orthori,2)
%     figure()
%     plot(dff_resp_orthori(:,i));
% end    

mean_PYY_orth = mean(PYY_orth(f_int,:));
mean_PYY_orth_resp = mean(PYY_orth(f_int,idx_resp));

archive_power {1,e} = f;
archive_power {2,e} = PYY_pref;
archive_power {3,e} = mean_PYY;
archive_power {4,e} = PYY_orth;
archive_power {5,e} = mean_PYY_orth;

archive_resp {1,e} = Is_resp;
archive_resp {2,e} = dff_resp_prefori;
archive_resp {3,e} = mean_PYY_resp;
archive_resp {4,e} = mean_PYY_orth_resp;
archive_resp {5,e} = percent_cells_resp;
archive_resp {6,e} = dff_resp_orthori;

clear f1Hz f_int i idx idx2 idx_resp

sum(Is_active)
sum (Is_resp)

%% synchrony %Check variables!!!!!
% clear dff_resp_prefori_2 dff_resp_prefori_s
% clear r lags max_r norm_r syn
% clear r2 lags2 max_r2 norm_r2 syn2
% 
% dff_resp_prefori_2 = reshape (dff_resp_prefori, 30,[],size(dff_resp_prefori,2));
% dff_resp_prefori_s = squeeze (mean(dff_resp_prefori_2));
% [r,lags] = xcorr(dff_resp_prefori_s);
% 
% max_r = max(r);
% norm_r = r./max_r;
% 
% syn= norm_r(112,:);
% syn = reshape(syn,size(dff_resp_prefori_s,2),[]);
% 
% figure()
% h=heatmap (syn);
% 
% 
% [r2,lags2] = xcorr(dff_resp_prefori_s);
% 
% max_r2 = max(r2);
% norm_r2 = r2./max_r2;
% 
% syn2= norm_r2(3360,:);
% figure()
% histogram (syn2, 'Normalization','probability');
% 
% syn2 = reshape(syn2,size(dff_resp_prefori,2),[]);
% 
% figure()
% h=heatmap (syn2);
% 
% clear dff_resp_prefori_2 dff_resp_prefori_s
% clear max_r max_r2
%%

figure (); plot (mean(archive_resp{2,1},2))
ylim([0 2.5])
figure (); plot (mean(archive_resp{2,2},2))
ylim([0 2.5])
figure (); plot (mean(archive_resp{2,3},2))
ylim([0 2.5])
%% longitudinal
clear Is_ori_select_long ori_pref_long dff_cells_oris_long Is_resp_long dff_long_prefori

cells_long =[];

%add cells ID before proceding

for i = 1:size (cells_long,1)
Is_ori_select_long(i) = Is_ori_select(1,cells_long(i));
ori_pref_long (i) = ori_pref (1, cells_long(i));
dff_cells_oris_long (:,i)= dff_cells_oris (:,cells_long(i));
Is_resp_long(i) = Is_resp(1,cells_long(i));
end

dff_cells_oris_long = reshape (dff_cells_oris_long, 3360,6,[]);

for i = 1:size (cells_long,1)
    if isnan(ori_pref_long (i))
        dff_long_prefori(:,i) = nan (3360,1);
    else
        dff_long_prefori(:,i) = dff_cells_oris_long (:, ori_pref_long(i),i);
    end
end    

archive_long_t {1,e} = cells_long;
archive_long_t {2,e} = Is_ori_select_long;
archive_long_t {3,e} = ori_pref_long;
archive_long_t {4,e} = dff_long_prefori;
archive_long_t {5,e} = Is_resp_long;


%add S1 before proceding

%%

[archive_long] = select_long_pairs(archive_long_t);

clear archive_long_t 
%% Study neurons that only respond in one session - S0 or S1

%S0
length(unique(archive_long{1,1}))
[a,b]=unique(archive_long{1,1});
duplicate_indices = setdiff( 1:numel(archive_long{1,1}), b );

num_cells = size (archive_activity{1,1},2);
cell_S0 = [1:num_cells];
cell_S0(1,archive_long{1,1})= nan;
cell_S0 (isnan (cell_S0)) = [];

%cell_S0_select = archive_select{1,1}(1,cell_S0);
cell_S0_prefori = archive_select{3,1}(1,cell_S0);
cell_S0_dff = archive_dff{1,1}(:,cell_S0);
cell_S0_resp = archive_resp{1,1}(1,cell_S0);

cell_S0_dff = reshape (cell_S0_dff, 3360,6,[]);

for i = 1:size (cell_S0,2)
    if isnan(cell_S0_prefori (i))
       cell_S0_dff_prefori(:,i) = nan (3360,1);
    else
        cell_S0_dff_prefori(:,i) = cell_S0_dff (:, cell_S0_prefori (i),i);
    end
end 

%select responsive neurons
col = find (cell_S0_resp ==0);
cell_S0_prefori (col)=[];
cell_S0_dff_prefori (:, col) = [];
cell_S0_resp (col) = [];

clear col

%S1
length(unique(archive_long{1,2}))
[a,b]=unique(archive_long{1,2});
duplicate_indices = setdiff( 1:numel(archive_long{1,2}), b );

num_cells = size (archive_activity{1,2},2);
cell_S1 = [1:num_cells];
cell_S1(1,archive_long{1,2})= nan;
cell_S1 (isnan (cell_S1)) = [];

%cell_S1_select = archive_select{1,2}(1,cell_S1);
cell_S1_prefori = archive_select{3,2}(1,cell_S1);
cell_S1_dff = archive_dff{1,2}(:,cell_S1);
cell_S1_resp = archive_resp{1,2}(1,cell_S1);

cell_S1_dff = reshape (cell_S1_dff, 3360,6,[]);

for i = 1:size (cell_S1,2)
    if isnan(cell_S1_prefori (i))
       cell_S1_dff_prefori(:,i) = nan (3360,1);
    else
        cell_S1_dff_prefori(:,i) = cell_S1_dff (:, cell_S1_prefori (i),i);
    end
end 

%select responsive neurons
col = find (cell_S1_resp ==0);
cell_S1_prefori (col)=[];
cell_S1_dff_prefori (:, col) = [];
cell_S1_resp (col) = [];

clear col

%save responsive cells
archive_nolong = {};
archive_nolong {1,1} = cell_S0;
archive_nolong {2,1} = cell_S0_resp;
archive_nolong {3,1} = cell_S0_dff_prefori;
archive_nolong {4,1} = cell_S0_prefori;

archive_nolong {1,2} = cell_S1;
archive_nolong {2,2} = cell_S1_resp;
archive_nolong {3,2} = cell_S1_dff_prefori;
archive_nolong {4,2} = cell_S1_prefori;

clear cell_S0 cell_S0_select cell_S0_prefori cell_S0_dff_prefori cell_S0_resp cell_S0_dff
clear cell_S1 cell_S1_select cell_S1_prefori cell_S1_dff_prefori cell_S1_resp cell_S1_dff

% In case there are duplicated values
%[a,b]=unique(archive_long{1,2});
%duplicate_indices = setdiff( 1:numel(archive_long{1,2}), b );


 %% Control for frequency: 
 %New phase reversal frequency 1 Hz. Signal should disappear at 1Hz and appear at 0.5 Hz
 %New phase reversal frequency 4 Hz. Signal should disappear at 1Hz and appear at 2 Hz
 
% f_int = 50; % position of ~0.5Hz in frequency vector
 f_int = 188; % position of ~2Hz in frequency vector
 f1Hz = 98; % position of ~1Hz in frequency vector
%  
 [mean_PYY, mean_PYY_resp, dff_resp_G,idx_resp,mean_PYY_1Hz, mean_PYY_1Hz_resp] = response_stimfreq_v2 (PYY_pref,PYY_G, f_int, f1Hz, dff_select_prefori);
% 
 percent_cells_resp = size (dff_resp_G, 2)/sum (Is_active);

%% Power analysis over gray screen and gratings in 1s bins of ROIs that respond at 1Hz
clear Pyy
clear PYY_t
clear t
clear f

dff_resp_prefori_2 = reshape (dff_resp_G,30,112,[]);

for j = 1:size (dff_resp_prefori_2,3)
    for i = 1:size (dff_resp_prefori_2,2)
        Y = fft(dff_resp_prefori_2(:,i,j),30);
        Pyy = Y.*conj(Y)/30;
        f = 30.04/30*(0:16);
        

%         figure;plot(f(100:205),Pyy(100:205))
%         title('Power spectral density')
%         xlabel('Frequency (Hz)')
%         ylim([0 0.025])

        PYY_t(:,i,j)= Pyy;
    end
end

PYY_t_mean(:,:)= mean (PYY_t,3);

t= [1:112];
figure; h = heatmap (t,f(1:5),PYY_t_mean(1:5,:),'Colormap',jet);

%% Power analysis over gray screen and gratings in 10s bins of ROIs that respond at 1Hz
clear Pyy
clear PYY_t
clear PYY_t_mean
clear t
clear f

dff_resp_prefori_2 = dff_resp_G(1:3300, :);
dff_resp_prefori_2 = reshape (dff_resp_prefori_2,300,11,[]);

for j = 1:size (dff_resp_prefori_2,3)
    for i = 1:size (dff_resp_prefori_2,2)
        Y = fft(dff_resp_prefori_2(:,i,j),300);
        Pyy = Y.*conj(Y)/300;
        f = 30.04/300*(0:151);
        

%         figure;plot(f(100:205),Pyy(100:205))
%         title('Power spectral density')
%         xlabel('Frequency (Hz)')
%         ylim([0 0.025])

        PYY_t(:,i,j)= Pyy;
    end
end

PYY_t_mean(:,:)= mean (PYY_t,3);

figure; plot(f(6:152),PYY_t_mean([6:152],2))
axis ([0.5 15 -0.5 3])
ylabel('Power spectral density')
xlabel('Frequency (Hz)')

t= [10:10:110];
% figure; h = heatmap (t,f(11),PYY_t_mean(11,:),'Colormap',jet);

PYY_t_mean_1Hz = PYY_t_mean(11,:);

figure; h = heatmap (t,f(11),PYY_t_mean_1Hz,'Colormap',flipud(gray));
caxis (h, [0 2.5])
ylabel('Power spectral density at 1Hz')
xlabel('Time (s)')


% figure; plot(t,PYY_t_mean_1Hz)
% % axis ([0.5 15 -0.5 10])
% ylabel('Power spectral density')
% xlabel('Frequency (Hz)')
% hold

%% Select responsive cells to a particualr orientation

clear S0_cell_num S1_cell_num S2_cell_num resp_prefori resp_oris

for e = 1:size(archive_select,2)
resp_prefori = archive_select{3,e}(find(archive_resp{1,e}));
resp_oris{e,1} = archive_resp{2,e}(:,find(resp_prefori==1)); %select cells that respond to 30degrees
resp_oris{e,2} = archive_resp{2,e}(:,find(resp_prefori==2)); %select cells that respond to 60degrees
resp_oris{e,3} = archive_resp{2,e}(:,find(resp_prefori==3)); %select cells that respond to 90degrees
resp_oris{e,4} = archive_resp{2,e}(:,find(resp_prefori==4)); %select cells that respond to 1200degrees
resp_oris{e,5} = archive_resp{2,e}(:,find(resp_prefori==5)); %select cells that respond to 150degrees
resp_oris{e,6} = archive_resp{2,e}(:,find(resp_prefori==6)); %select cells that respond to 180degrees
end

for l = 1:size (resp_oris,2)
    S0_cell_num (1,l) = size (resp_oris{1,l},2);
    S1_cell_num (1,l) = size (resp_oris{2,l},2);
    S2_cell_num (1,l) = size (resp_oris{3,l},2);
end

%% Study ROIs that respond at 1Hz during their prefered orientation - total trace
clear ori_pref_resp
clear summary_resp
clear summary_frames_resp

%Orientation preference
% ori_pref_resp = oris (idx_resp);
% 
% for i = 1:1:6
%    ori_pref_resp_dist (i)=sum (ori_pref_resp == i); 
% end  

for i = 1 : size(archive_resp,2)

    dff_resp_prefori = archive_resp {2,i}; %preferred orientation
    % dff_resp_prefori = archive_resp {6,i}; %orthogonal orientation

    % Mean and SD of responsive ROIs during their prefered orientation
    summary_resp (i) = mean(mean(dff_resp_prefori));
    % summary_resp (2,:) = std (dff_resp_prefori);
    % summary_resp (3,:) = max (dff_resp_prefori);
    %resp_threshold = resp_summary (1,:) + 2* resp_summary (2,:);

    %Study frames above threshold of responsive ROIs during their prefered
    %orientation
    [summary_frames, ~] =  select_Ca_transients (dff_resp_prefori);
    summary_frames_resp (:,i) = summary_frames(:,1);
    
     clear dff_resp_prefori summary_frames
end    

%% Study ROIs response to gray screen - total trace
% Analysis of Ca traces from all active neurons during gray screen period
clear ori_pref_resp
clear summary_resp_G
clear summary_frames_resp_G

a=1;

for i = 1 : size(archive_dff,2)
    % i=1;

    dff_resp_G = rmmissing(archive_dff {2,a},2); 

    % Mean and SD of responsive ROIs during their prefered orientation
    summary_resp_G (i) = mean(mean(dff_resp_G));
    % summary_resp (2,:) = std (dff_resp_prefori);
    % summary_resp (3,:) = max (dff_resp_prefori);
    %resp_threshold = resp_summary (1,:) + 2* resp_summary (2,:);

    %Study frames above threshold of responsive ROIs during their prefered
    %orientation
    [summary_frames_G, ~] =  select_Ca_transients (dff_resp_G);
    summary_frames_resp_G (:,i) = summary_frames_G(:,1);

    a = a+1;
    
    clear dff_resp_G summary_frames_G
end    


%% Put together all mice
%  PYY_poll (1,:) = PYY_t_mean(1,:);
% t= [1:112]; m = [1:7];
% figure; h = heatmap (t,1,PYY_poll,'Colormap',jet);
% caxis (h, [0 100])
% ylabel('Day 1')
% xlabel('Time (s)')

% % Nomalized signal with baseline
% PYY_poll_2=PYY_poll';
% PYY_poll_2 = PYY_poll_2(1:110,:);
% PYY_poll_2 = reshape (PYY_poll_2, 10, []);
% PYY_poll_2 = mean (PYY_poll_2);
% PYY_poll_2 = reshape (PYY_poll_2, 11, []);
% 
% baseline = PYY_poll_2 (1,:);
% B = repmat(baseline,11,1)
% PYY_poll_norm = PYY_poll_2./B;

%% Put together responsive cells from all mice_new

dff_resp_S0_T = [];
dff_resp_S1_T = [];
dff_resp_S2_T = [];
% 
% %only for individual orientations
% dff_resp_S3_T = [];
% dff_resp_S4_T = [];
% dff_resp_S5_T = [];

s = 1; %change depending on session

% dff_resp_S0 = archive_resp{2,1}; %preferred orientation
% dff_resp_S0 = archive_resp{6,1}; %orthogonal orientation
dff_resp_S0 = resp_oris{s,1}; %particular orientation
dff_resp_S0_T = horzcat (dff_resp_S0_T, dff_resp_S0);

% dff_resp_S1 = archive_resp{2,2}; %preferred orientation
% dff_resp_S1 = archive_resp{6,2}; %orthogonal orientation
dff_resp_S1 = resp_oris{s,2}; %particular orientation
dff_resp_S1_T = horzcat (dff_resp_S1_T, dff_resp_S1);

% dff_resp_S2 = archive_resp{2,3}; %preferred orientation
% dff_resp_S2 = archive_resp{6,3}; %orthogonal orientation
dff_resp_S2 = resp_oris{s,3}; %particular orientation
dff_resp_S2_T = horzcat (dff_resp_S2_T, dff_resp_S2);

% only for  individual orientations
dff_resp_S3 = resp_oris{s,4}; %particular orientation
dff_resp_S3_T = horzcat (dff_resp_S3_T, dff_resp_S3);

dff_resp_S4 = resp_oris{s,5}; %particular orientation
dff_resp_S4_T = horzcat (dff_resp_S4_T, dff_resp_S4);

dff_resp_S5 = resp_oris{s,6}; %particular orientation
dff_resp_S5_T = horzcat (dff_resp_S5_T, dff_resp_S5);


clear dff_resp_S0 dff_resp_S1 dff_resp_S2 
clear dff_resp_S3 dff_resp_S4 dff_resp_S5

%% Study ROIs that respond at 1Hz during their prefered orientation - timeline (study habituation)
% j=1;
% for i = 1:30:3360
%    
% %     [summary_frames_resp_t, ~] =  select_Ca_transients_v2 (dff_resp_prefori (i:i+29,:));
%     [summary_frames_resp_t, ~] =  select_Ca_transients_v2 (dff_resp_S0_T (i:i+29,:)); %check session
%     summary_frames_number (j,:) = summary_frames_resp_t(1,:);
%     summary_frames_dff(:,j) = summary_frames_resp_t(2,:);
%      j = j+1;
% 
% end
% 
% [max_num,max_idx] = max(summary_frames_dff');
% [~, idx] = sort (max_idx);
% summary_frames_dff_2 = summary_frames_dff(idx, :);
% 
% t= [1:112]; 
% figure; h = heatmap (t,1:size(summary_frames_dff_2,1),summary_frames_dff_2,'Colormap',hot,'GridVisible','off');
% caxis (h, [0 3])
% % ylabel('Day 1')
% % xlabel('Time (s)')
% 
% summary_frames_dff_mean = mean (summary_frames_dff);
% summary_frames_dff_std = std (summary_frames_dff);
% summary_frames_dff_SEM = summary_frames_dff_std/sqrt(size(summary_frames_dff,1));
% figure; plot (summary_frames_dff_mean);
% 
% y = summary_frames_dff_mean; % your mean vector;
% x = 1:numel(y);
% % std_dev = summary_frames_dff_std;
% % curve1 = y + std_dev;
% % curve2 = y - std_dev;
% SEM = summary_frames_dff_SEM;
% curve1 = y + SEM;
% curve2 = y - SEM;
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% figure; fill(x2, inBetween, [0.9 0.9 0.9]);
% hold on;
% plot(x, y, 'k', 'LineWidth', 2);
% axis ([0 112 0 1.1])
% ylabel('dF/F')
% xlabel('Time (s)')
% 
% %divisions dff
% % s10 = mean(mean (summary_frames_dff(:, 1:10),2));
% % s20 = mean(mean (summary_frames_dff(:, 11:20),2));
% % s60 = mean(mean (summary_frames_dff(:, 51:60),2));
% % s110 = mean(mean (summary_frames_dff(:, 101:110),2));
% % summary_frames_dff_timeline = [s10 s20 s60 s110];
% 
% s9 = mean(mean (summary_frames_dff(:, 5:9),2));
% s16 = mean(mean (summary_frames_dff(:, 12:16),2));
% s23 = mean(mean (summary_frames_dff(:, 19:23),2));
% s30 = mean(mean (summary_frames_dff(:, 26:30),2));
% summary_frames_dff_timeline = [s9 s16 s23 s30];
% 
% clear summary_frames_resp_t summary_frames_number summary_frames_dff
%%
%% Study ROIs that respond at 1Hz during their prefered orientation - timeline (study habituation)(all frames)

% dff_resp_prefori_s = squeeze(mean(reshape (dff_resp_prefori, 30, 112, [])))';
dff_resp_prefori_s = squeeze(mean(reshape (dff_resp_S0_T, 30, 112, [])))'; %check session

[max_num,max_idx] = max(dff_resp_prefori_s');
[~, idx] = sort (max_idx);
dff_resp_prefori_s_2 = dff_resp_prefori_s(idx,:);

t= [1:112]; 
 figure; h = heatmap (t,1:size(dff_resp_prefori_s_2,1),dff_resp_prefori_s_2,'Colormap',hot,'GridVisible','off');
% figure; h = heatmap (t,1:size(dff_resp_prefori_s_2,1),dff_resp_prefori_s_2,'GridVisible','off');
caxis (h, [0 3])

% ylabel('Day 1')
% xlabel('Time (s)')

dff_mean = mean (dff_resp_prefori_s);
dff_std = std (dff_resp_prefori_s);
dff_SEM = dff_std/sqrt(size(dff_resp_prefori_s,1));
figure; plot (dff_mean);

y = dff_mean; % your mean vector;
x = 1:numel(y);
% std_dev = summary_frames_dff_std;
% curve1 = y + std_dev;
% curve2 = y - std_dev;
SEM = dff_SEM;
curve1 = y + SEM;
curve2 = y - SEM;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
figure; fill(x2, inBetween, [0.9 0.9 0.9]);
hold on;
plot(x, y, 'k', 'LineWidth', 2);
axis ([0 112 0 2.5])
ylabel('dF/F')
xlabel('Time (s)')

%divisions dff
% s10 = mean(mean (summary_frames_dff(:, 1:10),2));
% s20 = mean(mean (summary_frames_dff(:, 11:20),2));
% s60 = mean(mean (summary_frames_dff(:, 51:60),2));
% s110 = mean(mean (summary_frames_dff(:, 101:110),2));
% summary_frames_dff_timeline = [s10 s20 s60 s110];

s9 = mean(mean (dff_resp_prefori_s(:, 5:7),2));
s16 = mean(mean (dff_resp_prefori_s(:, 12:14),2));
s23 = mean(mean (dff_resp_prefori_s(:, 19:21),2));
s30 = mean(mean (dff_resp_prefori_s(:, 26:28),2));
s110 = mean(mean (dff_resp_prefori_s(:, 109:111),2));
summary_frames_dff_timeline = [s9 s16 s23 s30 s110];

% % frames_dff_timeline for all cells
% s6_all = mean (dff_resp_prefori_s(:,5:7),2); % baseline
% s13_all = mean (dff_resp_prefori_s(:, 12:14),2); % peak at onset
% s110_all = mean (dff_resp_prefori_s(:, 109:111),2); %peak at offset

% frames_dff_timeline for all cells
s13_all = max (dff_resp_prefori_s(:, 10:14),[],2); % peak at onset
s110_all = max (dff_resp_prefori_s(:, 107:111),[],2); %peak at offset

% %calculate Full Width at Half Maximum (FWHM) for onset peak
% clear fmax_1 fmax_idx_1 halfMax_1 fwhm_onset fmax_2 fmax_idx_2 halfMax_2 fwhm_offset
% 
% [fmax_1, fmax_idx_1] = max(mean (dff_resp_prefori_s(:, 12:14),1));
% fmax_idx_1 = fmax_idx_1+11;
% halfMax_1 = fmax_1/2; %Find the half max value.
% index1 = find(mean (dff_resp_prefori_s(:, 1:fmax_idx_1),1) >= halfMax_1, 1, 'first'); % Find where the data first rises above half the max.
% index2 = find(mean (dff_resp_prefori_s(:, fmax_idx_1:112),1) <= halfMax_1, 1, 'first') + fmax_idx_1-1; % Find where the data first drops below half the max.
% fwhm_onset = index2-index1; % FWHM in indexes
% clear index1 index2

% FWHM onset for individual cells
clear fmax_all fmax_idx_all_1 halfMax_all fwhm_onset_all
% [fmax_all, fmax_idx_all_1] = max (dff_resp_prefori_s(:, 12:14),[],2);
% fmax_idx_all_1 = fmax_idx_all_1 + 11;
[fmax_all, fmax_idx_all_1] = max (dff_resp_prefori_s(:, 10:14),[],2);
fmax_idx_all_1 = fmax_idx_all_1 + 9;
halfMax_all = fmax_all/2;

for i = 1:size (dff_resp_prefori_s,1)
    if fmax_all(i)<0
        fwhm_onset_all(i,1) = "NaN";
    else 
    index1_all = find(dff_resp_prefori_s(i, 1:fmax_idx_all_1(i)) >= halfMax_all(i), 1, 'first');
    index2_all = find(dff_resp_prefori_s(i, fmax_idx_all_1(i):112) <= halfMax_all(i), 1, 'first')+fmax_idx_all_1(i)-1;
    fwhm_onset_all(i,1)=index2_all - index1_all;
    end
    clear index1_all index2_all
end


% %calculate Full Width at Half Maximum (FWHM) for offset peak
% [fmax_2, fmax_idx_2] = max(mean (dff_resp_prefori_s(:, 109:111),1));
% fmax_idx_2 = fmax_idx_2+108;
% halfMax_2 = fmax_2/2; %Find the half max value.
% index1 = find(mean (dff_resp_prefori_s(:, 1:fmax_idx_2),1) <= halfMax_2, 1, 'last'); % Find where the data last rises above half the max.
% index2 = find(mean (dff_resp_prefori_s(:, fmax_idx_2:112),1) >= halfMax_2, 1, 'last') + fmax_idx_2-1; % Find where the data first drops below half the max.
% fwhm_offset = index2-index1; % FWHM in indexes
% clear index1 index2

% FWHM offset for individual cells
clear fmax_all fmax_idx_all_2 halfMax_all fwhm_offset_all
% [fmax_all, fmax_idx_all_2] = max (dff_resp_prefori_s(:, 109:111),[],2);
% fmax_idx_all_2 = fmax_idx_all_2 + 108;
[fmax_all, fmax_idx_all_2] = max (dff_resp_prefori_s(:, 107:111),[],2);
fmax_idx_all_2 = fmax_idx_all_2 + 106;
halfMax_all = fmax_all/2;

for i = 1:size (dff_resp_prefori_s,1)
    if fmax_all(i)<0
        fwhm_offset_all(i,1) = "NaN";
    else    
    index1_all = find(dff_resp_prefori_s(i, 1:fmax_idx_all_2(i)) <= halfMax_all(i), 1, 'last');
    index2_all = find(dff_resp_prefori_s(i, fmax_idx_all_2(i):112) >= halfMax_all(i), 1, 'last')+fmax_idx_all_2(i)-1;
    fwhm_offset_all(i,1)=index2_all - index1_all;
    end
    clear index1_all index2_all
end

%% synchrony
[r2,lags2] = xcorr(dff_resp_S1_T);

max_r2 = max(r2);
norm_r2 = r2./max_r2;

syn2= norm_r2(3360,:);
figure()
histogram (syn2, 'Normalization','probability');

syn_S1 = syn2;

syn2 = reshape(syn2,size(dff_resp_S1_T,2),[]);

figure()
h=heatmap (syn2);


clear dff_resp_prefori_2 dff_resp_prefori_s
clear max_r max_r2

 [h,p] = kstest2(syn_S0,syn_S1)
%% Study ROIS that respond at 1Hz during entire stimulation
clear dff_resp_oris
clear OSI
clear mean_resp_OSI

dff_resp_oris = dff_select_oris (:, idx_resp);

%study orientation selectivity index

[OSI] = tuning_index (dff_resp_oris, ori_pref_resp);
mean_resp_OSI = mean (OSI);
