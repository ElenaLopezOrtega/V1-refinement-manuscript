function [session_mat_file, session_tif_file]= GenerateAllFigurePanels_ScanimageWorkflow_Paninski4_SRP6ori(session_files, activity_channel)
%
% - Generate figure panels for given grating stim session.
% (Added resonant scanner support 2017/10)
% 
% Takes a series of Scanimage/focusstack TIF files and 
% 1) makes a focusstack and extracts metadata
% 2) aligns and cropps image 
% 3) performs session-wise normalization of total fluorescence signal
% 4) uses metadata from focusstack to analyze visual responses to gratings
% 5) provides detailed analysis of visual responsiveness & orientation
% selectivity etc
% 6) saves aligned trial-average tif, aligned cropped tif, mat-file, 
% figures(PPTX) etc.
% 
% Ingie Hong, Johns Hopkins Medical Institute, 2017

%% Error catch
%try
    
%% Load data and experimental specifications
if nargin < 1 || isempty(session_files)
    session_files=select_files();
end
[fl, fn, ~]=fileparts(session_files{1});
session_mat_file=[fl filesep fn '_' num2str(length(session_files)) '_trials'];
session_tif_file=[fl filesep fn '_' num2str(length(session_files)) '_trials'];
disp(['Starting analysis for session: ' session_tif_file])

warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning')
% - Create a stack
%session_files=select_files();
tic
fs = FocusStack( session_files );
toc
tic
%w = warning('query','last');
%warning('off',w.identifier)

if ~exist('activity_channel')
    activity_channel=1;
end

%% Ingie's parameters for v1 analysis
disp('Loading image metadata...')
nTrial = numel(fs.cstrFilenames); % One trial per file
nFramesPerTrial = size(fs,3)/nTrial;

% Check framerate 
if abs(fs.tFrameDuration-0.6349)<0.01
    nImPerBlock=12; %14
    BasePeriod=(3:6); % 3:7
    StimPeriod=(7:12); % 8:14
    tBaseTimeShift=0;
    tsub = 1;                                        % degree of downsampling (for 30Hz imaging rate you can try also larger, e.g. 8-10)
elseif abs(fs.tFrameDuration-0.4096)<0.01
    nImPerBlock=18; %14
    BasePeriod=(4:9); % 3:7
    StimPeriod=(10:18); % 8:14
    tBaseTimeShift=0;
    tsub = 1;                                        % degree of downsampling (for 30Hz imaging rate you can try also larger, e.g. 8-10)
elseif abs(fs.tFrameDuration-0.0333)<0.01
    disp('Detected 30hz resonant scanner image...');
    nImPerBlock=240; %14
    BasePeriod=(61:120); % 3:7
    StimPeriod=(121:240); % 8:14
    tBaseTimeShift=0.080;
    tsub = 10;                                        % degree of downsampling (for 30Hz imaging rate you can try also larger, e.g. 8-10)
else
    error('Find and set framerate-related parameters')
end

%% - Align stack -------------------- From run_pipeline of CMNF
disp('Loading image data and aligning...')
% complete pipeline for calcium imaging data pre-processing
gcp;        % start a parallel engine
foldername = fileparts(session_files{1});
        % folder where all the files are located. Currently supported .tif,
        % .hdf5, .raw, .avi, and .mat files
FOV = size(fs(:,:,1,1));
numFiles = length(session_files);
d3 = 1;

%% motion correct (and save registered h5 files as 2d matrices (to be used in the end)..)
% register files one by one. use template obtained from file n to
% initialize template of file n + 1; 

motion_correct = true;      % perform motion correction
non_rigid = false;           % flag for non-rigid motion correction

template = [];
options_rigid = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'d3', d3 ,'bin_width',100,'max_shift',32,...
    'output_type','h5','h5_filename',fullfile([session_tif_file,'_rig.h5']));

if isempty(dir(fullfile([session_tif_file,'_rig.h5'])))
    %disp('Loading image data..')
    %Y_all=fs(:,:,:,activity_channel);
    disp('Starting motion correction...')
    warning('off','MATLAB:imagesci:hdf5dataset:datatypeOutOfRange')
    
    [M,shifts,template] = normcorre_batch( permute( fs(:,:,:,activity_channel), [2 1 3 4 5]) ,options_rigid,template); 
    save(fullfile([session_tif_file,'_shifts_rig.mat']),'shifts','-v7.3');           % save shifts of each file at the respective subfolder
    
    figure;imagesc(template)
    axis image
    title('Template for motion correction')
else
    disp('Motion correction file detected. Using previous registration results...')
    load(fullfile([session_tif_file,'_shifts_rig.mat']))
end


%% Old alignment code
% tic
% fs.Align(activity_channel, false, 10, [], 1, [1/200 1/20], 15);
% toc
% disp('Loading image data...')
% Y=permute(fs(:,:,:,activity_channel), [2 1 3 4]);
% disp('Saving aligned&concatenated image data...')
% %saveastiff(Y,'registered.tif');
%disp('Finished loading/saving image data.')
%% Use alternative alignment (MOCO, NoRMCorre,SIMA etc)


%% Cropping and Multi-trial normalization

% Crop Y to remove blank edges
crop_pixels = ceil(max(abs(squeeze([shifts.shifts_up])))*1.2)
FOV = FOV-2*crop_pixels; % Crop all edges 

%max(abs(fs.mfFrameShifts))        % Check maximum alignment correction
% disp( [ 'Maxiumum alignment displacement: ' num2str(max(abs(fs.mfFrameShifts))) ] )
% crop_pixels=  ceil ( max(abs(fs.mfFrameShifts(:)))*1.2 );

% Y=Y(21:end-20,21:end-20,:,:); % Crop 20 pixels
%Y=Y(11:end-40,41:end-10,:);
% Y=Y(1+crop_pixels:end-crop_pixels,crop_pixels+1:end-crop_pixels,:);

% Normalize session-to-session intensity values
% [Y, Y_avg] = multi_trial_normalization4(Y, fs);
%multi_trial_normalization3;

%Y=Y_normalized;
%clear Y_normalized
%saveastiff(Y,'registered.tif');
% saveastiff(uint16(Y),[session_tif_file '.tif'] );
% save_tif(Y_avg, [session_tif_file '_trial_avg.tif'])

% Inpaint Y
%inpaintn_divide;

%% Save interim results (mat file)
save(session_mat_file, '-v7.3')

%% downsample h5 files and save into a single memory mapped matlab file
disp('Downsampling and saving aligned image data...')

h5_files = subdir(fullfile([session_tif_file '*_rig.h5'])); 
numFiles = length(h5_files);    
fr = 30;                                         % frame rate
%tsub = 10;                                        % degree of downsampling (for 30Hz imaging rate you can try also larger, e.g. 8-10)
ds_filename = [session_tif_file,'_ds_data.mat'];

if true %isempty(dir(fullfile(ds_filename)))
    data_type = class(read_file(h5_files(1).name,1,1));
    data = matfile(ds_filename,'Writable',true);
    data.Y  = zeros([FOV,0],data_type);
    data.Yr = zeros([prod(FOV),0],data_type);
    data.sizY = [FOV,0];
    Yall=zeros([FOV,size(fs,3)], data_type);
    F_dark = Inf;                                    % dark fluorescence (min of all data)
    batch_size = 2000;                               % read chunks of that size
    batch_size = round(batch_size/tsub)*tsub;        % make sure batch_size is divisble by tsub
    Ts = zeros(numFiles,1);                          % store length of each file
    cnt = 0;                                         % number of frames processed so far
    tt1 = tic;
    for i = 1:numFiles
        name = h5_files(i).name;
        info = h5info(name);
        dims = info.Datasets.Dataspace.Size;
        ndimsY = length(dims);                       % number of dimensions (data array might be already reshaped)
        Ts(i) = dims(end);
        Ysub = zeros(FOV(1),FOV(2),floor(Ts(i)/tsub),data_type);
        data.Y(FOV(1),FOV(2),sum(floor(Ts/tsub))) = zeros(1,data_type);
        data.Yr(prod(FOV),sum(floor(Ts/tsub))) = zeros(1,data_type);
        cnt_sub = 0;
        for t = 1:batch_size:Ts(i)
            Y = bigread2(name,t,min(batch_size,Ts(i)-t+1));    
            F_dark = min(nanmin(Y(:)),F_dark);
            ln = size(Y,ndimsY);
            Y = Y(1+crop_pixels(1):end-crop_pixels(1),crop_pixels(2)+1:end-crop_pixels(2),:); % Crop empty edges
            Y = reshape(Y,[FOV,ln]);
            Yall(:,:,t:t+size(Y,3)-1)=Y;   % To retain non-downsampled data
            Y = cast(downsample_data(Y,'time',tsub),data_type);
            ln = size(Y,3);
            Ysub(:,:,cnt_sub+1:cnt_sub+ln) = Y;
            cnt_sub = cnt_sub + ln;
        end
        data.Y(:,:,cnt+1:cnt+cnt_sub) = Ysub;
        data.Yr(:,cnt+1:cnt+cnt_sub) = reshape(Ysub,[],cnt_sub);
        toc(tt1);
        cnt = cnt + cnt_sub;
        data.sizY(1,3) = cnt;
    end
    data.F_dark =  F_dark;
else
    data = matfile(ds_filename);
    i=1;
    name = h5_files(i).name;
    info = h5info(name);
    dims = info.Datasets.Dataspace.Size;
    ndimsY = length(dims);                       % number of dimensions (data array might be already reshaped)
    Ts(i) = dims(end);
end
Y=data.Y;                % Potentially large memory allocation - downsampled Yall
clear data
clear Ysub
%% Flip sign for data that is inverted
if false
    Yall=-Yall+max(Y(:));
    Y=-Y;
    Y=Y-min(Y(:));
end

%% Save downsampled tif file.
disp('Saving downsampled tif file...')
    
figure;imagesc(mean(Y(:,:,:),3))
axis image
title('Mean raw image')

if isempty(dir(fullfile([session_tif_file '.tif'])))

    Yuint16=uint16(   ( (2^15)/(max(Y(:))-min(Y(:))))  *   (Y - min(Y(:)))  );   % Normalize to uint16 range
    saveastiff(Yuint16, [session_tif_file '.tif'])
    %clear Y                 % Keep Y for CNMF
    clear Yuint16
else
    disp('Downsampled tif file already exists. Skipping..')
end

%% Save registered non-downsampled Yall bin file
disp('Saving registered non-downsampled Yall bin file...')
if isempty(dir(fullfile([session_mat_file '_Yall.bin' ])))
    fileID = fopen([session_mat_file '_Yall.bin' ] ,'w');
    fwrite(fileID,Yall,'int16');
    fclose(fileID);
else
    disp('Yall bin file already exists. Skipping..')
end

%% Apply registration to other channel 
other_channel=setdiff([1 2], activity_channel); 
if isempty(dir(fullfile([session_tif_file '_GR.tif'])))   
    im2=permute(fs(:,:,1:1000,other_channel), [2 1 3 4 5]);
    options_apply = NoRMCorreSetParms('d1',size(im2,1),'d2',size(im2,2),'d3', size(im2,3) ,'bin_width',100,'max_shift',32, 'output_type', 'mat');
    im2_reg = apply_shifts( im2 ,shifts(1:1000),options_apply);
    im2_reg = mean(im2_reg(1+crop_pixels(1):end-crop_pixels(1),crop_pixels(2)+1:end-crop_pixels(2),:),3);
    im2 = mean(im2(1+crop_pixels(1):end-crop_pixels(1),crop_pixels(2)+1:end-crop_pixels(2),:),3);

    figure;imagesc(im2)
    title('Un-aligned image')
    axis image
    figure;imagesc(im2_reg)
    title('Other channel aligned image')
    axis image

    im=mean(loadtiff([session_tif_file '.tif'],1,1000),3);
    figure;imagesc(im);title('Activity channel aligned image');axis image

    saveastiff(uint16(cat(3,im2_reg, im)), [session_tif_file '_GR.tif'])
else
    disp('GR dual channel image file already exists. Skipping..')
    im = loadtiff([session_tif_file '_GR.tif']);
    im2_reg = double(im(:,:,1));
    im = double(im(:,:,2));
    figure;imagesc(im);title('Activity channel aligned image');axis image
    figure;imagesc(im2_reg)
    title('Other channel aligned image')
    axis image
end    
%% Paninski analysis

% Ca demixing & deconvol
%warning('off','MATLAB:nargchk:deprecated')
%Paninski_script4; % Uses Paninski CNMF algorithm for spatial segmentation/deconvolution


%% now run CNMF on patches on the downsampled file, set parameters first
disp('Starting CNMF analysis on downsampled aligned image data in memmap...')
sizY = size(Y);                       % size of data matrix
patch_size = [128,128];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [8,8];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 30;                                           % number of components to be found (Ingie: use 10 for sparse samples and 20 for dense samples.)
%K = 20;
tau = 16;                                          % std of gaussian kernel (size of neuron) 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'cluster_pixels',false,...
    'ssub',2,...                                % spatial downsampling when processing 
    'tsub',4,...                                % further temporal downsampling when processing
    'fudge_factor',0.96,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',tau,... 
    'max_size_thr',300,'min_size_thr',80,...   % max/min acceptable size for each component (default: 300,9)
    'space_thresh', 0.4,'time_thresh', 0.4,...  % rval_space and rval_time threshold for each component (default: 0.4,0.4) 
    'max_pr_thr', 0.5,...                       % max_pr threshold for each component (default: 0.9)  
    'spatial_method','regularized',...          % method for updating spatial components
    'df_prctile',50,...                         % take the median of background fluorescence to compute baseline fluorescence 
    'fr',fr/tsub...
    );

%% Run on patches (around 15 minutes)

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(Y,K,patches,tau,p,options);

%% merge found components
%[Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components([],A,b,C,f,P,S,options); % Tried but didn't work 

%% compute correlation image on a small sample of the data (optional - for visualization purposes) 
Cn = correlation_image_max(single(Y),8);

%% classify components
[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(Y,A,C,b,f,YrA,options);

%keep(:) = true; %% TURNING OFF SELECTION OF ROIS!
%% run GUI for modifying component selection (optional, close twice to save values)
 run_GUI = false;
%run_GUI = true;
figure; Coor = plot_contours_ingie2(A,Cn,options,1); 
title('Initial segmented ROIs')
if run_GUI
    
    GUIout = ROI_GUI(Y,A,P,options,Cn,C,b,f);
    options = GUIout{2};
    keep = GUIout{3};    
end

%% view contour plots of selected and rejected components (optional)
throw = ~keep;
figure;
    ax1 = subplot(121); plot_contours_ingie2(A(:,:),Cn,options,1,[],Coor,1,find(keep)); title('Selected components','fontweight','bold','fontsize',14);
    ax2 = subplot(122); plot_contours_ingie2(A(:,:),Cn,options,1,[],Coor,1,find(throw));title('Rejected components','fontweight','bold','fontsize',14);
    linkaxes([ax1,ax2],'xy')

%% manually update components 
prompt = 'Which selected components should be eliminated? ';
x = input(prompt);
if ~isempty(x)
    keep (x)=false;
end
 
prompt = 'Which rejected components should be added? ';
y = input(prompt);
if ~isempty(y)
    keep(y)=true;
end    

%% view contour plots of selected and rejected components (optional)
throw = ~keep;
figure;
    ax1 = subplot(121); plot_contours_ingie2(A(:,:),Cn,options,1,[],Coor,1,find(keep)); title('Manually Selected components','fontweight','bold','fontsize',14);
    ax2 = subplot(122); plot_contours_ingie2(A(:,:),Cn,options,1,[],Coor,1,find(throw));title('Manually Rejected components','fontweight','bold','fontsize',14);
    linkaxes([ax1,ax2],'xy')
%% keep only the active components    
A_keep = A(:,keep);
C_keep = C(keep,:);

%% deconvolve (downsampled) temporal components plot GUI with components (optional)

%tic;
%[C_keep,f_keep,Pk,Sk,YrAk] = update_temporal_components_fast(data,A_keep,b,C_keep,f,P,options);
%toc

%plot_components_GUI(data,A_keep,C_keep,b,f,Cn,options)

%% extract fluorescence and DF/F on native temporal resolution
% C is deconvolved activity, C + YrA is non-deconvolved fluorescence 
% F_df is the DF/F computed on the non-deconvolved fluorescence

P.p = 2;                    % order of dynamics. Set P.p = 0 for no deconvolution at the moment
C_us = cell(numFiles,1);    % cell array for thresholded fluorescence
f_us = cell(numFiles,1);    % cell array for temporal background
P_us = cell(numFiles,1);  
S_us = cell(numFiles,1);
YrA_us = cell(numFiles,1);  % 
b_us = cell(numFiles,1);    % cell array for spatial background
for i = 1:numFiles    
    int = sum(floor(Ts(1:i-1)/tsub))+1:sum(floor(Ts(1:i)/tsub));
    Cin = imresize([C_keep(:,int);f(:,int)],[size(C_keep,1)+size(f,1),Ts(i)]);
    [C_us{i},f_us{i},P_us{i},S_us{i},YrA_us{i}] = update_temporal_components_fast(Yall,A_keep,b,Cin(1:end-1,:),Cin(end,:),P,options);
    b_us{i} = max(mm_fun(f_us{i},Yall) - A_keep*(C_us{i}*f_us{i}'),0)/norm(f_us{i})^2;
end

prctfun = @(data) prctfilt(data,30,1000,300);       % first detrend fluorescence (remove 20th percentile on a rolling 1000 timestep window)
F_us = cellfun(@plus,C_us,YrA_us,'un',0);           % cell array for projected fluorescence
Fd_us = cellfun(prctfun,F_us,'un',0);               % detrended fluorescence

Ab_d = cell(numFiles,1);                            % now extract projected background fluorescence
for i = 1:numFiles
    Ab_d{i} = prctfilt((bsxfun(@times, A_keep, 1./sum(A_keep.^2))'*b_us{i})*f_us{i},30,1000,300,0);
end
    
F0 = cellfun(@plus, cellfun(@(x,y) x-y,F_us,Fd_us,'un',0), Ab_d,'un',0);   % add and get F0 fluorescence for each component
F_df = cellfun(@(x,y) x./y, Fd_us, F0 ,'un',0);                            % DF/F value
%% detrend each segment and then deconvolve


% %% perform deconvolution
% Cd = cellfun(@(x) zeros(size(x)), Fd_us, 'un',0);
% Sp = cellfun(@(x) zeros(size(x)), Fd_us, 'un',0);
% bas = zeros(size(Cd{1},1),numFiles);
% c1 = bas;
% sn = bas;
% gn = cell(size(bas));
% options.p = 2;
% tt1 = tic;
% for i = 1:numFiles
%     c_temp = zeros(size(Cd{i}));
%     s_temp = c_temp;
%     f_temp = Fd_us{i};
%     parfor j = 1:size(Fd_us{i},1)
%         [c_temp(j,:),bas(j,i),c1(j,i),gn{j,i},sn(j,i),s_temp(j,:)] = constrained_foopsi(f_temp(j,:),[],[],[],[],options);
%     end
%     Cd{i} = c_temp;
%     Sp{i} = s_temp;
%     toc(tt1);
% end
figure; Coor = plot_contours_ingie2(A_keep,Cn,options,1);
%plot_components_GUI(double(Yall),A_keep,Cin(1:end-1,:),b_us{1},f_us{1},Cn,options)
title('Final segmented ROIs')
disp('Finished CNMF analysis.')

%% Select neurons


%% Use metadata to de-randomize calcium data
mfRegionTraces_denoised_df_f=F_df{1};
save ([session_tif_file,'_mfRegionTraces_df_f.mat'], 'mfRegionTraces_denoised_df_f');
%mfRegionTraces_denoised_df_f = full(C_df);
%mfRegionTraces_denoised_df_f = C;
mfWholeFOV_df_f= reshape(mean(mean(Yall,1),2), 1, [])/mean(mean(mean(Yall,1),2),3);
save ([session_tif_file,'_mfWholeFOV_df_f.mat'], 'mfWholeFOV_df_f');

% - Make a matrix of trial responses

% vnGratingStimuli = 1:13;
% nCellToUse = 82;
% nCellToUse_denoised = 1;
% fXScale = 3/50;
% fYScale = 1/2.5;
% nTraceWidth = 1;
% fdFFScaleBarLength = 200;
% tTimeScaleBarLength = 2;
% tTimeShift = 0;
% tTickLength = 0.2;

figure;plot(mfRegionTraces_denoised_df_f')
set(gcf, 'Position',  [1, 41, 1920, 964])
title('mfRegionTraces_denoised_df_f','Interpreter', 'none')
legend('Location', 'eastoutside')
axis tight



%% Export results to Powerpoint
%sortfigs % now included in exportFigsToPPTX_ScanimageWorkflow
exportFigsToPPTX_ScanimageWorkflow

clear Yall
clear Y
clear Yr
save(session_mat_file, '-v7.3')

end
