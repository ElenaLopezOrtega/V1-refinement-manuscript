function [frames_summary_mouse, cell_summary, Is_active] =  select_Ca_transients (dff_cells_oris)

%Select frames and local peaks that are above a threshold (2 ZScores).
%Individual traces can be represented

clear frames_summary
clear frames_summary_mouse
clear cell_summary
clear active_cells_500_idx
clear active_cells_650_idx
clear active_cells_800_idx


% dff_norm values:
dff_cells_mean = mean (dff_cells_oris);
dff_cells_stdv = std (dff_cells_oris);
threshold_cells = dff_cells_mean + 3*dff_cells_stdv; %Z-score of 3
% threshold_cells = dff_cells_mean + 1*dff_cells_stdv; %Z-score of 1


active_neurons = 0;
for i = 1:size(dff_cells_oris, 2)
    clear frames
    clear lk
    clear pk 
    
    frames_idx = dff_cells_oris(:,i)> threshold_cells(i); %find frames above threshold
    frames = dff_cells_oris(frames_idx,i);
    if isempty(frames)
        frames_summary (1,i)= NaN;
        frames_summary (2,i)= NaN;
        frames_summary (3,i)= NaN;
        Is_active(i) = 0;
    else
    frames_summary (1,i)= size (frames,1);
    frames_summary (2,i)= mean (frames);
    frames_summary (3,i) = max (frames);
    Is_active(i)=1;
    end

%     [pk, lk] = findpeaks(dff_cells_oris(:,i)); %find local peaks above threshold
%     below_thres = find (pk<=threshold_cells(i));
%     pk(below_thres) =0; lk(below_thres) =0;
%     pk = nonzeros(pk);lk = nonzeros(lk);
%     peaks_summary(1,i) = size (pk,1);
%     peaks_summary(2,i) = mean(pk);
%     peaks_summary (3,i) = max (pk);


% %     Plot:
% %   plot cell traces including mean response, threshold and local peaks
%     figure; plot (1:size (dff_cells_oris,1),dff_cells_oris(:,i))
%     xlim([1 size(dff_cells_oris,1)]);
% %     figure; plot (1:size(dff_cells_oris,1),dff_cells_oris(:,i),lk,pk, 'o')
% %     figure; plot (1:3360,dff(:,i)-baseline(i,1),lk,amp, 'o')
%     line([1,size(dff_cells_oris,1)],[dff_cells_mean(i),dff_cells_mean(i)],'Color','blue')
%     line([1,size(dff_cells_oris,1)],[threshold_cells(i),threshold_cells(i)],'Color','red','LineStyle','--') 


%     figure; plot (1:3360,dff_cells_oris(:,i))
%     xlim([1 3360]);
%     figure; plot (1:20160,dff_cells_oris(:,i),lk,pk, 'o')
%     figure; plot (1:3360,dff(:,i)-baseline(i,1),lk,amp, 'o')
%     line([1,20160],[dff_cells_mean(i),dff_cells_mean(i)],'Color','blue')
%     line([1,20160],[threshold_cells(i),threshold_cells(i)],'Color','red','LineStyle','--') 
end

frames_summary= rmmissing(frames_summary,2); %eliminate cells with NaN values

prct25 = prctile(frames_summary (1,:),25); %max number of frames above threshold for lower 25% of the population of cells
prct50 = prctile(frames_summary (1,:),50); %max number of frames above threshold for lower 50% of the population of cells
prct75 = prctile(frames_summary (1,:),75); %max number of frames above threshold for lower 50% of the population of cells

% Determine different thersholds for activity based on number of frames
% that are above 2ZScores

active_cells_500_idx = find (frames_summary(1,:)>500);
frames_summary_500 = frames_summary (:, active_cells_500_idx);

active_cells_650_idx = find (frames_summary(1,:)>650);
frames_summary_650 = frames_summary (:, active_cells_650_idx);

active_cells_800_idx = find (frames_summary(1,:)>800);
frames_summary_800 = frames_summary (:, active_cells_800_idx);

frames_summary_mouse(:,1) = mean (frames_summary,2);
frames_summary_mouse(:,2) = mean (frames_summary_500,2);
frames_summary_mouse(:,3) = mean (frames_summary_650,2);
frames_summary_mouse(:,4) = mean (frames_summary_800,2);

% Number of cells
cell_summary (1) = sum(Is_active);
cell_summary (2) = size (active_cells_500_idx, 2);
cell_summary (3) = size (active_cells_650_idx, 2);
cell_summary (4) = size (active_cells_800_idx, 2);