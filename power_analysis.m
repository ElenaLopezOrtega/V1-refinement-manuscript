function [PYY, f] = power_analysis (dff)

clear Pyy
clear PYY
clear t
clear f

for i = 1:size (dff,2)
    Y = fft(dff(:,i),3005);
    Pyy = Y.*conj(Y)/3005;
    f = 30.04/3005*(0:1503);
    

%     figure; plot(f,Pyy(1:1504))
%     title('Power spectral density')
%     xlabel('Frequency (Hz)')
%     
%     figure;plot(f(100:205),Pyy(100:205))
%     title('Power spectral density')
%     xlabel('Frequency (Hz)')
% %     ylim([0 0.025])
    
    PYY(:,i)= Pyy;
end

% Heatmap 1 Hz signal
t= [1:size(dff,2)]; %value in x is # of cells
figure; h = heatmap (t,f(82:122),PYY(82:122,:),'Colormap',jet);
%caxis (h, [0 29])

% % Heatmap 0.5 Hz signal
% t= [1:size(dff,2)];
% figure; h = heatmap (t,f(32:72),PYY(32:72,:),'Colormap',jet);

% % Heatmap 2 Hz signal
% t= [1:size(dff,2)];
% figure; h = heatmap (t,f(182:222),PYY(182:222,:),'Colormap',jet);

PYY_mean = mean (PYY');
PYY_std = std (PYY');
PYY_SEM = PYY_std/sqrt(size(PYY,2));

 figure; plot(f(51:1504),PYY_mean(51:1504))
% figure; plot(f(41:1504),PYY_mean(41:1504))
ylim([0 6])
title('Power spectral density')
xlabel('Frequency (Hz)')

 y = PYY_mean(51:1504); % your mean vector;
 x = f(51:1504);
% y = PYY_mean(41:1504); % your mean vector;
% x = f(41:1504);
% std_dev = summary_frames_dff_std;
% curve1 = y + std_dev;
% curve2 = y - std_dev;
SEM = PYY_SEM (51:1504);
% SEM = PYY_SEM (41:1504);
curve1 = y + SEM;
curve2 = y - SEM;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
figure; fill(x2, inBetween, [0.9 0.9 0.9]);
hold on;
plot(x, y, 'k', 'LineWidth', 2);
ylim([0 6])
title('Power spectral density')
xlabel('Frequency (Hz)')
