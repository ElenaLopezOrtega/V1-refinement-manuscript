function [s14,s14_T] = plot_SEM (dff)

dff_s = squeeze(nanmean(reshape (dff, 30, 112, [])))';

dff_mean = nanmean (dff_s);
dff_std = std (dff_s);
dff_SEM = dff_std/sqrt(size(dff_s,1));

figure; plot (dff_mean);
ylim([0 3])


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

s14 = nanmean (dff_s(:, 12:14),2);
s14_T = nanmean(nanmean (dff_s(:, 12:14),2));