function [fval] = long_term_changes (val, sessions)


fval = val(:, sessions);

fval(isnan(fval)) = 0;
fval_mean = mean (fval, 2);
del = fval_mean ==0;
fval(del, :) = []; %eliminate segments that are not present in any of the sessions considered in this section

fval_mean_exp = mean(fval(:,[3 4]),2);
del2 = fval_mean_exp ==0;
fval(del2, :) = []; %eliminate segments that are not present in either experimetal sessions (D3 D5)

%  effect = fval(:,3);
effect = mean(fval(:,[3 4]),2);
[effect,ix] = sort (effect, 'descend');
fval = fval (ix,:); %organize segments according to the average of day 5 and day 7
fval(fval ==0) = NaN;


t= [0 1 5 7];
segm = (1:size (fval, 1));
map = [0 0 1; 0.01 0.01 0.99; 0.02 0.02 0.98; 0.03 0.03 0.97; 0.04 0.04 0.96; 0.05 0.05 0.95; 0.1 0.1 0.95; 0.15 0.15 0.94; 0.2 0.2 0.93; 0.25 0.25 0.92; 0.3 0.3 0.91; 0.35 0.35 0.90; 0.4 0.4 0.89; 0.45 0.45 0.88; 0.5 0.5 0.87; 0.55 0.55 0.86; 0.6 0.6 0.85; 0.65 0.65 0.84; 0.7 0.7 0.83; 0.75 0.75 0.82; 0.8 0.8 0.81; 0.8 0.8 0.8; 0.81 0.8 0.8; 0.82 0.75 0.75; 0.83 0.7 0.7; 0.84 0.65 0.65; 0.85 0.6 0.6; 0.86 0.55 0.55; 0.87 0.5 0.5; 0.88 0.45 0.45; 0.89 0.4 0.4; 0.9 0.35 0.35; 0.91 0.3 0.3; 0.92 0.25 0.25; 0.93 0.2 0.2; 0.94 0.15 0.15; 0.95 0.1 0.1; 0.95 0.05 0.05; 0.96 0.04 0.04; 0.97 0.03 0.03; 0.98 0.02 0.02; 0.99 0.01 0.01; 1 0 0];
figure; h = heatmap (t,segm,fval,'Colormap',map, 'MissingDataColor','1.00,1.00,1.00');
caxis (h, [0 2])
h.GridVisible = 'off';
