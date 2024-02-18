function [neuron_total,ix, segm_total] = segments_neuron (val, cell_number, sessions)

fval = val(:, sessions);

fval(isnan(fval)) = 0;
fval_mean = mean (fval, 2);
del = fval_mean ==0;
fval(del, :) = []; 
cell_number (del) = [];

fval_mean_exp = mean(fval(:,[3 4]),2);
del2 = fval_mean_exp ==0;
fval(del2, :) = [];
cell_number (del2) = [];

% Averaged neuronal values
s = size (cell_number);

for i = 1 : cell_number (s(1))
    idx2 = cell_number ==i;
    neuron = mean (fval(idx2,:),1);
    neuron_total (i,:) = neuron;
end    

effect = neuron_total(:,3);
[effect,ix] = sort (effect, 'descend');
neuron_total = neuron_total (ix,:); %organize neurons according to day 5 signal
neuron_total(neuron_total ==0) = NaN;

t= [0 1 5 7];
n = (1:size (neuron_total, 1));
map = [0 0 1; 0.01 0.01 0.99; 0.02 0.02 0.98; 0.03 0.03 0.97; 0.04 0.04 0.96; 0.05 0.05 0.95; 0.1 0.1 0.95; 0.15 0.15 0.94; 0.2 0.2 0.93; 0.25 0.25 0.92; 0.3 0.3 0.91; 0.35 0.35 0.90; 0.4 0.4 0.89; 0.45 0.45 0.88; 0.5 0.5 0.87; 0.55 0.55 0.86; 0.6 0.6 0.85; 0.65 0.65 0.84; 0.7 0.7 0.83; 0.75 0.75 0.82; 0.8 0.8 0.81; 0.8 0.8 0.8; 0.81 0.8 0.8; 0.82 0.75 0.75; 0.83 0.7 0.7; 0.84 0.65 0.65; 0.85 0.6 0.6; 0.86 0.55 0.55; 0.87 0.5 0.5; 0.88 0.45 0.45; 0.89 0.4 0.4; 0.9 0.35 0.35; 0.91 0.3 0.3; 0.92 0.25 0.25; 0.93 0.2 0.2; 0.94 0.15 0.15; 0.95 0.1 0.1; 0.95 0.05 0.05; 0.96 0.04 0.04; 0.97 0.03 0.03; 0.98 0.02 0.02; 0.99 0.01 0.01; 1 0 0];
figure; h = heatmap (t,n,neuron_total,'Colormap',map, 'MissingDataColor','1.00,1.00,1.00');
caxis (h, [0 2])
h.GridVisible = 'off';


% Organize segments within each individual neuron
segm_total =[];
for i = 1:size (ix,1)
    idx3 = cell_number == ix(i);
    segm = fval(idx3,:);
    segm_eff = segm (:,3);
    [segm_effect,ix2]= sort (segm_eff, 'descend');
    segm = segm (ix2,:);
    
    segm_total= vertcat (segm_total, segm);
end   
segm_total(segm_total ==0) = NaN;

t= [0 1 5 7];
segm = (1:size (fval, 1));
map = [0 0 1; 0.01 0.01 0.99; 0.02 0.02 0.98; 0.03 0.03 0.97; 0.04 0.04 0.96; 0.05 0.05 0.95; 0.1 0.1 0.95; 0.15 0.15 0.94; 0.2 0.2 0.93; 0.25 0.25 0.92; 0.3 0.3 0.91; 0.35 0.35 0.90; 0.4 0.4 0.89; 0.45 0.45 0.88; 0.5 0.5 0.87; 0.55 0.55 0.86; 0.6 0.6 0.85; 0.65 0.65 0.84; 0.7 0.7 0.83; 0.75 0.75 0.82; 0.8 0.8 0.81; 0.8 0.8 0.8; 0.81 0.8 0.8; 0.82 0.75 0.75; 0.83 0.7 0.7; 0.84 0.65 0.65; 0.85 0.6 0.6; 0.86 0.55 0.55; 0.87 0.5 0.5; 0.88 0.45 0.45; 0.89 0.4 0.4; 0.9 0.35 0.35; 0.91 0.3 0.3; 0.92 0.25 0.25; 0.93 0.2 0.2; 0.94 0.15 0.15; 0.95 0.1 0.1; 0.95 0.05 0.05; 0.96 0.04 0.04; 0.97 0.03 0.03; 0.98 0.02 0.02; 0.99 0.01 0.01; 1 0 0];
figure; h = heatmap (t,segm,segm_total,'Colormap',map, 'MissingDataColor','1.00,1.00,1.00');
caxis (h, [0 2])
h.GridVisible = 'off';

