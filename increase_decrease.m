function [fval_prc, cell_inc, cell_dec] = increase_decrease (threshold_inc, threshold_dec, val, sessions, cell_number)

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
fval(fval ==0) = NaN;

fval_increase = find (fval (:, 3)> threshold_inc);
fval_prc (1) = size (fval_increase, 1)/ size (fval,1);
cell_number_inc = cell_number(fval_increase);
cell_inc (:,1) = unique(cell_number_inc); % which will give you the unique elements of A in array B
cell_inc (:,2)= histc(cell_number_inc, cell_inc (:,1));

fval_decrease = find (fval (:, 3)<threshold_dec);
fval_prc (2) = size (fval_decrease, 1)/ size (fval,1);
cell_number_dec = cell_number(fval_decrease);
cell_dec (:,1) = unique(cell_number_dec); % which will give you the unique elements of A in array B
cell_dec (:,2)= histc(cell_number_dec, cell_dec (:,1));