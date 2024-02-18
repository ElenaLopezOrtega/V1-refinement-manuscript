function [fval_prc, s] = increase_decrease_spines (threshold_inc, threshold_dec, val, sessions)

if size (val,2)>4
    fval = val(:, sessions);
else
    fval=val;
end

fval(isnan(fval)) = 0;
fval_mean = mean (fval, 2);
del = fval_mean ==0;
fval(del, :) = []; 

fval_mean_exp = mean(fval(:,[3 4]),2);
del2 = fval_mean_exp ==0;
fval(del2, :) = [];
fval(fval ==0) = NaN;

s = size (fval,1);

fval_increase = find (fval (:, 3)> threshold_inc);
fval_prc (1) = size (fval_increase, 1)/ size (fval,1);

fval_decrease = find (fval (:, 3)<threshold_dec);
fval_prc (2) = size (fval_decrease, 1)/ size (fval,1);
