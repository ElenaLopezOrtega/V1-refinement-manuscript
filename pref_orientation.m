function [Is_ori_select, ori_pref, ori_orth, dff_cells_oris_select] = pref_orientation (dff_cells_oris, Is_active)

%Study ROIs orientation selectivity
%Individual traces of orientation selective ROIs can be represented

clear dff_cells_oris_select
clear Is_ori_select

dff_cells_oris_dsmp = reshape (dff_cells_oris, 30,672,[]);
dff_cells_oris_dsmp = squeeze(mean (dff_cells_oris_dsmp(:,:,:))); % downsample for anova
dff_cells_oris_dsmp = reshape (dff_cells_oris_dsmp, 112,6,[]);

group = [1 2 3 4 5 6];
for i = 1:size (dff_cells_oris_dsmp,3)
    [p(i)] = anova1(dff_cells_oris_dsmp(:,:,i),group, 'off');
     
    if p(i) < 0.01 && Is_active (i)==1
        Is_ori_select (i) = 1;
        m(i)= max (mean(dff_cells_oris_dsmp(:,:,i)));
        ori_pref (i) = find (m(i) == mean(dff_cells_oris_dsmp(:,:,i)));
    else
        Is_ori_select (i) = 0;
        ori_pref (i) = NaN;
    end
    
end     

ori_orth = NaN (size (ori_pref));
ori_orth (ori_pref<4) = ori_pref(ori_pref<4)+3;
ori_orth (ori_pref>3) = ori_pref(ori_pref>3)-3;

for i = 1:size (dff_cells_oris_dsmp,3)
    if Is_ori_select(i) ==1
        [h(i), p2(i)]=ttest2 (dff_cells_oris_dsmp (:, ori_pref(i),i), dff_cells_oris_dsmp (:, ori_orth(i),i));
    else
        p2(i)=NaN;
    end    
end

Is_ori_select(p2>=0.01)=0;
ori_pref (p2>=0.01) = NaN;
ori_orth (p2>=0.01) = NaN;


% Select orientation selective cells
dff_cells_oris_select = dff_cells_oris (:, Is_ori_select==1);

%plot individual traces of orientation selective cells
% for i = 1:size (dff_cells_oris_select,2)
%     figure; plot (1:20160,dff_cells_oris_select(:,i));
% end    