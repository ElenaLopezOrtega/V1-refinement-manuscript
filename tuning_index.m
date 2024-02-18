function [OSI] = tuning_index (dff_oris, ori_pref)

mean_resp_oris =squeeze(mean((reshape (dff_oris, 3360,6,[]))));

if size (mean_resp_oris, 1)==1
    mean_resp_oris = mean_resp_oris';
end  

for i = 1: size (ori_pref)
    if ori_pref(i) <=3 
        ori_orth(i) = ori_pref(i)+3;
    else
        ori_orth(i) = ori_pref(i)-3;
    end
end    

ori_orth= ori_orth';

for i = 1:size (mean_resp_oris,2)
    pref = mean_resp_oris (ori_pref(i),i);
    orth = mean_resp_oris (ori_orth(i),i);
    OSI (i) = (pref-orth)/orth;
end   
