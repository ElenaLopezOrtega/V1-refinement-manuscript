% cond1 = 'noreach';
% cond2 = 'train';
% cond3 = 'vctx';
% cond4 = 'vcdark';
% cond5 = 'vcnoreach';
% cond6 = 'hcctrl';

%V1 analysis
 cond1 = 'SRP';
 cond2 = 'CNT';

% cond1 = 'SRP3';
% cond2 = 'SRP4';

%SAPs analysis
% cond1 = 'WTCNTSt';
% cond2 = 'WTCNTNS';

%SessionConditions = {'base1', 'base2', 'train1', 'train2', 'train3', 'train4', 'train5', 'train6', 'train7', 'train8'};

%V1 analysis
% SessionConditions = {'base1','base2','postSRP1','preSRP5','postSRP5'};
  SessionConditions = {'base1','base2','preSRP5','h48'};

%SAPs analysis
%SessionConditions = {'base1','base2','hhalf','h1','h2','h3','h24','h48'};


ts = size (SessionConditions,2);
% StatName = 'ntbutsmdodmch1'; %Spine AMPA / Dendrite dsRed mean
StatName = 'pDist';
% StatName = 'userBool1'
channel = 1;
PersCond = 1; %Minimum number of sessions that spine needs to be present in
MinSpines = 10; %Minimum number of spines per dendrite for dendrite to be included

myPool_cond1 = getmyMapListValuesCondDyn_v2(myMapList, mySegmentTable, cond1, SessionConditions, StatName, channel, PersCond, MinSpines);
myPool_cond2 = getmyMapListValuesCondDyn_v2(myMapList, mySegmentTable, cond2, SessionConditions, StatName, channel, PersCond, MinSpines);
% myPool_cond3 = getmyMapListValuesCondDyn_v2(myMapList, mySegmentTable, cond3, SessionConditions, StatName, channel, PersCond, MinSpines);
% myPool_cond4 = getmyMapListValuesCondDyn_v2(myMapList, mySegmentTable, cond4, SessionConditions, StatName, channel, PersCond, MinSpines);
% myPool_cond5 = getmyMapListValuesCondDyn_v2(myMapList, mySegmentTable, cond5, SessionConditions, StatName, channel, PersCond, MinSpines);
% myPool_cond6 = getmyMapListValuesCondDyn_v2(myMapList, mySegmentTable, cond6, SessionConditions, StatName, channel, PersCond, MinSpines);

% plotMean = 'SegmentsSpines';
% plotMyPool(myPool, plotMean)
% spines=removeNanRows(myPool.values.Spines);

%%
myPool=myPool_cond1;
Condition = cond1;
normadd = (myPool.additions./myPool.additions(:,2));
normadd(find(normadd==Inf))=NaN;
normelim = (myPool.eliminations./myPool.eliminations(:,2));
normelim(find(normelim==Inf))=NaN;
normcountall = (myPool.countall./nanmean(myPool.countall(:,1:2),2));
normcountall(find(normcountall==Inf))=NaN;
density = myPool.countall./myPool.dist; %added by ELO on 2/27/21
density(find(density==Inf))=NaN;
normdensity_SRP = (density./nanmean(density(:,1:2),2));
normdensity_SRP(find(normdensity_SRP==Inf))=NaN;

output.(Condition)(1,:) = nanmean(myPool.additions)*100;
output.(Condition)(2,:) = nanmean(myPool.eliminations)*100;
output.(Condition)(3,:) = nanmean(normadd);
output.(Condition)(4,:) = nanmean(normelim);
output.(Condition)(5,:) = nanmean(normcountall);
output.(Condition)(6,:) = nanmean(myPool.survivalrate)*100; %added by ELO on 2/27/21
output.(Condition)(7,:) = nanmean(density); %added by ELO on 2/27/21
output_survivalD1.(Condition) = nanmean(myPool.survivalD1rate);

myPool=myPool_cond2;
Condition = cond2;
normadd = (myPool.additions./myPool.additions(:,2))*100;
normadd(find(normadd==Inf))=NaN;
normelim = (myPool.eliminations./myPool.eliminations(:,2))*100;
normelim(find(normelim==Inf))=NaN;
normcountall = (myPool.countall./nanmean(myPool.countall(:,1:2),2))*100;
normcountall(find(normcountall==Inf))=NaN;
density = myPool.countall./myPool.dist; %added by ELO on 2/27/21
density(find(density==Inf))=NaN;
normdensity_CNT = (density./nanmean(density(:,1:2),2));
normdensity_CNT(find(normdensity_CNT==Inf))=NaN;

output.(Condition)(1,:) = nanmean(myPool.additions)*100;
output.(Condition)(2,:) = nanmean(myPool.eliminations)*100;
output.(Condition)(3,:) = nanmean(normadd);
output.(Condition)(4,:) = nanmean(normelim);
output.(Condition)(5,:) = nanmean(normcountall);
output.(Condition)(6,:) = nanmean(myPool.survivalrate)*100; %added by ELO on 2/27/21
output.(Condition)(7,:) = nanmean(density); %added by ELO on 2/27/21
output_survivalD1.(Condition) = nanmean(myPool.survivalD1rate);


% myPool=myPool_cond3;
% Condition = cond3;
% normadd = (myPool.additions./myPool.additions(:,2))*100;
% normadd(find(normadd==Inf))=NaN;
% normelim = (myPool.eliminations./myPool.eliminations(:,2))*100;
% normelim(find(normelim==Inf))=NaN;
% normcountall = (myPool.countall./nanmean(myPool.countall(:,1:2),2))*100;
% normcountall(find(normcountall==Inf))=NaN;
% 
% output.(Condition)(1,:) = nanmean(myPool.additions)*100;
% output.(Condition)(2,:) = nanmean(myPool.eliminations)*100;
% output.(Condition)(3,:) = nanmean(normadd);
% output.(Condition)(4,:) = nanmean(normelim);
% output.(Condition)(5,:) = nanmean(normcountall);
% output_survivalD1.(Condition) = nanmean(myPool.survivalD1rate);
% 
% myPool=myPool_cond4;
% Condition = cond4;
% normadd = (myPool.additions./myPool.additions(:,2))*100;
% normadd(find(normadd==Inf))=NaN;
% normelim = (myPool.eliminations./myPool.eliminations(:,2))*100;
% normelim(find(normelim==Inf))=NaN;
% normcountall = (myPool.countall./nanmean(myPool.countall(:,1:2),2))*100;
% normcountall(find(normcountall==Inf))=NaN;
% 
% output.(Condition)(1,:) = nanmean(myPool.additions)*100;
% output.(Condition)(2,:) = nanmean(myPool.eliminations)*100;
% output.(Condition)(3,:) = nanmean(normadd);
% output.(Condition)(4,:) = nanmean(normelim);
% output.(Condition)(5,:) = nanmean(normcountall);
% output_survivalD1.(Condition) = nanmean(myPool.survivalD1rate);
% 
% myPool=myPool_cond5;
% Condition = cond5;
% normadd = (myPool.additions./myPool.additions(:,2))*100;
% normadd(find(normadd==Inf))=NaN;
% normelim = (myPool.eliminations./myPool.eliminations(:,2))*100;
% normelim(find(normelim==Inf))=NaN;
% normcountall = (myPool.countall./nanmean(myPool.countall(:,1:2),2))*100;
% normcountall(find(normcountall==Inf))=NaN;
% 
% output.(Condition)(1,:) = nanmean(myPool.additions)*100;
% output.(Condition)(2,:) = nanmean(myPool.eliminations)*100;
% output.(Condition)(3,:) = nanmean(normadd);
% output.(Condition)(4,:) = nanmean(normelim);
% output.(Condition)(5,:) = nanmean(normcountall);
% output_survivalD1.(Condition) = nanmean(myPool.survivalD1rate);
% 
% myPool=myPool_cond6;
% Condition = cond6;
% normadd = (myPool.additions./myPool.additions(:,2))*100;
% normadd(find(normadd==Inf))=NaN;
% normelim = (myPool.eliminations./myPool.eliminations(:,2))*100;
% normelim(find(normelim==Inf))=NaN;
% normcountall = (myPool.countall./nanmean(myPool.countall(:,1:2),2))*100;
% normcountall(find(normcountall==Inf))=NaN;
% 
% output.(Condition)(1,:) = nanmean(myPool.additions)*100;
% output.(Condition)(2,:) = nanmean(myPool.eliminations)*100;
% output.(Condition)(3,:) = nanmean(normadd);
% output.(Condition)(4,:) = nanmean(normelim);
% output.(Condition)(5,:) = nanmean(normcountall);
% output_survivalD1.(Condition) = nanmean(myPool.survivalD1rate);
%%
% addition (number of spines added normalized to number of spines on previous session: %addition)
figure
x= [1:ts];
% plot(x',output.noreach(1,:)','k-')
plot(x',output.SRP(1,:)','k-')
hold on
% plot(x',output.train(1,:)','b-')
plot(x',output.CNT(1,:)','b-')
xlim([0 6]) 
%ylim([0 20]) 
title('addition')

% elimination (number of spines eliminated normalized to number of spines on previous session: %elimination)
figure
x= [1:ts];
% plot(x',output.noreach(2,:)','k-')
plot(x',output.SRP(2,:)','k-')
hold on
plot(x',output.CNT(2,:)','b-')
xlim([0 6]) 
% ylim([0 20]) 
title('elimination')

%normalized addition (% addition normalized to % addition on baseline2)
figure
x= [1:ts];
plot(x',output.SRP(3,:)','k-')
hold on
plot(x',output.CNT(3,:)','b-')
hold on
plot([0:10],repmat(100,11,1),'k--')
xlim([0 10]) 
% ylim([0 150]) 
title('normalized addition')

% normalized elimination (% elimination normalized to % elimination on baseline2)
figure
x= [1:ts];
plot(x',output.SRP(4,:)','k-')
hold on
plot(x',output.CNT(4,:)','b-')
hold on
plot([0:10],repmat(100,11,1),'k--')
xlim([0 10]) 
% ylim([80 150]) 
title('normalized elimination')

%normalized count (number of spines normalized to average baseline numbers)
figure
x= [1:ts];
plot(x',output.SRP(5,:)','k-')
hold on
plot(x',output.CNT(5,:)','b-')
hold on
plot([0:10],repmat(100,11,1),'k--')
xlim([0 10]) 
ylim([80 120]) 

% figure
% x= [1:10];
% plot(x',output.vcdark(5,:)','k-')
% hold on
% plot(x',output.vctx(5,:)','b-')
% % hold on
% % plot(x',output.vcnoreach(5,:)','r-')
% hold on
% plot([0:10],repmat(100,11,1),'k--')
% 
% xlim([0 10]) 
% ylim([80 110]) 
% 
% 
% figure
% x= [1:10];
% plot(x',output_survivalD1.noreach(1,:)','k-')
% hold on
% plot(x',output_survivalD1.train(1,:)','b-')
% hold on
% plot(x',output_survivalD1.hcctrl(1,:)','r-')
% xlim([0 12]) 
% ylim([0 1]) 
% % 
% % 
% % figure
% hold on
% x= [1:10];
% plot(x',output_survivalD1.vcdark(1,:)','k--')
% hold on
% plot(x',output_survivalD1.vctx(1,:)','b--')
% hold on
% plot(x',output_survivalD1.vcnoreach(1,:)','r--')
% xlim([0 12]) 
% ylim([0 1]) 

%%
myout=NaN(10,600);

% myout(1:10,1:size(myPool_cond1.survivalD1rate',2)) = myPool_cond1.survivalD1rate' ;
myout(1:6,1:size(myPool_cond1.survivalD1rate',2)) = myPool_cond1.survivalD1rate' ;

myout(1:10,101:100+size(myPool_cond2.survivalD1rate',2)) = myPool_cond2.survivalD1rate' ;
myout(1:10,201:200+size(myPool_cond3.survivalD1rate',2)) = myPool_cond3.survivalD1rate' ;
myout(1:10,301:300+size(myPool_cond4.survivalD1rate',2)) = myPool_cond4.survivalD1rate' ;
myout(1:10,401:400+size(myPool_cond5.survivalD1rate',2)) = myPool_cond5.survivalD1rate' ;
myout(1:10,501:500+size(myPool_cond6.survivalD1rate',2)) = myPool_cond6.survivalD1rate' ;

%% density
clear segm_den
segm_den = normdensity_CNT(:,3);
segm_den(isnan (segm_den)) = [];

den_distr_CNT = [];
den_distr_CNT (1,1) = sum (segm_den>1.13);
den_distr_CNT (1,2)  = den_distr_CNT (1,1)/size(segm_den,1);
den_distr_CNT (2,1) = sum (segm_den < 0.75);
den_distr_CNT (2,2)  = den_distr_CNT (2,1)/size(segm_den,1);
den_distr_CNT (3,1) = size(segm_den,1) - den_distr_CNT (1,1) - den_distr_CNT (2,1);
den_distr_CNT (3,2)  = den_distr_CNT (3,1)/size(segm_den,1);

clear segm_den
segm_den = normdensity_SRP(:,3);
segm_den(isnan (segm_den)) = [];

den_distr_SRP = [];
den_distr_SRP (1,1) = sum (segm_den>1.13);
den_distr_SRP (1,2)  = den_distr_SRP (1,1)/size(segm_den,1);
den_distr_SRP (2,1) = sum (segm_den < 0.75);
den_distr_SRP (2,2)  = den_distr_SRP (2,1)/size(segm_den,1);
den_distr_SRP (3,1) = size(segm_den,1) - den_distr_SRP (1,1) - den_distr_SRP (2,1);
den_distr_SRP (3,2)  = den_distr_SRP (3,1)/size(segm_den,1);

%% stability Vs initial AMPAR levels 
clear stable disap spines bl bl_stable bl_disap

%calcula baseline AMPAR levels for each spine
newStatPair=[{'utsmdodm'},{1},{'baselinech1'};...
             ];
         
numMaps=length(myMapList);

for i=1:numMaps
    
ps = mmMap.defaultPlotStruct();
% ps.plotbad = true;
% ps.plotintbad = true;
% ps.ploterrorwarning = true;
% ps.plotuserflag = true;

ps.stat = newStatPair{1}; 
ps.channel  = newStatPair{2};
ps = myMapList(i).GetMapValues(ps);
newStatName = newStatPair{3};
[m,n] = size(ps.val);
newStatValues = NaN(m,n);
% newStatValues = ps.val ./ mean(ps.val(~isnan(ps.val))); % ubssSum / mean(ubssSum);

 baseline = nanmean([myMapList(i).GetMapValuesCond(ps, 'base1') myMapList(i).GetMapValuesCond(ps, 'base2')],2);
% baseline =  myMapList(i).GetMapValuesCond(ps, 'base2');

newStatValues = repmat(baseline,1,n);


myMapList(i).addUserStat(newStatName, newStatValues);

% and then plot the new stat
% ps = myMapList(i).GetMapValues(ps);
% ps.stat = newStatName;
% myMapList(i).plotStat(ps);

% or use mmPlot class directly (most functions are static)
% mmPlot.plotStat(myMap, ps);

% and then save
%myMapList(i).save();

% now, load the map again and you will have your new stat

clearvars ps baseline n m newStatValues
end
clearvars s newStatPair

%Calculate baseline AMPAR levels for each spine group

cond1 = 'SRP';
cond2 = 'CNT';

SessionConditions = {'base1','base2','preSRP5','h48'};

StatName = 'baselinech1'; %Green signal values in second baseline

PersCond = 1; %Minimum number of sessions that spine needs to be present in
MinSpines = 10; %Minimum number of spines per dendrite for dendrite to be included

%VS
bl_cond1 = getmyMapListValuesCond(myMapList, mySegmentTable, cond1, SessionConditions, StatName, PersCond, MinSpines); %extract baseline values 

spines_cond1 =[];
spines_cond1 = myPool_cond1.values.Spines;

stable_cond1 = find(~isnan(spines_cond1(:,2))& ~isnan(spines_cond1(:,3))); %D1 spines that are maintained on D5
disap_cond1 = find(~isnan(spines_cond1(:,2))& isnan(spines_cond1(:,3))); %D1 spines that disappear before D5

bl_stable_cond1 = bl_cond1.values.Spines (stable_cond1,1);
bl_stable_cond1 (isnan (bl_stable_cond1))=[];
bl_disap_cond1 = bl_cond1.values.Spines (disap_cond1,1);
bl_disap_cond1 (isnan (bl_disap_cond1))=[];

[h,p] = kstest2(bl_stable_cond1, bl_disap_cond1)

 % histogram
% figure; histogram (VS(:,3), 50);

binRange = 0:0.05:2;
hcx = histcounts(bl_stable_cond1,[binRange Inf]);
hcy = histcounts(bl_disap_cond1,[binRange Inf]);
figure
bar(binRange,[hcx;hcy]')

%cumulative histogram
figure;
[h1,stats1] = cdfplot(bl_stable_cond1)
hold on
[h2.stats2] = cdfplot(bl_disap_cond1)

%CNT
%need to modify segment 24D2 in mySegmentTable because it does not contain
%spines with AMPAR values at baseline ('CNTNO')

bl_cond2 = getmyMapListValuesCond(myMapList, mySegmentTable, cond2, SessionConditions, StatName, PersCond, MinSpines); %extract baseline values 

spines_cond2 =[];
spines_cond2 = myPool_cond2.values.Spines;

stable_cond2 = find(~isnan(spines_cond2(:,2))& ~isnan(spines_cond2(:,3))); %D1 spines that are maintained on D5
disap_cond2 = find(~isnan(spines_cond2(:,2))& isnan(spines_cond2(:,3))); %D1 spines that disappear before D5

bl_stable_cond2 = bl_cond2.values.Spines (stable_cond2,1);
bl_stable_cond2 (isnan (bl_stable_cond2))=[];
bl_disap_cond2 = bl_cond2.values.Spines (disap_cond2,1);
bl_disap_cond2 (isnan (bl_disap_cond2))=[];

[h,p] = kstest2(bl_stable_cond2, bl_disap_cond2)

 % histogram
% figure; histogram (VS(:,3), 50);

binRange = 0:0.05:2;
hcx = histcounts(bl_stable_cond2,[binRange Inf]);
hcy = histcounts(bl_disap_cond2,[binRange Inf]);
figure
bar(binRange,[hcx;hcy]')

%cumulative histogram
figure;
[h1,stats1] = cdfplot(bl_stable_cond2)
hold on
[h2.stats2] = cdfplot(bl_disap_cond2)

% Compare VS and CNT
%eliminated spines
[h,p] = kstest2(bl_disap_cond1, bl_disap_cond2)

 % histogram
% figure; histogram (VS(:,3), 50);

binRange = 0:0.05:2;
hcx = histcounts(bl_disap_cond1,[binRange Inf]);
hcy = histcounts(bl_disap_cond2,[binRange Inf]);
figure
bar(binRange,[hcx;hcy]')

%cumulative histogram
figure;
[h1,stats1] = cdfplot(bl_disap_cond1)
hold on
[h2.stats2] = cdfplot(bl_disap_cond2)

%stable spines
[h,p] = kstest2(bl_stable_cond1, bl_stable_cond2)

 % histogram
% figure; histogram (VS(:,3), 50);

binRange = 0:0.05:2;
hcx = histcounts(bl_stable_cond1,[binRange Inf]);
hcy = histcounts(bl_stable_cond2,[binRange Inf]);
figure
bar(binRange,[hcx;hcy]')

%cumulative histogram
figure;
[h1,stats1] = cdfplot(bl_stable_cond1)
hold on
[h2.stats2] = cdfplot(bl_stable_cond2)