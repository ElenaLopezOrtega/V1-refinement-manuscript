%% Add new analysis to all maps in a pool for flag checking

newStatPair=[{'ubsdMean'},{2},{'ntbubsdMeanch2'}];

for s = 1:size(newStatPair,1)
         
numMaps=length(myMapList);

for i=1:numMaps
    
ps = mmMap.defaultPlotStruct();
ps.plotbad = false;
ps.plotintbad = false;
ps.ploterrorwarning = false;
ps.plotuserflag = true;

ps.stat = newStatPair{s,1}; 
ps.channel  = newStatPair{s,2};
ps = myMapList(i).GetMapValues(ps);
newStatName = newStatPair{s,3};
[m,n] = size(ps.val);
newStatValues = NaN(m,n);
% newStatValues = ps.val ./ mean(ps.val(~isnan(ps.val))); % ubssSum / mean(ubssSum);

baseline = nanmean([myMapList(i).GetMapValuesCond(ps, 'base1') myMapList(i).GetMapValuesCond(ps, 'base2')],2);

newStatValues = bsxfun(@rdivide, ps.val, baseline);


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

clear ps baseline n m newStatValues
end
end

clearvars s newStatPair

%% Add new analysis to all maps in a pool normalize to baseline per spine
% I am using mean mostly

newStatPair=[{'utsmdodm'},{1},{'ntbutsmdodmch1'};...
             {'utssdods'},{1},{'ntbutssdodsch1'};...
             {'utsmdtdm'},{2},{'ntbutsmdtdmch2'};...
             {'utssdtds'},{2},{'ntbutssdtdsch2'};...
             {'utdmdodm'},{1},{'ntbutdmdodmch1'};...
             {'utdsdods'},{1},{'ntbutdsdodsch1'};...
             {'uSEV'},{1},{'ntbuSEVch1'};...
             {'sLen2d'},{1},{'ntbsLen2d'};...
             {'sLen3d'},{1},{'ntbsLen3d'};...
             ];

for s = 1:size(newStatPair,1)
         
numMaps=length(myMapList);

for i=1:numMaps
    
ps = mmMap.defaultPlotStruct();
% ps.plotbad = true;
% ps.plotintbad = true;
% ps.ploterrorwarning = true;
% ps.plotuserflag = true;

ps.stat = newStatPair{s,1}; 
ps.channel  = newStatPair{s,2};
ps = myMapList(i).GetMapValues(ps);
newStatName = newStatPair{s,3};
[m,n] = size(ps.val);
newStatValues = NaN(m,n);
% newStatValues = ps.val ./ mean(ps.val(~isnan(ps.val))); % ubssSum / mean(ubssSum);

baseline = nanmean([myMapList(i).GetMapValuesCond(ps, 'base1') myMapList(i).GetMapValuesCond(ps, 'base2')],2);
%baseline =  myMapList(i).GetMapValuesCond(ps, 'base2');

newStatValues = bsxfun(@rdivide, ps.val, baseline);


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
end

clearvars s newStatPair

%% Add new analysis to all maps in a pool normalize to previous imaging session per spine
% I am using mean mostly

newStatPair=[{'utsmdodm'},{1},{'ntpsutsmdodmch1'}];

for s = 1:size(newStatPair,1)
         
numMaps=length(myMapList);

for i=1:numMaps
    
ps = mmMap.defaultPlotStruct();

ps.stat = newStatPair{s,1}; 
ps.channel  = newStatPair{s,2};
ps = myMapList(i).GetMapValues(ps);
newStatName = newStatPair{s,3};
[m,n] = size(ps.val);
newStatValues = NaN(m,n);

SessionConditions = myMapList(i).mapNV('condStr',:);

for k = 2:n
previoussessionname = SessionConditions{1,k-1}{1,1};
currentsessionname = SessionConditions{1,k}{1,1};
previoussession = [myMapList(i).GetMapValuesCond(ps, previoussessionname)]; 
currentsession = [myMapList(i).GetMapValuesCond(ps, currentsessionname)]; 

newStatValues(:,k) = currentsession./previoussession;

end



myMapList(i).addUserStat(newStatName, newStatValues);


clearvars ps baseline n m newStatValues
end
end

clearvars s newStatPair




%% Add new analysis to all maps in a pool normalize to segment baseline
% newStatPair=[{'utsmdodm'},{1},{'ntsbutsmdodmch1'};...
%              {'utssdods'},{1},{'ntsbutssdodsch1'}];

newStatPair=[{'utsmdodm'},{1},{'ntsbutsmdodmch1'};...
             {'utssdods'},{1},{'ntsbutssdodsch1'};...
             {'utsmdtdm'},{2},{'ntsbutsmdtdmch2'};...
             {'utssdtds'},{2},{'ntsbutssdtdsch2'};...
             {'utdmdodm'},{1},{'ntsbutdmdodmch1'};...
             {'utdsdods'},{1},{'ntsbutdsdodsch1'};...
             ];

for s = 1:size(newStatPair,1)
         
numMaps=length(myMapList);

for i=1:numMaps
    
ps = mmMap.defaultPlotStruct();
% ps.plotbad = true;
% ps.plotintbad = true;
% ps.ploterrorwarning = true;
% ps.plotuserflag = true;

ps.stat = newStatPair{s,1}; 
ps.channel  = newStatPair{s,2};
ps = myMapList(i).GetMapValues(ps);
[m,n] = size(ps.val);
newStatValues3D = NaN(m,n,myMapList(i).numMapSegments);
newStatValues = NaN(m,n);
newStatName = newStatPair{s,3};

for j = 1:myMapList(i).numMapSegments
ps.mapsegment=j;
ps = myMapList(i).GetMapValues(ps);
baseline = [myMapList(i).GetMapValuesCond(ps, 'base1') myMapList(i).GetMapValuesCond(ps, 'base2')];
baselineVec = reshape(baseline,[],1);
newStatValues3D(:,:,j) =  ps.val./nanmean(baselineVec);
end


newStatValues = nanmean(newStatValues3D,3);

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
end

clearvars s newStatPair

%% Remove alluser stats
numMaps=length(myMapList);
for i=1:numMaps
myMapList(i).removeUserStats();
end

disp('user stats removed')
clearvars i numMaps

%% Reset error flags user 3
numMaps=length(myMapList);
for i=1:numMaps
    
ps = mmMap.defaultPlotStruct();
ps.plotbad = true;
ps.plotintbad = true;
ps.ploterrorwarning = true;
ps.plotuserflag = true;
ps.stat = 'user3'; 
ps = myMapList(i).GetMapValues(ps);


[m,n] = size(ps.val);
newUser3 = NaN(m,n);

myMapList(i).toggleUser3(newUser3);
disp(['toggleUser3() map:' ps.mapName ]);

clearvars ps n m nnewUser3
end

disp('user3 reset')
%% Add error flag for certain spines
% SessionConditions = {'base1', 'base2', 'train1', 'train2', 'train3', 'train4', 'train5', 'train6', 'train7', 'train8', 't+1', 'ltm1', 'ltm2'};
% SessionConditions = {'base1','base2','postSRP1','postSRP2','preSRP3','postSRP3','postSRP4','preSRP5','postSRP5', 'h48','postSRP6','1week','postSRP7' };
%SessionConditions = {'base1','base2','preSRP3','preSRP5', 'h48','1week'};
SessionConditions = {'base1','base2','preSRP5', 'h48'};
BaselineConditions = {'base1', 'base2'};

numMaps=length(myMapList);

% Add flag in user3 for drastic changes in dendrite dsred fluorescence
% (exclude entire spinerun)
for i=1:numMaps
    
ps = mmMap.defaultPlotStruct();

ps.plotbad = false;
ps.plotintbad = false;
ps.ploterrorwarning = false;
ps.plotuserflag = true;

ps.stat = 'user3'; 
ps = myMapList(i).GetMapValues(ps);
newUser3 = ps.val;

ps.stat = 'ntbubsdMeanch2'; 
ps = myMapList(i).GetMapValues(ps);

[m,n] = size(ps.val);

% Applying criteria only on sessions that we want to analyze/plot later
currentconditions = cellstr(myMapList(i).GetValue_NV('condStr',1:n))'; % conditions of current map
condmembers = ismember(currentconditions,SessionConditions); % check which sessions are part of SessionConditions that we want to analyze/plot later

newUser3(repmat(min(ps.val(:,condmembers),[],2)<0.4,1,n))=1;
newUser3(repmat(max(ps.val(:,condmembers),[],2)>2.5,1,n))=1;



myMapList(i).toggleUser3(newUser3);
disp(['toggleUser3() map:' ps.mapName ]);

clear ps baseline n m newUser3 condmembers currentconditions
end

% Add flag in user3 for drastic changes in dendrite dsred fluorescence at
% baseline
% (exclude entire spinerun)
for i=1:numMaps
    
ps = mmMap.defaultPlotStruct();

ps.plotbad = false;
ps.plotintbad = false;
ps.ploterrorwarning = false;
ps.plotuserflag = true;

ps.stat = 'user3'; 
ps = myMapList(i).GetMapValues(ps);
newUser3 = ps.val;

ps.stat = 'ntbubsdMeanch2'; 
ps = myMapList(i).GetMapValues(ps);

[m,n] = size(ps.val);

% Applying criteria only on sessions that we want to analyze/plot later
currentconditions = cellstr(myMapList(i).GetValue_NV('condStr',1:n))'; % conditions of current map
condmembers = ismember(currentconditions,BaselineConditions); % check which sessions are part of SessionConditions that we want to analyze/plot later
tocheck = max(ps.val(:,condmembers),[],2) ./ min(ps.val(:,condmembers),[],2);
newUser3(repmat(tocheck > 2.5,1,n))=1;


myMapList(i).toggleUser3(newUser3);
disp(['toggleUser3() map:' ps.mapName ]);

clear ps baseline n m newUser3 condmembers currentconditions tockeck
end


% Add flag in user3 for low S/N (only exclude spine in session)
for i=1:numMaps
    
ps = mmMap.defaultPlotStruct();

ps.plotbad = true;
ps.plotintbad = true;
ps.ploterrorwarning = true;
ps.plotuserflag = true;

ps.stat = 'user3'; 
ps = myMapList(i).GetMapValues(ps);
newUser3 = ps.val;

[m,n] = size(ps.val);


ps.stat = 'sMean';
ps.channel = 2;
ps = myMapList(i).GetMapValues(ps);
currentsMean = ps.val;

ps.stat = 'sbMean';
ps.channel = 2;
ps = myMapList(i).GetMapValues(ps);
currentsbMean = ps.val;

ps.stat = 'sbSD';
ps.channel = 2;
ps = myMapList(i).GetMapValues(ps);
currentsbSD = ps.val;

currentSN = currentsMean ./ (currentsbMean + currentsbSD);

% Applying criteria only on sessions that we want to analyze/plot later
currentconditions = cellstr(myMapList(i).GetValue_NV('condStr',1:n))'; % conditions of current map
condmembers = ismember(currentconditions,SessionConditions); % check which sessions are part of SessionConditions that we want to analyze/plot later

newUser3(currentSN < 1) = 1;

myMapList(i).toggleUser3(newUser3);
disp(['toggleUser3() map:' ps.mapName ]);

clear ps baseline n m newUser3 condmembers currentconditions currentsbSD currentSN currentsMean currentsbMean
end


% Add flag in user3 for large differences in sN over dN (only exclude spine in session)
% Important for using sSum and dSum values
for i=1:numMaps
    
ps = mmMap.defaultPlotStruct();

ps.plotbad = true;
ps.plotintbad = true;
ps.ploterrorwarning = true;
ps.plotuserflag = true;

ps.stat = 'user3'; 
ps = myMapList(i).GetMapValues(ps);
newUser3 = ps.val;

[m,n] = size(ps.val);


ps.stat = 'sN';
ps.channel = 1;
ps = myMapList(i).GetMapValues(ps);
currentsN = ps.val;

ps.stat = 'dN';
ps.channel = 1;
ps = myMapList(i).GetMapValues(ps);
currentdN = ps.val;

currentsNdN = currentsN ./ currentdN;

% Applying criteria only on sessions that we want to analyze/plot later
currentconditions = cellstr(myMapList(i).GetValue_NV('condStr',1:n))'; % conditions of current map
condmembers = ismember(currentconditions,SessionConditions); % check which sessions are part of SessionConditions that we want to analyze/plot later

newUser3(currentsNdN > 10) = 1;

myMapList(i).toggleUser3(newUser3);
disp(['toggleUser3() map:' ps.mapName ]);

clear ps baseline n m newUser3 condmembers currentconditions currentsbSD currentSN currentsMean currentsbMean
end

clear currentdN currentsN currentsNdN i 

%% Add error flag for certain spines in one map NEEDS UPDATE SEE ABOVE
myMapList = myMapList;
MapName = 'rr107a';
i=find(strcmp({myMapList(:).mapName},MapName));
SessionConditions = {'base1', 'base2', 'train1', 'train2', 'train3', 'train4', 'train5', 'train6', 'train7', 'train8', 't+1', 'ltm1', 'ltm2'};

ps = mmMap.defaultPlotStruct();
ps.stat = 'ntbubsdMeanch2'; 
ps.channel  = 2;

ps.plotbad = true;
ps.plotintbad = true;
ps.ploterrorwarning = true;
ps.plotuserflag = true;
ps = myMapList(i).GetMapValues(ps);
 


[m,n] = size(ps.val);
newUser3 = NaN(m,n);

% Applying criteria only on sessions that we want to analyze/plot later
currentconditions = cellstr(myMapList(i).GetValue_NV('condStr',1:n))'; % conditions of current map
condmembers = ismember(currentconditions,SessionConditions); % check which sessions are part of SessionConditions that we want to analyze/plot later

newUser3(repmat(min(ps.val(:,condmembers),[],2)<0.4,1,n))=1;
newUser3(repmat(max(ps.val(:,condmembers),[],2)>2.5,1,n))=1;

myMapList(i).toggleUser3(newUser3);

clear ps baseline n m nnewUser3



%% Plot an individual segment

myMapList = myMapList;
MapName = 'V1_030_E';
SegmentNumber = 2;
PersCond = 1; %Minimum number of sessions that spine needs to be present in
%SessionConditions = {'base1','base2','postSRP1','postSRP2','preSRP3','postSRP3','postSRP4','preSRP5','postSRP5', 'h48','postSRP6','1week','postSRP7'};
%SessionConditions = {'base1','base2','postSRP1','preSRP5','postSRP5', 'h48','postSRP6','1week','postSRP7'};
% SessionConditions = {'base1','base2','postSRP3','postSRP5', 'h48','1week'};
SessionConditions = {'base1','base2','postSRP5', 'h48','1week'};
%StatName = 'ntsbutsmdodmch1'; %background subtracted dendrite sum;
StatName = 'ntbutsmdodmch1'; %Spine AMPA / Dendrite dsRed mean
MinSpines = 10;
%StatName = 'user3'; %background subtracted dendrite sum;

clearvars values
values = segmentCall(myMapList, MapName, SegmentNumber, SessionConditions, StatName, PersCond, MinSpines);
figure()
plot (values.mean)

%% Plot all segments of one condition
 
 Condition = 'SRP'; %V1
% Condition = 'S102KOSt'; %SAPs

% V1 analysis
   SessionConditions = {'base1','base2','postSRP1','preSRP5','postSRP5','h48'};
%   SessionConditions = {'base1','base2','postSRP1','preSRP3','preSRP5','postSRP5','h48'};
%   SessionConditions = {'base1','base2','preSRP5','h48'};

%  SessionConditions = {'base1','base2','preSRP3','preSRP5', 'h48', '1week'};
%SessionConditions = {'base1','base2','postSRP1','postSRP2','preSRP3','postSRP3','postSRP4','preSRP5','postSRP5', 'h48','postSRP6','1week','postSRP7'};
%SessionConditions = {'base1','base2','preSRP3','preSRP5', 'h48' };
%SessionConditions = {'base1','base2','postSRP1','postSRP2','postSRP3','postSRP4','postSRP5', 'h48' };

% SAPs analysis
% SessionConditions = {'base1','base2','hhalf','h1','h2','h3','h24','h48'};

% Normalized to spine baseline
    StatName = 'ntbutsmdodmch1'; %Spine AMPA / Dendrite dsRed mean
%    StatName = 'ntbutdmdodmch1'; %Dendrite AMPA / Dendrite dsRed mean
%    StatName = 'ntbutsmdtdmch2'; %Spine dsRed / Dendrite dsRed mean

% Normalized to segment baseline
%   StatName = 'ntsbutsmdodmch1'; %Spine AMPA / Dendrite dsRed mean
%   StatName = 'ntsbutdmdodmch1'; %Dendrite AMPA / Dendrite dsRed mean
%   StatName = 'ntsbutsmdtdmch2'; %Spine dsRed / Dendrite dsRed mean

% StatName = 'ntpsutsmdodmch1'; 

% StatName = 'pDist';
% StatName = 'ntbuSEVch1';
% StatName = 'ntsbutsmdodmch1'; 
% StatName = 'ntbubsdMeanch2';

PersCond = 3; %Minimum number of sessions that spine needs to be present in
MinSpines = 10; %Minimum number of spines per dendrite for dendrite to be included

myPool = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);


% Plot results
  plotMean = 'SegmentsSpines';
%  plotMean = 'Spines';

 plotMyPool(myPoolSRPCN, plotMean) %V1
% plotMyPool(myPool_WTCNT_NS, plotMean) %SAPs

%%
% V1 analysis
% Condition = 'SRPV1CN';
Condition = 'SRP';
myPoolSRPCN = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);

Condition = 'CNT';
myPoolCNTCN = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);

% % SAPs analysis
% Condition = 'S102KOSt';
% myPool_S102KO_St = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);
% 
% Condition = 'S102KONS';
% myPool_S102KO_NS = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);
% 
% Condition = 'S102CNTSt';
% myPool_S102CNT_St = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);
% 
% Condition = 'S102CNTNS';
% myPool_S102CNT_NS = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);
% 
% Condition = 'WTCNTSt';
% myPool_WTCNT_St = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);
% 
% Condition = 'WTCNTNS';
% myPool_WTCNT_NS = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);

% Condition = 'S97KOSt';
% myPool_S97KO_St = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);
% 
% Condition = 'S97KONS';
% myPool_S97KO_NS = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);
% 
% Condition = 'S97CNTSt';
% myPool_S97CNT_St = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);
% 
% Condition = 'S97CNTNS';
% myPool_S97CNT_NS = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);

%More V1 conditions
%  Condition = 'SRPV1ED';
%  myPoolSRPED = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);
% 
% Condition = 'CNTV1ED';
% myPoolCNTED = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);
% 
% Condition = 'SRPVCtx';
% myPoolSRPVCtx = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);

%Condition = 'CNTVCtx';
%myPoolCNTVCtx = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);

% Condition = 'SRPCntr';
% myPoolSRPCntr = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);
% 
% Condition = 'CNTCntr';
% myPoolCNTCntr = getmyMapListValuesCond(myMapList, mySegmentTable, Condition, SessionConditions, StatName, PersCond, MinSpines);


%% Plot one stat and pDist
clearvars myPoolStat_cond1 myPoolSecStat_cond1 myPoolStat_cond2 myPoolSecStat_cond2
clearvars myPoolStat_cond3 myPoolSecStat_cond3 myPoolStat_cond4 myPoolSecStat_cond4

% Conditions for all segments
% Condition = 'train';
  Condition1 = 'CNT'; 
  Condition2 = 'SRP';
 
% Conditions for specific segments (load mySegmentTable_022521_ss)
%  Condition1 = 'SRP1'; % SRP segments that decrease spine density and increase sAMPAR levels
%  Condition2 = 'SRP2'; % SRP segments that maintain spine density and increase sAMPAR levels
%  Condition3 = 'SRP3'; % SRP segments that decrease spine density and maintain sAMPAR levels
%  Condition4 = 'SRP4'; % SRP segments that maintain spine density and maintain sAMPAR levels

% SessionConditions = {'base1', 'base2', 'train1', 'train2', 'train3', 'train4', 'train5', 'train6', 'train7', 'train8'};
 SessionConditions = {'base1','base2','postSRP1','preSRP5','postSRP5','h48'};
% SessionConditions = {'base1','base2','preSRP5','h48'};
N_sess= size (SessionConditions,2);
%SessionConditions = {'train'};

StatName = 'ntbutsmdodmch1'; %Spine AMPA / Dendrite dsRed mean
%StatName = 'ntbutdmdodmch1'; %Dendrite AMPA / Dendrite dsRed mean
%StatName = 'z';
%SecStatName = 'ntsbutsmdodmch1';

SecStatName = 'pDist';
% SecStatName = 'ntbutdmdodmch1'; %Dendrite AMPA / Dendrite dsRed mean
%SecStatName = 'ntbutsmdtdmch2'; %Spine dsRed / Dendrite dsRed mean
%SecStatName = 'z';

PersCond = 3; %Minimum number of sessions that spine needs to be present in
MinSpines = 10; %Minimum number of spines per dendrite for dendrite to be included 
% PersCond = 1; %Minimum number of sessions that spine needs to be present in
% MinSpines = 1; %Minimum number of spines per dendrite for dendrite to be included 

%Only spines that have values in main Stat are considered in  SecStat
myPoolStat_cond1 = getmyMapListValuesCond(myMapList, mySegmentTable, Condition1, SessionConditions, StatName, PersCond, MinSpines);
myPoolSecStat_cond1 = getmyMapListValuesCondDist_v2(myMapList, mySegmentTable, Condition1, SessionConditions, StatName, SecStatName, PersCond, MinSpines);
[segment_cond1] = plotincrdecr_onsegment_v2 (myPoolStat_cond1, myPoolSecStat_cond1);

myPoolStat_cond2 = getmyMapListValuesCond(myMapList, mySegmentTable, Condition2, SessionConditions, StatName, PersCond, MinSpines);
myPoolSecStat_cond2 = getmyMapListValuesCondDist_v2(myMapList, mySegmentTable, Condition2, SessionConditions, StatName, SecStatName, PersCond, MinSpines);
[segment_cond2] = plotincrdecr_onsegment_v2 (myPoolStat_cond2, myPoolSecStat_cond2);

% myPoolStat_cond3 = getmyMapListValuesCond(myMapList, mySegmentTable, Condition3, SessionConditions, StatName, PersCond, MinSpines);
% myPoolSecStat_cond3 = getmyMapListValuesCondDist_v2(myMapList, mySegmentTable, Condition3, SessionConditions, StatName, SecStatName, PersCond, MinSpines);
% [segment_cond3] = plotincrdecr_onsegment_v2 (myPoolStat_cond3, myPoolSecStat_cond3);
% 
% myPoolStat_cond4 = getmyMapListValuesCond(myMapList, mySegmentTable, Condition4, SessionConditions, StatName, PersCond, MinSpines);
% myPoolSecStat_cond4 = getmyMapListValuesCondDist_v2(myMapList, mySegmentTable, Condition4, SessionConditions, StatName, SecStatName, PersCond, MinSpines);
% [segment_cond4] = plotincrdecr_onsegment_v2 (myPoolStat_cond4, myPoolSecStat_cond4);


prop_inc_cond1= segment_cond1(:,1)./segment_cond1(:,3);
prop_inc_cond2 = segment_cond2(:,1)./segment_cond2(:,3);
[h,p] = ttest2 (prop_inc_cond1,prop_inc_cond2);

prop_dec_cond1= segment_cond1(:,2)./segment_cond1(:,3);
prop_dec_cond2 = segment_cond2(:,2)./segment_cond2(:,3);
[h,p] = ttest2 (prop_dec_cond1,prop_dec_cond2);

prop_stb_cond1 = (segment_cond1(:,3) - segment_cond1(:,1) - segment_cond1(:,2))./segment_cond1(:,3);
prop_stb_cond2 = (segment_cond2(:,3) - segment_cond2(:,1) - segment_cond2(:,2))./segment_cond2(:,3);
[h,p] = ttest2 (prop_stb_cond1,prop_stb_cond2);

%% Proportion of spines that change in the same direction per segment

distrib_cond1 = abs ((segment_cond1(:,1) - segment_cond1(:,2))./(segment_cond1(:,1) + segment_cond1(:,2)));
distrib_cond1(isnan(distrib_cond1)) = 0;

distrib_cond2 = abs ((segment_cond2(:,1) - segment_cond2(:,2))./(segment_cond2(:,1) + segment_cond2(:,2)));
distrib_cond2(isnan(distrib_cond2)) = 0;

% histogram

binRange = 0:0.1:1;
hcx = histcounts(distrib_cond1,[binRange Inf]);
hcy = histcounts(distrib_cond2,[binRange Inf]);
figure
bar(binRange,[hcx;hcy]')

binRange = 0:0.1:1;
hcx = histcounts(distrib_cond1,[binRange Inf], 'Normalization', 'probability');
hcy = histcounts(distrib_cond2,[binRange Inf], 'Normalization', 'probability');
figure
bar(binRange,[hcx;hcy]', 1)
ylim ([0 1])

[h1,p1] = kstest2(distrib_cond1,distrib_cond2);

%randomization

%cond1

spines_CNT = myPoolStat_cond1.values.Spines (:,4);
spines_CNT (isnan (spines_CNT)) = [];

clear idx shuffle

numofrand = 10000;
sp_per_segm = segment_cond1 (:,3);
segment_cond1_perm = [];
shuffle_distrib_cond1 = [];

for b=1:numofrand
idx = randperm(size(spines_CNT,1));
spines_CNT_perm(idx,1)= spines_CNT;

    i=1;
    for a = 1:size (segment_cond1,1)
        f = sp_per_segm (a);
        segm = spines_CNT_perm (i:i+f-1);
        i = f+1;
        segment_cond1_perm (a,1) =  sum(segm > 1.43);
        segment_cond1_perm (a,2) =  sum(segm < 0.52);
        segment_cond1_perm (a,3) =  size(segm,1);
    end
        

shuffle = abs ((segment_cond1_perm(:,1) - segment_cond1_perm(:,2))./(segment_cond1_perm(:,1) + segment_cond1_perm(:,2)));
shuffle (isnan(shuffle)) = 0;
shuffle_distrib_cond1= [shuffle_distrib_cond1 shuffle]; 
end

%  shuffle_distrib_cond1 = mean (shuffle_distrib_cond1,2);
% shuffle_distrib = shuffle_distrib_cond1 (:,1000);
% shuffle_distrib_cond1 = [];
% shuffle_distrib_cond1 = shuffle_distrib;

% histogram

binRange = 0:0.1:1;
hcx = histcounts(distrib_cond1,[binRange Inf]);
hcy = histcounts(shuffle_distrib_cond1,[binRange Inf]);
figure
bar(binRange,[hcx;hcy]')

binRange = 0:0.1:1;
hcx = histcounts(shuffle_distrib_cond1,[binRange Inf], 'Normalization', 'probability');
hcy = histcounts(distrib_cond1,[binRange Inf], 'Normalization', 'probability');
figure
bar(binRange,[hcx;hcy]')
ylim ([0 1])

shuffle_distrib_cond1 = shuffle_distrib_cond1(:);
[h2,p2] = kstest2(distrib_cond1,shuffle_distrib_cond1);

%cond2

spines_SRP = myPoolStat_cond2.values.Spines (:,4);
spines_SRP (isnan (spines_SRP)) = [];

clear idx shuffle a f i sp_per_segm segm

numofrand = 10000;
sp_per_segm = segment_cond2 (:,3);
segment_cond2_perm = [];
shuffle_distrib_cond2 = [];

for b=1:numofrand
idx = randperm(size(spines_SRP,1));
spines_SRP_perm(idx,1)= spines_SRP;

    i=1;
    for a = 1:size (segment_cond2,1)
        f = sp_per_segm (a);
        segm = spines_SRP_perm (i:i+f-1);
        i = f+1;
        segment_cond2_perm (a,1) =  sum(segm > 1.43);
        segment_cond2_perm (a,2) =  sum(segm < 0.52);
        segment_cond2_perm (a,3) =  size(segm,1);
    end
        

shuffle = abs ((segment_cond2_perm(:,1) - segment_cond2_perm(:,2))./(segment_cond2_perm(:,1) + segment_cond2_perm(:,2)));
shuffle (isnan(shuffle)) = 0;
shuffle_distrib_cond2= [shuffle_distrib_cond2 shuffle]; 
end

% shuffle_distrib_cond2 = mean (shuffle_distrib_cond2,2);

% histogram

binRange = 0:0.1:1;
hcx = histcounts(distrib_cond2,[binRange Inf]);
hcy = histcounts(shuffle_distrib_cond2,[binRange Inf]);
figure
bar(binRange,[hcx;hcy]')

binRange = 0:0.1:1;
hcx = histcounts(shuffle_distrib_cond2,[binRange Inf], 'Normalization', 'probability');
hcy = histcounts(distrib_cond2,[binRange Inf], 'Normalization', 'probability');
figure
bar(binRange,[hcx;hcy]')
ylim ([0 1])

shuffle_distrib_cond2 = shuffle_distrib_cond2(:);
[h3,p3] = kstest2(distrib_cond2,shuffle_distrib_cond2);


% clearvars val xval yval valInc xvalInc yvalInc
% 
% xval = nanmean(myPoolSecStat.values.Spines(:,1:N_sess),2);
% yval = nanmean(myPool.values.Spines(:,3:N_sess),2);
% 
% % xval = nanmean(myPoolSecStat.values.Spines(:,1:10),2);
% % yval = nanmean(myPool.values.Spines(:,1:10),2);
% 
% % yval = nanmean(myPool.values.Spines(:,9:10),2);
% % xval = nanmean(myPool.values.Spines(:,3:4),2);
% 
% % xval = nanmean(myPoolSecStat.values.Spines(:,1:1),2);
% % yval = nanmean(myPool.values.Spines(:,1:1),2);
% 
% 
% xval = xval(~isnan(yval));
% yval = yval(~isnan(yval));
% 
% % yval = yval(~isnan(xval));
% % xval = xval(~isnan(xval));
% 
% % yvalInc = yval(yval < 1);
% % xvalInc = xval(yval < 1);
% 
% yvalInc = yval(xval > 1);
% xvalInc = xval(xval > 1);
% 
% figure
% plot(xval,yval,'o')
% 
% val(:,1)=xval;
% val(:,2)=yval;
% 
% valInc=[];
% valInc(:,1)=xvalInc;
% valInc(:,2)=yvalInc;

%% Create model data set
list_SRP = [];
list_CNT = [];

dataset_SRP = spine_dist_segm_dataset (myPoolStat_cond2, myPoolSecStat_cond2, Condition2, list_SRP);
dataset_CNT = spine_dist_segm_dataset (myPoolStat_cond1, myPoolSecStat_cond1, Condition1, list_CNT);

dataset_T = vertcat (dataset_SRP, dataset_CNT);

T = table(dataset_T(:,1), dataset_T(:,2), dataset_T(:,3), dataset_T(:,4), dataset_T(:,5), dataset_T(:,6), dataset_T(:,7), 'VariableNames', { 'spine', 'position', 'change', 'dendrite', 'cell', 'mouse', 'condition'} );
% Write data to text file
writetable(T, 'dataset.csv')
%% Plot SADS changes along dendrite and plot changes in 5um bins

randomize = 0;
toplotAll = {};
figure
allFieldNames = fieldnames(myPool.values.SpinePerSegment);
p=1;
for rows = 1:length(allFieldNames)
%for rows = 50
    
current_mapname = allFieldNames{rows};
[myPool.values.SpinePerSegmentNaNRemoved.(allFieldNames{rows}),currIdxs] = removeNanRows(myPool.values.SpinePerSegment.(allFieldNames{rows}));
myPoolSecStat.values.SpinePerSegmentNaNRemoved.(allFieldNames{rows}) = myPoolSecStat.values.SpinePerSegment.(allFieldNames{rows})(currIdxs,:);

currentpDist = nanmean(myPoolSecStat.values.SpinePerSegmentNaNRemoved.(allFieldNames{rows})(:,3:10),2);
currentSADS = nanmean(myPool.values.SpinePerSegmentNaNRemoved.(allFieldNames{rows})(:,3:10),2);
currentSADSincr = currentSADS;
currentSADSincr(currentSADS < 1.3) = NaN;
currentSADSdecr = currentSADS;
currentSADSdecr(currentSADS > 0.7) = NaN;
% randomize currentpDist
if randomize == 1
currentpDist = currentpDist(randperm(size(currentpDist,1)));
end

toplot=[currentpDist currentSADS currentSADSincr currentSADSdecr];
toplot_sorted=sortrows(toplot,1);

toplot_norm(:,1) = [-100:5:100]; 
for l = 2:length(toplot_norm(:,1))

    toplot_norm(l,2)=nanmean(toplot_sorted(toplot_sorted(:,1)<toplot_norm(l,1) & toplot_sorted(:,1)>toplot_norm(l-1,1),2));
end

if size(toplot_sorted,1) <= 3
else

[acf,lags]=autocorr(toplot_sorted(:,2));
autocorrfctn(1:max(lags)+1,p)=acf;
autocorrfctnlgd(1,p)={current_mapname};
subplot(16,4,p)
%plot(toplot_sorted(:,1),toplot_sorted(:,2),'LineWidth',2)
scatter(toplot_sorted(:,1),toplot_sorted(:,3),10,'filled','r')
hold on
scatter(toplot_sorted(:,1),toplot_sorted(:,4),10,'filled','b')
 hold on
 plot([min(toplot_sorted(:,1)),max(toplot_sorted(:,1))],[1,1],'k-', 'LineWidth',2)

ylim([-0 2])
xlim([-100 100])
%title(current_mapname)
hold on

% segment(p,1)=length(current_incr_pos);
% segment(p,2)=length(current_decr_pos);
% segment(p,3)=length(current_allspine_pos);


% SADS=[SADS;currentSADS];
% pdist=[pdist;currentpdist];
% incr_decr_ratio=[incr_decr_ratio;current_incr_decr_ratio];
toplotNormAll(:,p)=toplot_norm(:,2);
toplotAll{p,1}=toplot_sorted;
p=p+1;
%axis off


end
end
nanmean(nanmean(autocorrfctn(2:7,:)))
figure
pcolor(toplotNormAll')
shading flat;
set(gca, 'ydir', 'reverse');
colormap jet;
caxis([0.0 2]);
colorbar;
%% Plot median SADS around a spine

allFieldNames = fieldnames(myPool.values.SpinePerSegment);

randomize = 1;
repeats = 100;
incordec = 3;

clearvars binsY Xall Yall

for rep = 1:repeats
inccentplotall_pdist=[];
inccentplotall_SADS=[];
nnball3D={};
for rows = 1:length(allFieldNames)
%for rows = 50   
current_mapname = allFieldNames{rows};
[myPool.values.SpinePerSegmentNaNRemoved.(allFieldNames{rows}),currIdxs] = removeNanRows(myPool.values.SpinePerSegment.(allFieldNames{rows}));
myPoolSecStat.values.SpinePerSegmentNaNRemoved.(allFieldNames{rows}) = myPoolSecStat.values.SpinePerSegment.(allFieldNames{rows})(currIdxs,:);


% figure
% x = nanmean(myPoolpDist.values.SpinePerSegmentNaNRemoved.(allFieldNames{rows}),2);
% y = nanmean(myPool.values.SpinePerSegmentNaNRemoved.(allFieldNames{rows})(:,3:10),2);
% scatter(x,y)

currentpDist = nanmean(myPoolSecStat.values.SpinePerSegmentNaNRemoved.(allFieldNames{rows})(:,3:10),2);
currentSADS = nanmean(myPool.values.SpinePerSegmentNaNRemoved.(allFieldNames{rows})(:,3:10),2);

% currentSADS = currentSADS./nanmean(currentSADS); %Normalize to average change of dendrite???

currentSADSincr = currentSADS;
currentSADSincr(currentSADS < 1.3) = NaN;
currentSADSdecr = currentSADS;
currentSADSdecr(currentSADS > 1.0) = NaN;
% randomize currentpDist
if randomize == 1
currentpDist = currentpDist(randperm(size(currentpDist,1)));
end

toplot=[currentpDist currentSADS currentSADSincr currentSADSdecr];
toplot_sorted=sortrows(toplot,1);


if size(toplot_sorted,1) <= 3
else
    toplot_padded=padarray(toplot_sorted,5,NaN,'both');
    
inccentplot_pdist=[];    
inccentplot_SADS=[];    
col=1;
for row=6:size(toplot_sorted,1)+5
    
if isnan(toplot_padded(row,incordec)) == 0
    %row=5+randi(size(toplot_sorted,1));
    
    inccentplot_pdist(1:11,col)=toplot_padded(row-5:row+5,1)-toplot_padded(row,1);
    inccentplot_SADS(1:11,col)=toplot_padded(row-5:row+5,2);
    
    col=col+1;
    
end    
end

inccentplotall_pdist=[inccentplotall_pdist inccentplot_pdist];
inccentplotall_SADS=[inccentplotall_SADS inccentplot_SADS];

% t=size(nnball3D,1);
% inccentplot3D_pdist{t+1,1}=inccentplot_pdist;
% inccentplot3D_SADS{t+1,1}=inccentplot_SADS;
% inccentplot3D_pdist{t+1,2}=current_mapname;
% inccentplot3D_SADS{t+1,2}=current_mapname;

end



end

Xall(:,rep)=reshape(inccentplotall_pdist,[11*size(inccentplotall_pdist,2) 1]);
Yall(:,rep)=reshape(inccentplotall_SADS,[11*size(inccentplotall_pdist,2) 1]);

% inccentconc=[X Y];


% %%
% figure
% plot(-5:5, nanmean(inccentplotall_SADS,2), 'b-o')
% hold on
% plot(-5:5, ones(11), 'k--')
% 
% xlim([-6 6])
% ylim([0.8 1.8])
% 
% 
% figure
% for i=1:size(inccentplotall_SADS,2)
%     
% scatter(inccentplotall_pdist(:,i), inccentplotall_SADS(:,i))
% hold on    
% end




% X=nanmean(X,2);
% Y=nanmean(Y,2);

X = Xall(:,rep);
Y = Yall(:,rep);

X=abs(X);



%bins=[-47:2:59];
binmin=0;
binint=2;
binmax=59;

bins=[binmin:binint:binmax];
bincount = length(bins);
j=1;
binsY(j,1,rep)=bins(j);
%binsY(j,2)=median(Y((X < bins(j)),1));
binsY(j,2,rep)=median(Y((X ==0),1),'omitnan');
%binsY(j,3)=sum((X < bins(j)),1);
binsY(j,3,rep)=sum((X == 0),1);
binsY(j,4,rep)=std(Y((X < bins(j))),'omitnan');
binsY(j,5,rep)=binsY(j,4)./sqrt(binsY(j,3));


for j=2:bincount


binsY(j,1,rep)=bins(j)-1;
binsY(j,2,rep)=median(Y((X < bins(j) & X > bins(j-1)),1),'omitnan');
binsY(j,3,rep)=sum((X < bins(j) & X > bins(j-1)),1);
binsY(j,4,rep)=std(Y((X < bins(j) & X > bins(j-1))),'omitnan');
binsY(j,5,rep)=binsY(j,4)./sqrt(binsY(j,3));

end


end
binsY=nanmean(binsY,3);

figure
errorbar(binsY(:,1),binsY(:,2),binsY(:,5) , 'b-o')
hold on
plot([-1 50], ones(2), 'k--')

xlim([-1 50])
ylim([0.8 1.5])
binsY(1,2:3)