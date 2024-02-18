function retPool = getmyMapListValuesCondDyn_v2(myMapList, mySegmentTable, Condition, SessionConditions, StatName, channel, PersCond, MinSpines)

% Plot all segments of one condition
% 
% Condition = 'train';
% SessionConditions = {'base1', 'base2', 'train1', 'train2', 'train3', 'train4', 'train5', 'train6', 'train7', 'train8', 't+1', 'ltm1', 'ltm2'};
% StatName = 'ntbutsmdodmch1'; %background subtracted dendrite sum
% PersCond = 6; %Minimum number of sessions that spine needs to be present in
% MinSpines = 5;


% only include segments that are not excluded
rows=mySegmentTable.Include == 1;
myGoodSegmentTable = mySegmentTable(rows,:);
clearvars rows

% only include segments that are part of condition to plot
rows=strcmp(myGoodSegmentTable.MouseCondition,Condition);
myCondSegmentTable = myGoodSegmentTable(rows,:);
clearvars SegmentNumber i rows values valuesgeo meanvalues stdev val valneurons valneuronsSegmentsmeans valmice valmiceSegmentsmeans valmiceNeuronsmeans valmiceNeuronsSegmentsmeans valneuronsSegmentsgeomeans
clearvars valmiceSegmentsgeomeans valmiceNeuronsSegmentsgeomeans valmiceNeuronsgeomeans se retPool NeuronNames NeuronName geomeanvalues currentvalues currentmousename count MapName MouseName

%val=struct('values',[],'mean',[],'geomean',[],'count',[],'std',[],'se',[],'segment',[],'addpercent',[],'subpercent',[],'survival',[],'survivalD1',[],'countall',[]);
val=struct('values',[],'mean',[],'geomean',[],'count',[],'std',[],'se',[],'segment',[],'addpercent',[],'subpercent',[],'survival',[],'survivalD1',[],'countall',[], 'dist',[]); %added by ELO on 2/27/21
valneurons=struct();
valneuronsSegmentsmeans=struct();
valneuronsSegmentsgeomeans=struct();
valmice = struct();
valmiceSegmentsmeans = struct();
valmiceSegmentsgeomeans = struct();
valmiceNeuronsmeans = struct();
valmiceNeuronsgeomeans = struct();
valmiceNeuronsSegmentsmeans = struct();
valmiceNeuronsSegmentsgeomeans = struct();

for i = 1:size(myCondSegmentTable,1)
   
MapName = myCondSegmentTable.MapName{i};
MouseName = myCondSegmentTable.MouseName{i};
% NeuronName = myCondSegmentTable.neuron{i}; %I dont have info about neuron
SegmentNumber = myCondSegmentTable.SegmentNumber(i);

currentvalues = segmentCallDyn_v2(myMapList, MapName, SegmentNumber, SessionConditions, StatName, channel, PersCond, MinSpines);

% % Use this bit of code to normalize currentvaluesPDist.values to 0-1
% currentvalues.values = currentvalues.values - min(currentvalues.values);
% currentvalues.values = currentvalues.values ./ max(currentvalues.values);


val.values = [val.values; currentvalues.values];
val.valuesSpinePerSegment.([MapName num2str(SegmentNumber)]) = currentvalues.values;
val.mean = [val.mean; currentvalues.mean];
val.geomean = [val.geomean; currentvalues.geomean];
val.count = [val.count; currentvalues.count];
val.std = [val.std; currentvalues.std];
val.se = [val.se; currentvalues.se];
val.segment = [val.segment; {[MapName '_' num2str(SegmentNumber)]}];

val.addpercent = [val.addpercent; currentvalues.addpercent];
val.subpercent = [val.subpercent; currentvalues.subpercent];
val.survival = [val.survival; currentvalues.survival];
val.survivalD1 = [val.survivalD1; currentvalues.survivalD1];
val.countall = [val.countall; currentvalues.countall];
val.dist = [val.dist; currentvalues.dist]; %added by ELO on 2/27/21

% if isempty(NeuronName) == 0 % skip if NeuronName is empty (as in not defined in mySegmentTable)
%     if isfield (valneurons,NeuronName) == 0 % make new field in valneuron structure if it's first time this neuron comes up
%         valneurons.(NeuronName) = [];
%         valneuronsSegmentsmeans.(NeuronName) = [];
%         valneuronsSegmentsgeomeans.(NeuronName) = [];
%         valneurons.(NeuronName) = [valneurons.(NeuronName); currentvalues.values]; % concatenate values to this new field
%         valneuronsSegmentsmeans.(NeuronName) = [valneuronsSegmentsmeans.(NeuronName); currentvalues.mean]; % concatenate values to this new field
%         valneuronsSegmentsgeomeans.(NeuronName) =  [valneuronsSegmentsgeomeans.(NeuronName); currentvalues.geomean]; % concatenate values to this new field
%     else 
%         valneurons.(NeuronName) = [valneurons.(NeuronName); currentvalues.values]; % concatenate values to existing field
%         valneuronsSegmentsmeans.(NeuronName) = [valneuronsSegmentsmeans.(NeuronName); currentvalues.mean]; % concatenate values to this new field
%         valneuronsSegmentsgeomeans.(NeuronName) =  [valneuronsSegmentsgeomeans.(NeuronName); currentvalues.geomean]; % concatenate values to this new field
%     end
% end

if isempty(MouseName) == 0 % skip if MouseName is empty (as in not defined in mySegmentTable)
    if isfield (valmice,MouseName) == 0 % make new field in valmice structure if it's first time this mouse comes up
        valmice.(MouseName) = [];
        valmiceSegmentsmeans.(MouseName) = [];
        valmiceSegmentsgeomeans.(MouseName) = [];
        valmice.(MouseName) = [valmice.(MouseName); currentvalues.values]; % concatenate values to this new field
        valmiceSegmentsmeans.(MouseName) = [valmiceSegmentsmeans.(MouseName); currentvalues.mean]; % concatenate values to this new field
        valmiceSegmentsgeomeans.(MouseName) = [valmiceSegmentsgeomeans.(MouseName); currentvalues.geomean]; % concatenate values to this new field
    else 
        valmice.(MouseName) = [valmice.(MouseName); currentvalues.values]; % concatenate values to existing field
        valmiceSegmentsmeans.(MouseName) = [valmiceSegmentsmeans.(MouseName); currentvalues.mean]; % concatenate values to this new field
        valmiceSegmentsgeomeans.(MouseName) = [valmiceSegmentsgeomeans.(MouseName); currentvalues.geomean]; % concatenate values to this new field
    end
end

end



% Values at spine level
values.Spines = val.values;
values.SpinePerSegment = val.valuesSpinePerSegment;
meanvalues.Spines = nanmean(val.values,1);
geomeanvalues.Spines = nangeomean(val.values,1);
count.Spines = nansum(val.count,1);
stdev.Spines = std(values.Spines,0,1,'omitnan');
se.Spines = stdev.Spines ./ sqrt(count.Spines);

% Values at segment level
values.SegmentsSpines = val.mean;
meanvalues.SegmentsSpines = nanmean(val.mean,1);
geomeanvalues.SegmentsSpines = nangeomean(val.geomean,1);
count.SegmentsSpines = sum(~isnan(val.mean),1);
stdev.SegmentsSpines = std(values.SegmentsSpines,0,1,'omitnan');
se.SegmentsSpines = stdev.SegmentsSpines ./ sqrt(count.SegmentsSpines);

additions =  val.addpercent;
eliminations = val.subpercent;
survivalrate = val.survival;
survivalD1rate = val.survivalD1;
countall = val.countall;
dist = val.dist; %added by ELO on 2/27/21

% Values at neuron level
% values.NeuronsSpines = structfun(@(x) nanmean(x,1),valneurons,'UniformOutput',false);
% meanvalues.NeuronsSpines = nanmean(cell2mat(struct2cell(values.NeuronsSpines)),1);
% valuesgeo.NeuronsSpines = structfun(@(x) nangeomean(x,1),valneurons,'UniformOutput',false);
% geomeanvalues.NeuronsSpines = nangeomean(cell2mat(struct2cell(valuesgeo.NeuronsSpines)),1);
% count.NeuronsSpines = sum(~isnan(cell2mat(struct2cell(values.NeuronsSpines))),1);
% stdev.NeuronsSpines = std(cell2mat(struct2cell(values.NeuronsSpines)),0,1,'omitnan');
% se.NeuronsSpines = stdev.NeuronsSpines ./ sqrt(count.NeuronsSpines);
% 
% values.NeuronsSegmentsSpines = structfun(@(x) nanmean(x,1),valneuronsSegmentsmeans,'UniformOutput',false);
% meanvalues.NeuronsSegmentsSpines = nanmean(cell2mat(struct2cell(values.NeuronsSegmentsSpines)),1);
% valuesgeo.NeuronsSegmentsSpines = structfun(@(x) nangeomean(x,1),valneuronsSegmentsgeomeans,'UniformOutput',false);
% geomeanvalues.NeuronsSegmentsSpines = nangeomean(cell2mat(struct2cell(valuesgeo.NeuronsSegmentsSpines)),1);
% count.NeuronsSegmentsSpines = sum(~isnan(cell2mat(struct2cell(values.NeuronsSegmentsSpines))),1);
% stdev.NeuronsSegmentsSpines = std(cell2mat(struct2cell(values.NeuronsSegmentsSpines)),0,1,'omitnan');
% se.NeuronsSegmentsSpines = stdev.NeuronsSegmentsSpines ./ sqrt(count.NeuronsSegmentsSpines);
% 
% 
% NeuronNames=fieldnames(valneurons);
% for i=1:length(NeuronNames)
% currentmousename = NeuronNames{i};
% currentmousename = currentmousename(1:end-2);
%     if isfield (valmiceNeuronsmeans,currentmousename) == 0
%         valmiceNeuronsmeans.(currentmousename) = [];
%         valmiceNeuronsmeans.(currentmousename) = [valmiceNeuronsmeans.(currentmousename) ; values.NeuronsSpines.(NeuronNames{i})];
%         valmiceNeuronsgeomeans.(currentmousename) = [];
%         valmiceNeuronsgeomeans.(currentmousename) = [valmiceNeuronsgeomeans.(currentmousename) ; valuesgeo.NeuronsSpines.(NeuronNames{i})];
%         valmiceNeuronsSegmentsmeans.(currentmousename) = [];
%         valmiceNeuronsSegmentsmeans.(currentmousename) = [valmiceNeuronsSegmentsmeans.(currentmousename) ; values.NeuronsSegmentsSpines.(NeuronNames{i})];
%         valmiceNeuronsSegmentsgeomeans.(currentmousename) = [];
%         valmiceNeuronsSegmentsgeomeans.(currentmousename) = [valmiceNeuronsSegmentsgeomeans.(currentmousename) ; valuesgeo.NeuronsSegmentsSpines.(NeuronNames{i})];
%     else
%         valmiceNeuronsmeans.(currentmousename) = [valmiceNeuronsmeans.(currentmousename) ; values.NeuronsSpines.(NeuronNames{i})];
%         valmiceNeuronsgeomeans.(currentmousename) = [valmiceNeuronsgeomeans.(currentmousename) ; valuesgeo.NeuronsSpines.(NeuronNames{i})];
%         valmiceNeuronsSegmentsmeans.(currentmousename) = [valmiceNeuronsSegmentsmeans.(currentmousename) ; values.NeuronsSegmentsSpines.(NeuronNames{i})];
%         valmiceNeuronsSegmentsgeomeans.(currentmousename) = [valmiceNeuronsSegmentsgeomeans.(currentmousename) ; valuesgeo.NeuronsSegmentsSpines.(NeuronNames{i})];
%     end
% end


% Values at mouse level
values.MiceSpines = structfun(@(x) nanmean(x,1),valmice,'UniformOutput',false);
meanvalues.MiceSpines = nanmean(cell2mat(struct2cell(values.MiceSpines)),1);
valuesgeo.MiceSpines = structfun(@(x) nangeomean(x,1),valmice,'UniformOutput',false);
geomeanvalues.MiceSpines = nangeomean(cell2mat(struct2cell(valuesgeo.MiceSpines)),1);
count.MiceSpines = sum(~isnan(cell2mat(struct2cell(values.MiceSpines))),1);
stdev.MiceSpines = std(cell2mat(struct2cell(values.MiceSpines)),0,1,'omitnan');
se.MiceSpines = stdev.MiceSpines ./ sqrt(count.MiceSpines);

values.MiceSegmentsSpines = structfun(@(x) nanmean(x,1),valmiceSegmentsmeans,'UniformOutput',false);
meanvalues.MiceSegmentsSpines = nanmean(cell2mat(struct2cell(values.MiceSegmentsSpines)),1);
valuesgeo.MiceSegmentsSpines = structfun(@(x) nangeomean(x,1),valmiceSegmentsgeomeans,'UniformOutput',false);
geomeanvalues.MiceSegmentsSpines = nangeomean(cell2mat(struct2cell(valuesgeo.MiceSegmentsSpines)),1);
count.MiceSegmentsSpines = sum(~isnan(cell2mat(struct2cell(values.MiceSegmentsSpines))),1);
stdev.MiceSegmentsSpines = std(cell2mat(struct2cell(values.MiceSegmentsSpines)),0,1,'omitnan');
se.MiceSegmentsSpines = stdev.MiceSegmentsSpines ./ sqrt(count.MiceSegmentsSpines);

values.MiceNeuronsSpines = structfun(@(x) nanmean(x,1),valmiceNeuronsmeans,'UniformOutput',false);
meanvalues.MiceNeuronsSpines = nanmean(cell2mat(struct2cell(values.MiceNeuronsSpines)),1);
valuesgeo.MiceNeuronsSpines = structfun(@(x) nangeomean(x,1),valmiceNeuronsgeomeans,'UniformOutput',false);
geomeanvalues.MiceNeuronsSpines = nangeomean(cell2mat(struct2cell(valuesgeo.MiceNeuronsSpines)),1);
count.MiceNeuronsSpines = sum(~isnan(cell2mat(struct2cell(values.MiceNeuronsSpines))),1);
stdev.MiceNeuronsSpines = std(cell2mat(struct2cell(values.MiceNeuronsSpines)),0,1,'omitnan');
se.MiceNeuronsSpines = stdev.MiceNeuronsSpines ./ sqrt(count.MiceNeuronsSpines);


values.MiceNeuronsSegmentsSpines = structfun(@(x) nanmean(x,1),valmiceNeuronsSegmentsmeans,'UniformOutput',false);
meanvalues.MiceNeuronsSegmentsSpines = nanmean(cell2mat(struct2cell(values.MiceNeuronsSegmentsSpines)),1);
valuesgeo.MiceNeuronsSegmentsSpines = structfun(@(x) nangeomean(x,1),valmiceNeuronsSegmentsmeans,'UniformOutput',false);
geomeanvalues.MiceNeuronsSpines = nangeomean(cell2mat(struct2cell(valuesgeo.MiceNeuronsSegmentsSpines)),1);
count.MiceNeuronsSegmentsSpines = sum(~isnan(cell2mat(struct2cell(values.MiceNeuronsSegmentsSpines))),1);
stdev.MiceNeuronsSegmentsSpines = std(cell2mat(struct2cell(values.MiceNeuronsSegmentsSpines)),0,1,'omitnan');
se.MiceNeuronsSegmentsSpines = stdev.MiceNeuronsSegmentsSpines ./ sqrt(count.MiceNeuronsSegmentsSpines);



values.segmentnames = val.segment;

retPool.values = values;
retPool.meanvalues = meanvalues;
retPool.valuesgeo = valuesgeo;
retPool.geomeanvalues = geomeanvalues;
retPool.count = count;
retPool.stdev = stdev;
retPool.se = se;
retPool.stat = StatName;
retPool.condition = Condition;
retPool.SessionConditions = SessionConditions;
retPool.additions = additions;
retPool.eliminations = eliminations;
retPool.survivalrate = survivalrate;
retPool.survivalD1rate = survivalD1rate;
retPool.countall = countall;
retPool.dist = dist; %added by ELO on 2/27/21
%clearvars currentvalues i MapName myCondSegmentTable myGoodSegmentTable PersCond rows SessionConditions SegmentNumber meanvalues values
end