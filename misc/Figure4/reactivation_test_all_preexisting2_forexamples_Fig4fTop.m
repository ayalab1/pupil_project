clc
close all
clear all;
%% add code to path
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';
animalprefix = 'PPP4';
day = 11;
daystring = num2str(day);
dir = [dir,animalprefix,'/day',daystring,'/'];
prefix = ['day',daystring];
eps_RUN = [2,4];
eps_postSleep = [3,5]; 
eps_preSleep = 1;
PRE_component = 2; % example assembly number
%% define parameter
speedthresh = 2; % cm/s, running speed threshold
%% load spikes
if strcmp(animalprefix,'PPP4') || (strcmp(animalprefix,'PPP8') && day == 8)
    spikeStructure = importSpikes('basepath',dir,'CellType',"Pyramidal Cell");
elseif strcmp(animalprefix,'PPP7')
    if day == 18 || day == 17
        spikeStructure = importSpikes('basepath',dir,'CellType',"Pyramidal Cell");
    else
        spikeStructure = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"dCA1");
    end
else
    spikeStructure = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"CA1");
end
spikes = [];
for i = 1:length(spikeStructure.times)
    timestamps = spikeStructure.times{i};
    id = ones(size(timestamps))*i; % the id for that unit
    spikes = [spikes; timestamps id]; % add [timestamps id] for this unit to the matrix
end
spikes = sortrows(spikes); % sort spikes accoring to their time    
%% load data
load(fullfile(dir,[prefix,'.animal.behavior.mat'])) % behavior
vel  = behavior.speed_smooth; % running speed

load(fullfile(dir,[prefix,'.MergePoints.events.mat'])) % Session info
postSleep_times = MergePoints.timestamps(eps_postSleep,:);
preSleep_times = MergePoints.timestamps(eps_preSleep,:);

% Load NREM sleep epochs:
load(fullfile(dir,[prefix,'.SleepState.states.mat'])) 
NREM = bsxfun(@plus,SleepState.ints.NREMstate,[-1 0]);
postSleep = postSleep_times;
preSleep = preSleep_times;

RUN = vec2list(vel(:,2) > speedthresh,vel(:,1)); % generate [start end] list of RUN epochs
taskIntervals = SplitIntervals(RUN,'pieceSize',0.05); % 50 ms bins in behavior

if strcmp(prefix,'day28') && strcmp(animalprefix,'PPP4')
    load(fullfile(dir,[prefix,'.ripples_task.events.mat'])) % ripple
elseif strcmp(prefix,'day8') && strcmp(animalprefix,'PPP7')
    load(fullfile(dir,[prefix,'.dorsalripples.events.mat'])) % ripple
else
    load(fullfile(dir,[prefix,'.ripples.events.mat'])) % ripple
end
%% calculate activity templates from RUN epochs
templates = ActivityTemplates(spikes,'bins',taskIntervals,'mode','ica');
binSize = 0.05; step = 0.01;
bins = Bins(0,MergePoints.timestamps(end),binSize,step); 
% We don't necessarity need the bins outside of pre-task sleep and post-task sleep
bins = Restrict(bins,[preSleep;postSleep]); % this saves memory when computing the strength
strength = ReactivationStrength(spikes,templates,'bins',bins);
%% calculate ripple-triggered reactivation strength
preRipples = Restrict(ripples.timestamps(:,1),preSleep);
postRipples = Restrict(ripples.timestamps(:,1),postSleep);

nComponents = size(templates,3);
nBins = 101; durations = [-1 1]*0.5;
[mPre,mPost] = deal(nan(nComponents,nBins)); % pre-allocate matrices
x = linspace(durations(1),durations(2),nBins);
baseline_ids = find(x < -0.1);
ripple_ids = find(x > 0 & x < 0.1);

% test for significant components
for i=1:nComponents
    % use "PETH" to compute the peri-event time histogram around ripples
    % for significant test
    temp_pre = PETH(strength(:,[1 1+i]),preRipples,'durations',durations,'nBins',nBins);
    basline = nanmean(temp_pre(:,baseline_ids),2);
    rips = nanmean(temp_pre(:,ripple_ids),2);
    pre_pval = signrank(basline,rips);% test significance

    temp_post = PETH(strength(:,[1 1+i]),postRipples,'durations',durations,'nBins',nBins);
    basline = nanmean(temp_post(:,baseline_ids),2);
    rips = nanmean(temp_post(:,ripple_ids),2);
    post_pval = signrank(basline,rips);% test significance

    if  post_pval < 0.05 || pre_pval < 0.05 % if ripple-triggered RS passes significance threshold
        Components_sign(i) = 1;
    else
        Components_sign(i) = 0;
    end
end

% calculate RS difference from pre to post
for i=1:nComponents
    % use "mPETH" to compute the mean peri-event time histogram around ripples
    mPre(i,:) = mPETH(strength(:,[1 1+i]),preRipples,'durations',durations,'nBins',nBins);
    mPost(i,:) = mPETH(strength(:,[1 1+i]),postRipples,'durations',durations,'nBins',nBins);
end
Post_Pre = sign(nanmean(mPost(:,ripple_ids),2) - nanmean(mPre(:,ripple_ids),2));
%% defined pre-existing and plastic cell assemblies
preexisting_Components = find(Components_sign & Post_Pre' <= 0);
plastic_Components = find(Components_sign & Post_Pre' > 0);
%% Plot the result of all preexisting and plastic assemblies from pre to post:
figure(3); clf;
set(gcf,'position',[300 200 1500 700]);
matrices = {mPre(preexisting_Components,:), mPost(preexisting_Components,:),mPost(preexisting_Components,:)-mPre(preexisting_Components,:)};
colors = {[0 0 0],[0.5 0 0],[1 0 0]};
titles = {'pre-task sleep','post-task sleep','difference'};
for i=1:3
    subplot(2,4,i); PlotColorMap(matrices{i},'x',x);
    set(gca,'box','off','TickDir','out','fontsize',12);
    ylabel('component ID');
    PlotHVLines(0,'v','w--','linewidth',2);
    title(titles{i});
end
for i=1:3
    subplot(2,4,i+4); 
    semplot(x,matrices{i},colors{i});
    xlabel('time from ripple start (s)');
    ylabel('mean activation strength');
end
clims % set the same color limits on the three panels
EquateScales(4,5,6); % set the same y-axis limits on the bottom three panels

% add a line at 0
for i=1:3
    subplot(2,4,i+4);
    PlotHVLines(0,'v','k--','linewidth',2);
    PlotIntervals([0 0.1]);
    set(gca,'box','off','TickDir','out','fontsize',12);
end

figure(4); clf;
set(gcf,'position',[300 200 1500 700]);
matrices = {mPre(plastic_Components,:), mPost(plastic_Components,:),mPost(plastic_Components,:)-mPre(plastic_Components,:)};
colors = {[0 0 0],[0.5 0 0],[1 0 0]};
titles = {'pre-task sleep','post-task sleep','difference'};
for i=1:3
    subplot(2,4,i); PlotColorMap(matrices{i},'x',x);
    set(gca,'box','off','TickDir','out','fontsize',12);
    ylabel('component ID');
    PlotHVLines(0,'v','w--','linewidth',2);
    title(titles{i});
end
for i=1:3
    subplot(2,4,i+4); 
    semplot(x,matrices{i},colors{i});
    xlabel('time from ripple start (s)');
    ylabel('mean activation strength');
end
clims % set the same color limits on the three panels
EquateScales(4,5,6); % set the same y-axis limits on the bottom three panels

% add a line at 0
for i=1:3
    subplot(2,4,i+4);
    PlotHVLines(0,'v','k--','linewidth',2);
    PlotIntervals([0 0.1]);
    set(gca,'box','off','TickDir','out','fontsize',12);
end
%% separate RS into two groups
strength_preexisting = [strength(:,1),strength(:,preexisting_Components+1)];
strength_plastic = [strength(:,1),strength(:,plastic_Components+1)];
%% get RS of the example assembly by pupil-tile ripple
PRE_component = preexisting_Components(PRE_component);

load([dir,prefix,'.PupilRipples','EP',num2str(3),'.events.mat']) % ripple info by pupil tile; use epoch3

for tile = 1:3 % small, middle, large
    validid  = find(ripples.pupil_tile == (tile*2-1) | ripples.pupil_tile == (tile*2));
    postRipples = ripples.timestamps(validid,1);
    temp_post_PRE = PETH(strength(:,[1 PRE_component+1]),postRipples,'durations',durations,'nBins',nBins);
   
    % zscore by baseline distribution
    baseline = temp_post_PRE(:,baseline_ids);
    baseline_mean = nanmean(baseline(:));  
    baseline_std = nanstd(baseline(:));  
    temp_post_PRE = (temp_post_PRE - baseline_mean)./baseline_std;
   
    PRE_tile{tile} = temp_post_PRE;
    clear temp_post_PRE
end
%% get the RS of the example assembly during PRE and POST sessions
load(fullfile(dir,[prefix,'.ripples.events.mat'])) % ripple across all sessions
preRipples = Restrict(ripples.timestamps(:,1),preSleep);
postRipples = Restrict(ripples.timestamps(:,1),postSleep);

RS_post_PRE = PETH(strength(:,[1 PRE_component+1]),postRipples,'durations',durations,'nBins',nBins);
% zscore by baseline distribution
baseline = RS_post_PRE(:,baseline_ids);
baseline_mean = nanmean(baseline(:));
baseline_std = nanstd(baseline(:));
RS_post_PRE = (RS_post_PRE - baseline_mean)./baseline_std;

RS_pre_PRE = PETH(strength(:,[1 PRE_component+1]),preRipples,'durations',durations,'nBins',nBins);
% zscore by baseline distribution
baseline = RS_pre_PRE(:,baseline_ids);
baseline_mean = nanmean(baseline(:));
baseline_std = nanstd(baseline(:));
RS_pre_PRE = (RS_pre_PRE - baseline_mean)./baseline_std;



