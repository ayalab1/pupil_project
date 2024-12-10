% this is the main script to calculate reactivation strength based on
% pipelineReactivation.m in AYALab neurocode
% Detecting task-related components with PCA/ICA and observing the
%     of those components in sleep epochs
clc
close all
clear all;
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';
animalprefix = 'PPP8';
day = 8;
daystring = num2str(day);
dir = [dir,animalprefix,'/day',daystring,'/'];
prefix = ['day',daystring];
eps_RUN = [2,4]; % behavior epochs
eps_postSleep = [3,5];  % POST sleep epochs
eps_preSleep = [1]; % PRE sleep epochs
%% define parameters
speedthresh = 2; % cm/s, running speed threshold
savedata = 1; % save results?
%%  load spikes
% importSpikes function from AYALab neurocode
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
%% load task info
load(fullfile(dir,[prefix,'.animal.behavior.mat'])) % behavior
vel  = behavior.speed_smooth; % running speed

load(fullfile(dir,[prefix,'.MergePoints.events.mat'])) % Session info
postSleep_times = MergePoints.timestamps(eps_postSleep,:);
preSleep_times = MergePoints.timestamps(eps_preSleep,:);

% Load NREM sleep epochs:
load(fullfile(dir,[prefix,'.SleepState.states.mat'])) 
NREM = bsxfun(@plus,SleepState.ints.NREMstate,[-1 0]);
postSleep = Restrict(NREM,postSleep_times);
preSleep = Restrict(NREM,preSleep_times);

RUN = vec2list(vel(:,2) > speedthresh,vel(:,1)); % generate [start end] list of immobile epochs
taskIntervals = SplitIntervals(RUN,'pieceSize',0.05); % 50 ms bins in behavior
%% load ripples
if strcmp(prefix,'day28') && strcmp(animalprefix,'PPP4')
    load(fullfile(dir,[prefix,'.ripples_task.events.mat'])) % ripple
elseif strcmp(prefix,'day8') && strcmp(animalprefix,'PPP7')
    load(fullfile(dir,[prefix,'.dorsalripples.events.mat'])) % ripple
else
    load(fullfile(dir,[prefix,'.ripples.events.mat'])) % ripple
end
%% calculate templates
templates = ActivityTemplates(spikes,'bins',taskIntervals,'mode','ica');
binSize = 0.05; step = 0.01;
bins = Bins(0,MergePoints.timestamps(end),binSize,step); 
% We don't necessarity need the bins outside of pre-task sleep and post-task sleep
bins = Restrict(bins,[preSleep;postSleep]); % this saves memory when computing the strength
strength = ReactivationStrength(spikes,templates,'bins',bins);

preRipples = Restrict(ripples.timestamps(:,1),preSleep);
postRipples = Restrict(ripples.timestamps(:,1),postSleep);

nComponents = size(templates,3);
nBins = 101; durations = [-1 1]*0.5;
[mPre,mPost] = deal(nan(nComponents,nBins)); % pre-allocate matrices
for i=1:nComponents
    % use "mPETH" to compute the mean peri-event time histogram around ripples
    mPre(i,:) = mPETH(strength(:,[1 1+i]),preRipples,'durations',durations,'nBins',nBins);
    mPost(i,:) = mPETH(strength(:,[1 1+i]),postRipples,'durations',durations,'nBins',nBins);
end
x = linspace(durations(1),durations(2),nBins);
%% Plot the result (plotting functions from AYALab neurocode) 
figure(3); clf;
set(gcf,'position',[300 200 1500 700]);
matrices = {mPre,mPost,mPost-mPre};
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
%% save data
strength = [strength(:,1),nanmean(strength(:,2:end),2)];
if savedata
    save(sprintf('%s/%sReactivationStrength.mat', dir, prefix), 'strength');
end





