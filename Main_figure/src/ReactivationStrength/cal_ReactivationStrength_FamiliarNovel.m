% this is the main script to calculate reactivation strength in familiar 
% vs. novel environment, normalized by PRE sleep
% based on pipelineReactivation.m in AYALab neurocode
% Detecting task-related components with PCA/ICA and observing the
%     of those components in sleep epochs
clc
close all
clear all;
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';
animalprefix = 'PPP7';
day = 14;
if  strcmp(animalprefix,'PVR4')|| strcmp(animalprefix,'PVR5')
    daystring = num2str(day);
    prefix = ['pupil',daystring];
    dir = [dir,animalprefix,'/pupil',daystring,'/'];
else
    daystring = num2str(day);
    dir = [dir,animalprefix,'/day',daystring,'/'];
    prefix = ['day',daystring];
end
eps_RUN = [2,4]; % Tmaze epochs
eps_postSleep = 5; %POST sleep
eps_preSleep = 1; %PRE sleep
%% define parameters
speedthresh = 2; % cm/s, running velocity threshold
savedata = 1; % save data?
%% load spikes from CA1 PYRs
if strcmp(animalprefix,'PPP4') ||  strcmp(animalprefix,'PPP15') || strcmp(animalprefix,'PPP13') ||... 
        (strcmp(animalprefix,'PPP8') && day == 8) || (strcmp(animalprefix,'PPP8') && (day == 15 || day == 14 || day == 17 || day == 19))||...
        strcmp(animalprefix,'PVR4') ||  strcmp(animalprefix,'PVR5') 
    spikeStructure = importSpikes('basepath',dir,'CellType',"Pyramidal Cell");
elseif strcmp(animalprefix,'PPP7') % dorsal CA1 cells only
    if day == 18 || day == 17
        spikeStructure = importSpikes('basepath',dir,'CellType',"Pyramidal Cell");
    else
        spikeStructure = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"dCA1");
    end
else
    spikeStructure = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"CA1");
end

spikes = [];
cell_count = 0;
for i = 1:length(spikeStructure.times)
    cell_count = cell_count + 1;
    timestamps = spikeStructure.times{i};
    id = ones(size(timestamps))*i; % the id for that unit
    spikes = [spikes; timestamps id]; % add [timestamps id] for this unit to the matrix
end
spikes = sortrows(spikes); % sort spikes accoring to their time    
%% load session info and sleep states
load(fullfile(dir,[prefix,'.animal.behavior.mat'])) % behavior
vel  = behavior.speed_smooth; % running speed

load(fullfile(dir,[prefix,'.MergePoints.events.mat'])) % Session info
postSleep_times = MergePoints.timestamps(eps_postSleep,:);
preSleep_times = MergePoints.timestamps(eps_preSleep,:);
%% Load NREM sleep epochs
load(fullfile(dir,[prefix,'.SleepState.states.mat']))
NREM = bsxfun(@plus,SleepState.ints.NREMstate,[-1 0]);
postSleep = Restrict(NREM,postSleep_times);
preSleep = Restrict(NREM,preSleep_times);
%% load ripples
if strcmp(prefix,'day28') && strcmp(animalprefix,'PPP4')
    load(fullfile(dir,[prefix,'.ripples_task.events.mat'])) % ripple
elseif strcmp(prefix,'day8') && strcmp(animalprefix,'PPP7')
    load(fullfile(dir,[prefix,'.dorsalripples.events.mat'])) % ripple
else
    load(fullfile(dir,[prefix,'.ripples.events.mat'])) % ripple
end
%% calculate the template for recent and remote memories separately
for ep = eps_RUN
    epochtimes = MergePoints.timestamps(ep,:);
    vel_epID = find(vel(:,1) >= epochtimes(1) & vel(:,1) <= epochtimes(2));
    vel_ep = vel(vel_epID,:);
    
    % caculate template activity
    RUN = vec2list(vel_ep(:,2) > speedthresh,vel_ep(:,1)); % generate [start end] list of immobile epochs
    validid = find((RUN(:,2)-RUN(:,1)) > 0.05);
    RUN = RUN(validid,:);
    taskIntervals = SplitIntervals(RUN,'pieceSize',0.05); % 50 ms bins in behavior
    [templates,~,weights,~] = ActivityTemplates(spikes,'bins',taskIntervals,'mode','ica');
    binSize = 0.05; step = 0.01;
    bins = Bins(0,MergePoints.timestamps(end),binSize,step); 
    
    % caculate reactivated activity
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

    %% Plot the result
    x = linspace(durations(1),durations(2),nBins);
    figure(ep); clf;
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
    %% normalized by the PRE sleep
    triggerid = find(x >=0 & x <= 0.1);
    preStrength = nanmean(mPre(:,triggerid),2);
    for i=1:nComponents
        strength(:,i+1) = (strength(:,i+1)-preStrength(i))/preStrength(i);
    end
    
    strength_ep = [strength(:,1),nanmean(strength(:,2:end),2)];
    strength_all{ep} = strength_ep;
end
%% save data
if savedata
    save(sprintf('%s/%sReactivationStrength_diffTemp_PREnorm.mat', dir, prefix), 'strength_all');
end





