clc
close all
clear all;
%% add codes to path
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';
animalprefix = 'PPP7';
day = 12;

daystring = num2str(day);
dir = [dir,animalprefix,'/day',daystring,'/'];
prefix = ['day',daystring];

eps_RUN = [3,5]; % Tmaze epochs
eps_postSleep = 6; % POST sleep
eps_preSleep = 2; % PRE sleep
example_assemble_num = [6,2]; % example ensembles shown in Fig. S4b-d
speedthresh = 2; % cm/s, running velocity threshold
cell_excluded = []; % cell excluded
%% load spikes
spikeStructure = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"dCA1");

spikes = [];
cell_count = 0;
for i = 1:length(spikeStructure.times)
    if ~ismember(i,cell_excluded)
        cell_count = cell_count + 1;
        timestamps = spikeStructure.times{i};
        id = ones(size(timestamps))*i; % the id for that unit
        spikes = [spikes; timestamps id]; % add [timestamps id] for this unit to the matrix
    end
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
load(fullfile(dir,[prefix,'.ripples.events.mat'])) % ripple
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
    %% normalized by the presleep
    triggerid = find(x >=0 & x <= 0.1);
    preStrength = nanmean(mPre(:,triggerid),2);
    for i=1:nComponents
        strength(:,i+1) = (strength(:,i+1)-preStrength(i))/preStrength(i);
    end
    %%
    example_assemble_weights{ep} = weights(:,example_assemble_num((ep-1)/2));
    
    strength_ep = [strength(:,1),nanmean(strength(:,1+example_assemble_num((ep-1)/2)),2)];
    strength_all{ep} = strength_ep;
end
%% get RS for the example assemblies
strength_remote = strength_all{eps_RUN(1)};
strength_remote(:,2) = (strength_remote(:,2)  - nanmean(strength_remote(:,2)))./nanstd(strength_remote(:,2));

strength_recent = strength_all{eps_RUN(2)};
strength_recent(:,2) = (strength_recent(:,2)  - nanmean(strength_recent(:,2)))./nanstd(strength_recent(:,2));

load(fullfile(dir,[prefix,'.MergePoints.events.mat'])) % Session info
epochtimes = MergePoints.timestamps(eps_postSleep,:);

valididx = find(strength_remote(:,1) >= epochtimes(1) & strength_remote(:,1) <= epochtimes(2));
strength_remote = strength_remote(valididx,:);
% smooth out jumping point
RS_remote = strength_remote(:,2);
RS_time = strength_remote(:,1);
RS_interp = interp1(RS_time(~isnan(RS_remote)),RS_remote(~isnan(RS_remote)),RS_time);
RS_srate = 1/nanmedian(diff(RS_time));
RS_remote_fit = locsmooth(RS_interp,RS_srate,5);

valididx = find(strength_recent(:,1) >= epochtimes(1) & strength_recent(:,1) <= epochtimes(2));
strength_recent = strength_recent(valididx,:);
% smooth out jumping point
RS_recent = strength_recent(:,2);
RS_time = strength_recent(:,1);
RS_interp = interp1(RS_time(~isnan(RS_recent)),RS_recent(~isnan(RS_recent)),RS_time);
RS_recent_fit = locsmooth(RS_interp,RS_srate,5);
%% sleep state
load(fullfile(dir,[prefix,'.SleepState.states.mat'])) % sleep states
Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM
sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
Sleep_state_ep = Sleep_state(sleepss_id,:);
SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
SWSid = find(Sleep_state_ep(:,2) == 3);
SWS_vec_ep(SWSid) = 1;
SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));
SWSdur = SWSlist(:,2) - SWSlist(:,1);
SWSlist = SWSlist(find(SWSdur > 10),:);

load(fullfile(dir,[prefix,'.animal.behavior.mat'])) % behavior
pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % time, normalized diameter
%% restrict RS to NREM
pupil_trace_SWS = [];
RS_recent_SWS = [];
RS_remote_SWS = [];
for SWSseg = 1:length(SWSlist(:,1))
    current_SWS = SWSlist(SWSseg,:);
    tempid  = find(pupil_trace(:,1) >= current_SWS(1) & pupil_trace(:,1) <= current_SWS(2));
    pupil_trace_SWS = [pupil_trace_SWS;pupil_trace(tempid,:)];

    tempid  = find(RS_time >= current_SWS(1) & RS_time <= current_SWS(2));
    RS_recent_SWS =[RS_recent_SWS;RS_time(tempid),RS_recent_fit(tempid)'];
    RS_remote_SWS =[RS_remote_SWS;RS_time(tempid),RS_remote_fit(tempid)'];
end
%% zscore RS
RS_recent_SWS_z = RS_recent_SWS;
RS_remote_SWS_z = RS_remote_SWS;
RS_recent_SWS_z(:,2) = (RS_recent_SWS(:,2) - nanmean(RS_recent_SWS(:,2)))./nanstd(RS_recent_SWS(:,2));
RS_remote_SWS_z(:,2) = (RS_remote_SWS(:,2) - nanmean(RS_remote_SWS(:,2)))./nanstd(RS_remote_SWS(:,2));
%% example RS segment
timerange = [23695, 23892];
tempid  = find(pupil_trace_SWS(:,1) >= timerange(1) & pupil_trace_SWS(:,1) <= timerange(2));
pupil_trace_example = pupil_trace_SWS(tempid,:);
pupil_d_interp = interp1(pupil_trace_example(~isnan(pupil_trace_example(:,2)),1),pupil_trace_example(~isnan(pupil_trace_example(:,2)),2),pupil_trace_example(:,1));
pupil_trace_example(:,2) = locsmooth(pupil_d_interp,30,2);


tempid  = find(RS_recent_SWS(:,1) >= timerange(1) & RS_recent_SWS(:,1) <= timerange(2));
RS_recent_example = RS_recent_SWS(tempid,:);
RS_remote_example = RS_remote_SWS(tempid,:);

RS_recent_example(:,2) = (RS_recent_example(:,2) - nanmean(RS_recent_example(:,2)))./nanstd(RS_recent_example(:,2));
RS_remote_example(:,2) = (RS_remote_example(:,2) - nanmean(RS_remote_example(:,2)))./nanstd(RS_remote_example(:,2));

figure,
plot(pupil_trace_example(:,1),pupil_trace_example(:,2));
hold on

plot(RS_recent_example(:,1),RS_recent_example(:,2));
plot(RS_remote_example(:,1),RS_remote_example(:,2));
%% sort ensemble weights
recent_weight = example_assemble_weights{5};
remote_weight = example_assemble_weights{3};

assemble_mean = nanmean(abs(recent_weight));
assemble_std = nanstd(abs(recent_weight));
threshold_recent = assemble_mean + 3*assemble_std; % van de Ven 2016
recent_memberid = find(recent_weight > threshold_recent);

assemble_mean = nanmean(abs(remote_weight));
assemble_std = nanstd(abs(remote_weight));
threshold_remote = assemble_mean + 3*assemble_std; % van de Ven 2016
remote_memberid = find(remote_weight > threshold_remote);

nonmemberids = setdiff(1:length(recent_weight),[remote_memberid;recent_memberid]');
recent_weight_sorted = recent_weight([recent_memberid;remote_memberid;nonmemberids']);
remote_weight_sorted = remote_weight([recent_memberid;remote_memberid;nonmemberids']);
%% plot place fields of the members across environments
load(fullfile(dir,  [prefix,'.firingMapsAvg_2D.cellinfo.mat']));% load firing rate maps
figure,
subplot(241)
imagesc(firingMaps_2D{3}.rateMaps{recent_memberid(1)})
axis square
mycolormap = inferno(100);
mycolormap = [1,1,1;mycolormap];
colormap(mycolormap)
caxis([-1,8])
subplot(242)
imagesc(firingMaps_2D{3}.rateMaps{recent_memberid(2)})
axis square
mycolormap = inferno(100);
mycolormap = [1,1,1;mycolormap];
colormap(mycolormap)
caxis([-1,8])
subplot(243)
imagesc(firingMaps_2D{3}.rateMaps{remote_memberid(1)})
axis square
mycolormap = inferno(100);
mycolormap = [1,1,1;mycolormap];
colormap(mycolormap)
caxis([-1,12])
subplot(244)
imagesc(firingMaps_2D{3}.rateMaps{remote_memberid(2)})
axis square
mycolormap = inferno(100);
mycolormap = [1,1,1;mycolormap];
colormap(mycolormap)
caxis([-1,12])

subplot(245)
imagesc(firingMaps_2D{5}.rateMaps{recent_memberid(1)})
axis square
mycolormap = inferno(100);
mycolormap = [1,1,1;mycolormap];
colormap(mycolormap)
caxis([-1,12])
subplot(246)
imagesc(firingMaps_2D{5}.rateMaps{recent_memberid(2)})
axis square
mycolormap = inferno(100);
mycolormap = [1,1,1;mycolormap];
colormap(mycolormap)
caxis([-1,5])
subplot(247)
imagesc(firingMaps_2D{5}.rateMaps{remote_memberid(1)})
axis square
mycolormap = inferno(100);
mycolormap = [1,1,1;mycolormap];
colormap(mycolormap)
caxis([-1,8])
subplot(248)
imagesc(firingMaps_2D{5}.rateMaps{remote_memberid(2)})
axis square
mycolormap = inferno(100);
mycolormap = [1,1,1;mycolormap];
colormap(mycolormap)
caxis([-1,8])


