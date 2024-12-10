clc
clear all
close all
%% add codes to path
% the chronux toolbox is needed
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
% directory saved all the replay decoding results
replaydir = '/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/wb_scripts/Hongyu_scripts/Replay_decoding_20230526/';
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';
animalprefix = 'PPP7';
day = 8;
daystring = num2str(day);
dir = [dir,animalprefix,'/day',daystring,'/'];
prefix = ['day',daystring];
animalname_prefix = [animalprefix,'-',prefix];
ep = 5; % epoch
%% load file
load(fullfile(dir,[prefix,'.MergePoints.events.mat'])) % Session info
epochtimes = MergePoints.timestamps(ep,:); % epoch start and end time
load(fullfile(dir,[prefix,'.animal.behavior.mat'])) % behavior
pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % time, normalized diameter

% load raw ripples
if strcmp(prefix,'day28') && strcmp(animalprefix,'PPP4')
    load(fullfile(dir,[prefix,'.ripples_task.events.mat'])) % ripple
elseif strcmp(prefix,'day8') && strcmp(animalprefix,'PPP7')
    load(fullfile(dir,[prefix,'.dorsalripples.events.mat'])) % ripple
else
    load(fullfile(dir,[prefix,'.ripples.events.mat'])) % ripple
end
riplist_all = [ripples.timestamps(:,1),ripples.timestamps(:,1) + ripples.duration]; % [start, end]

% load ripples, defined by HSEs
load(sprintf('%s/%srippleHSEs.events.mat', dir, prefix));% rippleHSE info
riptimes = [rippleHSEs{day}{ep}.starttime',rippleHSEs{day}{ep}.endtime'];

dur = 1000*(riptimes(:,2) - riptimes(:,1));
keepidx = find(dur >= 50);%at least 5 bins, 50 ms for 10ms bins
HSElist = riptimes(keepidx,:);

% load replay information
load(sprintf('%s%sreplayHSEseqdecode_CA1_%02d.mat', replaydir,animalname_prefix,ep));
eventidx = replaytrajactory{ep}.eventidx;
candidate_times = HSElist(eventidx,:);
replay_pval = 1 - min(replaytrajactory{ep}.pvalue_r,[],2); % 1 - p_shuf

HSElist_all = [candidate_times,replay_pval];
for event = 1:length(riplist_all(:,1))
    temp = find(HSElist_all(:,1) <= riplist_all(event,1) & HSElist_all(:,2) >= riplist_all(event,1));
    if ~isempty(temp)
        riplist_all(event,3) = HSElist_all(temp,3);
    else
        riplist_all(event,3) = 0; %pval = 0 for non-candidate events
    end
end

% sleep states
load(fullfile(dir,[prefix,'.SleepState.states.mat'])) % sleep states
Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM
sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
Sleep_state_ep = Sleep_state(sleepss_id,:);

% NREM
SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
SWSid = find(Sleep_state_ep(:,2) == 3);
SWS_vec_ep(SWSid) = 1;
SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));
SWSdur = SWSlist(:,2) - SWSlist(:,1);
SWSlist = SWSlist(find(SWSdur > 30),:);
%% Ripple pval trace
ripple_id = find(riplist_all(:,1) <= epochtimes(2) & riplist_all(:,1) >= epochtimes(1));
ripples_ep = riplist_all(ripple_id,:);
ripples_times = [];
for rip = 1:length(ripples_ep(:,1))
    replay_prob = ripples_ep(rip,3);
    ripples_times = [ripples_times; ripples_ep(rip,1)* ones(round(replay_prob*100),1)];
end

ripple_events(1).times = ripples_times;
[R,t,~] = psth(ripple_events,5); % psth function from chronux toolbox
%% pupil diameter trace
tempid  = find(pupil_trace(:,1) >= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
pupil_d = pupil_trace(tempid,2);
pupil_time = pupil_trace(tempid,1);
pupil_SR = 1/nanmedian(diff(pupil_trace(:,1)));
pupil_d_interp = interp1(pupil_time(~isnan(pupil_d)),pupil_d(~isnan(pupil_d)),pupil_time);
pupil_d_fit = locsmooth(pupil_d_interp,pupil_SR,5);

pupil_d_fit = (pupil_d_fit-nanmean(pupil_d_fit))/nanstd(pupil_d_fit);
pupil_d_fit(isnan(pupil_d_fit)) = 0;

rippleprob_interp = interp1(t,R,pupil_time);
rippleprob_interp(isnan(rippleprob_interp)) = 0;
%% plot the trace for the whole session
figure,
subplot(211)
plot(Sleep_state_ep(:,1),Sleep_state_ep(:,2))
subplot(212)
hold on
plot(pupil_time,(-1*((pupil_d_fit-0.4)*10))/1.1,'linewidth',1)
plot(pupil_time+5,rippleprob_interp/2,'k','linewidth',2)
%% restrict to the example period, find the replay event IDs
valiid = find(ripples_ep(:,2) <= 23000 & ripples_ep(:,1) >= 22300);
ripple_valid = ripples_ep(validid,:);
replay_IDs = find(ripple_valid(:,3) > 0.975);
replay_valid = ripple_valid(replay_IDs,:);

replaynum = length(replay_IDs);
replay_list = [];
for r = 1:replaynum
    replaytime = replay_valid(r,1);
    [~,i] = min(abs(pupil_time - replaytime));
    replay_list = [replay_list;pupil_time(i),rippleprob_interp(i)];
end
%% find the non-replay event ID
nonreplay_IDs = find(ripple_valid(:,3) < 0.8);
nonreplay_valid = ripple_valid(nonreplay_IDs,:);

nonreplaynum = length(nonreplay_IDs);
nonreplay_list = [];
for r = 1:nonreplaynum
    nonreplaytime = nonreplay_valid(r,1);
    [~,i] = min(abs(pupil_time - nonreplaytime));
    nonreplay_list = [nonreplay_list;pupil_time(i),rippleprob_interp(i)];
end

% %%
% ripple_pupil_corr_all = [];
% for seg = 1:length(SWSlist(:,1))
%     currentseg = SWSlist(seg,:);
%     segid = find(pupil_time >= currentseg(1) & pupil_time <= currentseg(2));
%     if ~isempty(segid)
%         validid = find(~isnan(pupil_d_fit(segid)));
%         pupil_d_fit_seg = pupil_d_fit(segid(validid));
%         rippleprob_seg = rippleprob_interp(segid(validid));
% 
%         % flip the pupil
%         [ripple_pupil_corr,lags] = xcorr(-1*pupil_d_fit_seg,rippleprob_seg,round(30*pupil_SR),'biased');
%         ripple_pupil_corr_all = [ripple_pupil_corr_all;ripple_pupil_corr];
%     end
% end
% %%
% lags = lags/pupil_SR;
% ripple_trigger_pupil_stats = [lags',nanmean(ripple_pupil_corr_all)',nanstd(ripple_pupil_corr_all)',sum(double((~isnan(ripple_pupil_corr_all))))'];
% figure,plot(ripple_trigger_pupil_stats(:,1),ripple_trigger_pupil_stats(:,2))
%     
    

    
    
    
    
