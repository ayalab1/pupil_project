clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
% directory saved all the replay decoding results
replaydir = '/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/wb_scripts/Hongyu_scripts/Replay_decoding_20230526/';
dir = '/Volumes/Extreme Pro/Pupil_data/';

% list of all data
animal_info = [{'PPP4'},{8},{3};...
    {'PPP4'},{10},{5};...
    {'PPP4'},{11},{1};...
    {'PPP4'},{11},{3};...
    {'PPP4'},{11},{5};...
    {'PPP7'},{8},{3};...
    {'PPP7'},{8},{5};...
    {'PPP7'},{12},{1};...
    {'PPP7'},{12},{4};...
    {'PPP7'},{14},{1};...
    {'PPP8'},{7},{1};...
    {'PPP8'},{7},{3};...
    {'PPP8'},{7},{5};...
    {'PPP8'},{8},{1};...
    {'PPP8'},{8},{3}];
%% define parameters
savedata = 1;
mindur = 10; % min NREM duration
velfilter = 1; 
win = -30:0.1:30;
%% session loop
for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    ep = animal_info{session_list,3};
    daystring = num2str(day);
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];
    animalname_prefix = [animalname,'-',prefix];

    disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])
    %% reset
    ripple_trigger_pupil = [];
    replayPval_trigger_pupil = [];
    %% load data
    load(fullfile(animaldir,[prefix,'.MergePoints.events.mat'])) % Session info
    epochtimes = MergePoints.timestamps(ep,:);

    % pupil
    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % time, normalized diameter

    % sleep state
    load(fullfile(animaldir,[prefix,'.SleepState.states.mat'])) % sleep states
    Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM
    sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
    Sleep_state_ep = Sleep_state(sleepss_id,:);
    SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
    SWSid = find(Sleep_state_ep(:,2) == 3);
    SWS_vec_ep(SWSid) = 1;
    SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));
    SWSdur = SWSlist(:,2) - SWSlist(:,1);
    SWSlist = SWSlist(find(SWSdur > mindur),:);

    % load ripple, defined by HSEs
    load(sprintf('%s/%srippleHSEs.events.mat', animaldir, prefix));% rippleHSE info
    riptimes = [rippleHSEs{day}{ep}.starttime',rippleHSEs{day}{ep}.endtime'];  

    dur = 1000*(riptimes(:,2) - riptimes(:,1));
    keepidx = find(dur >= 50);%at least 5 bins, 50 ms for 10ms bins
    riplist = riptimes(keepidx,:);
    
    % load replay information
    load(sprintf('%s%sreplayHSEseqdecode_CA1_%02d.mat', replaydir,animalname_prefix,ep));
    eventidx = replaytrajactory{ep}.eventidx;
    candidate_times = riplist(eventidx,:);
    replay_pval = 1 - min(replaytrajactory{ep}.pvalue_r,[],2); % 1 - p_shuf
    %% restrict to NREM
    if ~isempty(SWSlist)
        pupil_trace_id = find(pupil_trace(:,1)>= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
        pupil_trace = pupil_trace(pupil_trace_id,:);
        
        % confine ripples to sws
        riplist_ep = [];
        riplist_pval_ep = [];
        for swsseg = 1:length(SWSlist(:,1))
            validid = find(candidate_times(:,1) >= SWSlist(swsseg,1) & candidate_times(:,2) <= SWSlist(swsseg,2));  % get time stamps
            riplist_ep = [riplist_ep;candidate_times(validid,:)];
            riplist_pval_ep = [riplist_pval_ep;replay_pval(validid)];
        end
        
        %%
        for i = 1:length(riplist_ep)
            rippletime =  riplist_ep(i,:); 
            tempid  = find(pupil_trace(:,1) >= rippletime(1)-40 & pupil_trace(:,1) <= rippletime(1)+40);
            pupil_d = pupil_trace(tempid,2);
            pupil_time = pupil_trace(tempid,1)-rippletime(1);
            if length(find(~isnan(pupil_d))) > 10
                pupil_d_interp = interp1(pupil_time(~isnan(pupil_d)),pupil_d(~isnan(pupil_d)),win);
                ripple_trigger_pupil = [ripple_trigger_pupil;pupil_d_interp];
                replayPval_trigger_pupil = [replayPval_trigger_pupil;riplist_pval_ep(i)];
            end
        end
    end
    %% normalization
    baseline_mean = nanmean(ripple_trigger_pupil(:,1:300),2);
    baseline_std = nanstd(ripple_trigger_pupil(:,1:300)')';
    ripple_triggerZ_pupil = (ripple_trigger_pupil - repmat(baseline_mean,1,length(win)))./repmat(baseline_std,1,length(win));
    ripple_triggerZ_pupil_sum = sum(ripple_triggerZ_pupil,2);
    validid = find(~isnan(ripple_triggerZ_pupil_sum));

    ripple_triggerZ_pupil_mean = nanmean(ripple_triggerZ_pupil(validid,301:end),2);
    replayPval_trigger_pupil = replayPval_trigger_pupil(validid);

    %% save data
    save(fullfile(animaldir,[prefix,'.ReplayPval_Pupil-','EP',num2str(ep),'.mat']),'replayPval_trigger_pupil','ripple_triggerZ_pupil_mean');
end
