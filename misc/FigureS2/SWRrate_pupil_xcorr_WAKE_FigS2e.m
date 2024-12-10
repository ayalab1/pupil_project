clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';
% list of all data
animal_info = [{'HYC2'},{2},{2};...
    {'HYC2'},{2},{3};...
    {'HYC3'},{8},{1};...
    {'HYC3'},{9},{1};...
    {'PPP4'},{8},{3};...
    {'PPP4'},{9},{5};...
    {'PPP4'},{10},{5};...
    {'PPP4'},{11},{1};...
    {'PPP4'},{11},{3};...
    {'PPP4'},{11},{5};...
    {'PPP4'},{18},{5};...
    {'PPP7'},{8},{3};...
    {'PPP7'},{8},{5};...
    {'PPP7'},{12},{1};...
    {'PPP7'},{12},{4};...
    {'PPP7'},{14},{1};...
    {'PPP7'},{23},{3};...
    {'PPP7'},{23},{6};...
    {'PPP8'},{7},{1};...
    {'PPP8'},{7},{3};...
    {'PPP8'},{7},{5};...
    {'PPP8'},{8},{1};...
    {'PPP8'},{8},{3}];
%% define parameters
min_duration = 30; %s
%% gather data 
ripple_pupil_corr_all = [];
for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    ep = animal_info{session_list,3};
    daystring = num2str(day);
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];
    disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])
    
    % load data
    load(fullfile(animaldir,[prefix,'.MergePoints.events.mat'])) % Session info
    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % time, normalized diameter

    if strcmp(prefix,'day28') && strcmp(animalname,'PPP4')
        load(fullfile(animaldir,[prefix,'.ripples_task.events.mat'])) % ripple
    elseif strcmp(prefix,'day8') && strcmp(animalname,'PPP7')
        load(fullfile(animaldir,[prefix,'.dorsalripples.events.mat'])) % ripple
    else
        load(fullfile(animaldir,[prefix,'.ripples.events.mat'])) % ripple
    end
    
    load(fullfile(animaldir,[prefix,'.SleepState.states.mat'])) % sleep states
    Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM
    %% get pupil info
    epochtimes = MergePoints.timestamps(ep,:); % epoch start and end time

    % sleep state
    sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
    Sleep_state_ep = Sleep_state(sleepss_id,:);
    
    WAKE_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
    WAKEid = find(Sleep_state_ep(:,2) == 1);
    WAKE_vec_ep(WAKEid) = 1;
    WAKElist = vec2list(WAKE_vec_ep, Sleep_state_ep(:,1));
    WAKEdur = WAKElist(:,2) - WAKElist(:,1);
    WAKElist = WAKElist(find(WAKEdur > min_duration),:);
  
    % Ripple rate
    ripple_id = find(ripples.timestamps(:,1) <= epochtimes(2) & ripples.timestamps(:,1) >= epochtimes(1));
    ripples_ep = ripples.timestamps(ripple_id,1);
    ripple_events(1).times = ripples_ep;
    % smoothed PSTH
    [R,t,~] = psth(ripple_events,1);

    % pupil diameter
    tempid  = find(pupil_trace(:,1) >= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
    pupil_d = pupil_trace(tempid,2);
    pupil_time = pupil_trace(tempid,1);
    pupil_SR = 1/nanmedian(diff(pupil_trace(:,1)));
    pupil_d_interp = interp1(pupil_time(~isnan(pupil_d)),pupil_d(~isnan(pupil_d)),pupil_time);
    % smooth pupil trace
    pupil_d_fit=locsmooth(pupil_d_interp,pupil_SR,1);
    %% calculate cross-correlation
    R_WAKE = [];
    pupil_d_fit_WAKE = [];
    for seg = 1:length(WAKElist(:,1))
        currentseg = WAKElist(seg,:);
        segid = find(pupil_time >= currentseg(1) & pupil_time <= currentseg(2));
        pupil_d_fit_WAKE = [pupil_d_fit_WAKE, pupil_d_fit(segid)];
        segid = find(t >= currentseg(1) & t <= currentseg(2));
        R_WAKE = [R_WAKE, R(segid)];
    end

    % zscore
    R = (R-nanmean(R_WAKE))/nanstd(R_WAKE);
    pupil_d_fit = (pupil_d_fit-nanmean(pupil_d_fit_WAKE))/nanstd(pupil_d_fit_WAKE);

    pupil_d_fit(isnan(pupil_d_fit)) = 0;

    ripplerate_interp = interp1(t,R,pupil_time);
    ripplerate_interp(isnan(ripplerate_interp)) = 0;
    
    for seg = 1:length(WAKElist(:,1))
        currentseg = WAKElist(seg,:);
        segid = find(pupil_time >= currentseg(1) & pupil_time <= currentseg(2));
        if ~isempty(segid)
            validid = find(~isnan(pupil_d_fit(segid)));
            pupil_d_fit_seg = pupil_d_fit(segid(validid));
            ripplerate_seg = ripplerate_interp(segid(validid));
            
            [ripple_pupil_corr,lags] = xcorr(pupil_d_fit_seg,ripplerate_seg,round(30*30),'biased');
            ripple_pupil_corr_all = [ripple_pupil_corr_all;ripple_pupil_corr];
        end
    end
    close all
end
%%
lags = lags/30;
ripple_trigger_pupil_stats = [lags',nanmean(ripple_pupil_corr_all)',nanstd(ripple_pupil_corr_all)',sum(double((~isnan(ripple_pupil_corr_all))))'];
figure,plot(ripple_trigger_pupil_stats(:,1),ripple_trigger_pupil_stats(:,2))

    

    
    
    
    
