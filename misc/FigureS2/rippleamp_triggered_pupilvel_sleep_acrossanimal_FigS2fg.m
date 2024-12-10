clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';
% data list
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
%% set parameters
velfilter = 1; % cm/s; running velocity threshold
win = -30:0.1:30;
State = 1; %1 awake, 3 NREM, 5 REM
%% gather data
ripple_trigger_pupil = [];
rippleamp_trigger_pupil = [];

for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    ep = animal_info{session_list,3};
    daystring = num2str(day);
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];
    
    disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])

    % pupil info
    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % time, normalized diameter
    epochtimes = [behavior.epochs{ep}.startTime,behavior.epochs{ep}.stopTime];  

    % sleepstate info
    load(fullfile(animaldir,[prefix,'.SleepState.states.mat'])) % sleep states
    Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM
    sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
    Sleep_state_ep = Sleep_state(sleepss_id,:);
    SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
    SWSid = find(Sleep_state_ep(:,2) == State);
    SWS_vec_ep(SWSid) = 1;
    SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));
    SWSdur = SWSlist(:,2) - SWSlist(:,1);
    SWSlist = SWSlist(find(SWSdur > 10),:); % episodes longer than 10s
    
    % ripple info
    if strcmp(prefix,'day28') && strcmp(animalname,'PPP4')
        load(fullfile(animaldir,[prefix,'.ripples_task.events.mat'])) % ripple
    elseif strcmp(prefix,'day8') && strcmp(animalname,'PPP7')
        load(fullfile(animaldir,[prefix,'.dorsalripples.events.mat'])) % ripple
    else
        load(fullfile(animaldir,[prefix,'.ripples.events.mat'])) % ripple
    end
   
    riplist_all = [ripples.timestamps(:,1),ripples.timestamps(:,1) + ripples.duration,ripples.amplitude]; % [start, end, amplitude]

    %% restrict data to the sleep state
    if ~isempty(SWSlist)
        pupil_trace_id = find(pupil_trace(:,1)>= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
        pupil_trace = pupil_trace(pupil_trace_id,:);

        validid = find(riplist_all(:,1) >= pupil_trace(1,1) & riplist_all(:,2) <= pupil_trace(end,1));  % get time stamps
        riplist = riplist_all(validid,1:2);
        riplist_amp = riplist_all(validid,3);
        
        % restrict ripples
        riplist_ep = [];
        riplist_amp_ep = [];

        for swsseg = 1:length(SWSlist(:,1))
            validid = find(riplist(:,1) >= SWSlist(swsseg,1) & riplist(:,2) <= SWSlist(swsseg,2));  % get time stamps
            riplist_ep = [riplist_ep;riplist(validid,:)];
            riplist_amp_ep = [riplist_amp_ep;riplist_amp(validid)];
        end
        
        % zscore ripple amplitude
        riplist_ampZ_ep = (riplist_amp_ep - nanmean(riplist_amp_ep))./nanstd(riplist_amp_ep);
        %%
        for i = 1:length(riplist_ep)
            rippletime =  riplist_ep(i,:); 
            tempid  = find(pupil_trace(:,1) >= rippletime(1)-40 & pupil_trace(:,1) <= rippletime(1)+40);
            pupil_d = pupil_trace(tempid,2);
            pupil_time = pupil_trace(tempid,1)-rippletime(1);
            if length(find(~isnan(pupil_d))) > 10
                pupil_d_interp = interp1(pupil_time(~isnan(pupil_d)),pupil_d(~isnan(pupil_d)),win);
                ripple_trigger_pupil = [ripple_trigger_pupil;pupil_d_interp];
                rippleamp_trigger_pupil = [rippleamp_trigger_pupil;riplist_ampZ_ep(i)];
            end
        end
    end
end
%% gather data
baseline_mean = nanmean(ripple_trigger_pupil(:,1:300),2);
baseline_std = nanstd(ripple_trigger_pupil(:,1:300)')';
ripple_triggerZ_pupil = (ripple_trigger_pupil - repmat(baseline_mean,1,length(win)))./repmat(baseline_std,1,length(win));
ripple_triggerZ_pupil_sum = sum(ripple_triggerZ_pupil,2);
validid = find(~isnan(ripple_triggerZ_pupil_sum));

ripple_triggerZ_pupil = ripple_triggerZ_pupil(validid,:);
rippleamp_triggerZ_pupil = rippleamp_trigger_pupil(validid);
ripple_triggerZ_pupil_mean = nanmean(ripple_triggerZ_pupil(:,301:end),2);
%% correlation 
r = corr(ripple_triggerZ_pupil_mean,rippleamp_triggerZ_pupil,'type','spearman');
%% sorted pupil by ripple amplitude
Y = prctile(rippleamp_triggerZ_pupil,1/6*100);
tile1 = [min(rippleamp_triggerZ_pupil),Y];
tile2 = [Y,prctile(rippleamp_triggerZ_pupil,2/6*100)];
tile3 = [prctile(rippleamp_triggerZ_pupil,2/6*100),prctile(rippleamp_triggerZ_pupil,3/6*100)];
tile4 = [prctile(rippleamp_triggerZ_pupil,3/6*100),prctile(rippleamp_triggerZ_pupil,4/6*100)];
tile5 = [prctile(rippleamp_triggerZ_pupil,4/6*100),prctile(rippleamp_triggerZ_pupil,5/6*100)];
tile6 = [prctile(rippleamp_triggerZ_pupil,5/6*100),prctile(rippleamp_triggerZ_pupil,6/6*100)];
sixtiles = [tile1;tile2;tile3;tile4;tile5;tile6];
figure
for i = 1:6
    currentid = find(rippleamp_triggerZ_pupil >= sixtiles(i,1) & rippleamp_triggerZ_pupil <= sixtiles(i,2));
    ripple_triggerZ_pupil_tile = ripple_triggerZ_pupil(currentid,:);
    ripple_trigger_pupil_stats{i} = [win',nanmean(ripple_triggerZ_pupil_tile)',nanstd(ripple_triggerZ_pupil_tile)',sum(double((~isnan(ripple_triggerZ_pupil_tile))))'];
    plot(win',nanmean(ripple_triggerZ_pupil_tile)')
    hold on
end
%% sorted ripple amplitude by pupil
Y = prctile(ripple_triggerZ_pupil_mean,1/6*100);
tile1 = [min(ripple_triggerZ_pupil_mean),Y];
tile2 = [Y,prctile(ripple_triggerZ_pupil_mean,2/6*100)];
tile3 = [prctile(ripple_triggerZ_pupil_mean,2/6*100),prctile(ripple_triggerZ_pupil_mean,3/6*100)];
tile4 = [prctile(ripple_triggerZ_pupil_mean,3/6*100),prctile(ripple_triggerZ_pupil_mean,4/6*100)];
tile5 = [prctile(ripple_triggerZ_pupil_mean,4/6*100),prctile(ripple_triggerZ_pupil_mean,5/6*100)];
tile6 = [prctile(ripple_triggerZ_pupil_mean,5/6*100),prctile(ripple_triggerZ_pupil_mean,6/6*100)];
sixtiles = [tile1;tile2;tile3;tile4;tile5;tile6];

for i = 1:6
    currentid = find(ripple_triggerZ_pupil_mean >= sixtiles(i,1) & ripple_triggerZ_pupil_mean <= sixtiles(i,2));
    rippleamp_tile{i} = rippleamp_triggerZ_pupil(currentid);
end