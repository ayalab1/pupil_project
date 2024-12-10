clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wenbotang/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';

% list of all data
animal_info = [{'HYC3'},{8},{1};...
    {'HYC3'},{9},{1};...
    {'PPP4'},{8},{3};...
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
velfilter = 1; 
win = -30:0.1:30;
%% session loop
ripple_trigger_pupil = []; % reset
reactivation_trigger_pupil = [];
for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    ep = animal_info{session_list,3};
    daystring = num2str(day);
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];
    
    disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])

    % pupil info
    load(fullfile(animaldir,[prefix,'.MergePoints.events.mat'])) % Session info
    epochtimes = MergePoints.timestamps(ep,:);

    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % time, normalized diameter

    % sleep states
    load(fullfile(animaldir,[prefix,'.SleepState.states.mat'])) % sleep states
    Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM
    sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
    Sleep_state_ep = Sleep_state(sleepss_id,:);
    WAKE_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
    WAKEid = find(Sleep_state_ep(:,2) == 1); %WAKE state
    WAKE_vec_ep(WAKEid) = 1;
    WAKElist = vec2list(WAKE_vec_ep, Sleep_state_ep(:,1));
    WAKEdur = WAKElist(:,2) - WAKElist(:,1);
    WAKElist = WAKElist(find(WAKEdur > 0),:);

    % ripples 
    if strcmp(prefix,'day28') && strcmp(animalname,'PPP4')
        load(fullfile(animaldir,[prefix,'.ripples_task.events.mat'])) % ripple
    elseif strcmp(prefix,'day8') && strcmp(animalname,'PPP7')
        load(fullfile(animaldir,[prefix,'.dorsalripples.events.mat'])) % ripple
    else
        load(fullfile(animaldir,[prefix,'.ripples.events.mat'])) % ripple
    end
    
    riplist_all = [ripples.timestamps(:,1),ripples.timestamps(:,1) + ripples.duration,ripples.amplitude]; % [start, end, amplitude]
    
    % reactivation strength
    load(fullfile(animaldir,[prefix,'ReactivationStrength_all.mat']))
    %% restrict to WAKE
    if ~isempty(WAKElist)
        pupil_trace_id = find(pupil_trace(:,1)>= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
        pupil_trace = pupil_trace(pupil_trace_id,:);

        validid = find(riplist_all(:,1) >= pupil_trace(1,1) & riplist_all(:,2) <= pupil_trace(end,1));  % get time stamps
        riplist_ep = riplist_all(validid,1:2);
        riplist_RS = [];
        riplist = [];
        for rip = 1:length(validid)
            ripid = find(strength(:,1) >= riplist_ep(rip,1) & strength(:,1) <= riplist_ep(rip,2));
            if ~isempty(ripid) && ~isempty(MUAid)
                riplist_RS = [riplist_RS;nanmean(strength(ripid,2))];
                riplist = [riplist;riplist_ep(rip,:)];
            end
        end
                    
        % confine ripples to WAKE
        riplist_ep = [];
        riplist_RS_ep = [];
        for WAKEseg = 1:length(WAKElist(:,1))
            validid = find(riplist(:,1) >= WAKElist(WAKEseg,1) & riplist(:,2) <= WAKElist(WAKEseg,2));  % get time stamps
            riplist_ep = [riplist_ep;riplist(validid,:)];
            riplist_RS_ep = [riplist_RS_ep;riplist_RS(validid)];
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
                reactivation_trigger_pupil = [reactivation_trigger_pupil;riplist_RS_ep(i)];
            end
        end
    end
end
%% normalized to baseline
reactivation_triggerZ_pupil = (reactivation_trigger_pupil - nanmean(reactivation_trigger_pupil))./nanstd(reactivation_trigger_pupil);
baseline_mean = nanmean(ripple_trigger_pupil(:,1:300),2);
baseline_std = nanstd(ripple_trigger_pupil(:,1:300)')';
% zscore
ripple_triggerZ_pupil = (ripple_trigger_pupil - repmat(baseline_mean,1,length(win)))./repmat(baseline_std,1,length(win));

ripple_triggerZ_pupil_sum = sum(ripple_triggerZ_pupil,2);
validid = find(~isnan(ripple_triggerZ_pupil_sum));

ripple_triggerZ_pupil = ripple_triggerZ_pupil(validid,:);
reactivation_triggerZ_pupil = reactivation_triggerZ_pupil(validid);
ripple_triggerZ_pupil_mean = nanmean(ripple_triggerZ_pupil(:,301:end),2);
%% sorted RS by pupil size
Y = prctile(ripple_triggerZ_pupil_mean,1/6*100);
tile1 = [min(ripple_triggerZ_pupil_mean),Y];
tile2 = [Y,prctile(ripple_triggerZ_pupil_mean,2/6*100)];
tile3 = [prctile(ripple_triggerZ_pupil_mean,2/6*100),prctile(ripple_triggerZ_pupil_mean,3/6*100)];
tile4 = [prctile(ripple_triggerZ_pupil_mean,3/6*100),prctile(ripple_triggerZ_pupil_mean,4/6*100)];
tile5 = [prctile(ripple_triggerZ_pupil_mean,4/6*100),prctile(ripple_triggerZ_pupil_mean,5/6*100)];
tile6 = [prctile(ripple_triggerZ_pupil_mean,5/6*100),prctile(ripple_triggerZ_pupil_mean,6/6*100)];
sixtiles = [tile1;tile2;tile3;tile4;tile5;tile6];

figure
for i = 1:6
    currentid = find(ripple_triggerZ_pupil_mean >= sixtiles(i,1) & ripple_triggerZ_pupil_mean <= sixtiles(i,2));
    reactivation_tile{i} = reactivation_triggerZ_pupil(currentid);
end
