clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wenbotang/Documents/Remote_HMM_files/AYALab_Code/'))
%% setting for the plots
defaultGraphicsSetttings
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';

% list of all data, only post-experience sessions are used
animal_info = [{'PPP4'},{8},{3};...
    {'PPP4'},{10},{5};...
    {'PPP4'},{11},{3};...
    {'PPP4'},{11},{5};...
    {'PPP7'},{8},{3};...
    {'PPP7'},{8},{5};...
    {'PPP7'},{12},{4};...
    {'PPP8'},{7},{3};...
    {'PPP8'},{7},{5};...
    {'PPP8'},{8},{3}];
%% define parameters
velfilter = 1; 
win = -30:0.1:30;
mindur = 10; % min duration for NREM episodes
%% gather data
for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    ep = animal_info{session_list,3};
    daystring = num2str(day);
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];
    
    disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])

    %% load data
    load(fullfile(animaldir,[prefix,'.MergePoints.events.mat'])) % Session info
    
    % sleep state
    load(fullfile(animaldir,[prefix,'.SleepState.states.mat'])) % sleep states
    Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM

    % ripple
    if strcmp(prefix,'day28') && strcmp(animalname,'PPP4')
        load(fullfile(animaldir,[prefix,'.ripples_task.events.mat'])) % ripple
    elseif strcmp(prefix,'day8') && strcmp(animalname,'PPP7')
        load(fullfile(animaldir,[prefix,'.dorsalripples.events.mat'])) % ripple
    else
        load(fullfile(animaldir,[prefix,'.ripples.events.mat'])) % ripple
    end

    riplist_all = [ripples.timestamps(:,1),ripples.timestamps(:,1) + ripples.duration,ripples.amplitude]; % [start, end, amplitude]

    % spikes
    if strcmp(animalname,'PPP4') ||  (strcmp(animalname,'PPP8') && day == 8)
        spikes = importSpikes('basepath',animaldir,'CellType',"Pyramidal Cell");
    elseif strcmp(animalname,'PPP7')
        spikes = importSpikes('basepath',animaldir,'CellType',"Pyramidal Cell",'brainRegion',"dCA1");
    else
        spikes = importSpikes('basepath',animaldir,'CellType',"Pyramidal Cell",'brainRegion',"CA1");
    end
    %% get PRE-sleep SWRs
    epochtimes = MergePoints.timestamps(1,:);
    sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
    Sleep_state_ep = Sleep_state(sleepss_id,:);
    SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
    SWSid = find(Sleep_state_ep(:,2) == 3);
    SWS_vec_ep(SWSid) = 1;
    SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));
    if ~isempty(SWSlist)
        ripple_id = find(riplist_all(:,1) >= epochtimes(1) & riplist_all(:,2) <= epochtimes(2));  % get time stamps
        riplist = riplist_all(ripple_id,1:2);
        riplist_PRE = [];
        for swsseg = 1:length(SWSlist(:,1))
            validid = find(riplist(:,1) >= SWSlist(swsseg,1) & riplist(:,2) <= SWSlist(swsseg,2));  % get time stamps
            riplist_PRE = [riplist_PRE;riplist(validid,:)];
        end
    end
 
    %% get POST-sleep SWRs marked by pupil size
    ripple_trigger_pupil = [];
    ripple_id_pupil = [];
    %% restrict to NREM
    epochtimes = MergePoints.timestamps(ep,:);

    %pupil size
    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % time, normalized diameter

    sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
    Sleep_state_ep = Sleep_state(sleepss_id,:);
    SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
    SWSid = find(Sleep_state_ep(:,2) == 3);
    SWS_vec_ep(SWSid) = 1;
    SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));
    SWSdur = SWSlist(:,2) - SWSlist(:,1);
    SWSlist = SWSlist(find(SWSdur > mindur),:);
    %%
    if ~isempty(SWSlist)
        pupil_trace_id = find(pupil_trace(:,1)>= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
        pupil_trace = pupil_trace(pupil_trace_id,:);

        ripple_id = find(riplist_all(:,1) >= pupil_trace(1,1) & riplist_all(:,2) <= pupil_trace(end,1));  % get time stamps
        riplist = riplist_all(ripple_id,1:2);

        % confine ripples to sws
        riplist_ep = [];
        ripple_id_ep = [];
        riplist_amp_ep = [];
        riplist_dur_ep = [];
        for swsseg = 1:length(SWSlist(:,1))
            validid = find(riplist(:,1) >= SWSlist(swsseg,1) & riplist(:,2) <= SWSlist(swsseg,2));  % get time stamps
            ripple_id_ep = [ripple_id_ep; ripple_id(validid)];
            riplist_ep = [riplist_ep;riplist(validid,:)];
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
                ripple_id_pupil = [ripple_id_pupil;ripple_id_ep(i)];
            end
        end
    end
    %% get the pupil tile info for each ripple
    baseline_mean = nanmean(ripple_trigger_pupil(:,1:300),2);
    baseline_std = nanstd(ripple_trigger_pupil(:,1:300)')';
    ripple_triggerZ_pupil = (ripple_trigger_pupil - repmat(baseline_mean,1,length(win)))./repmat(baseline_std,1,length(win));
    ripple_triggerZ_pupil_sum = sum(ripple_triggerZ_pupil,2);
    validid = find(~isnan(ripple_triggerZ_pupil_sum));

    ripple_triggerZ_pupil = ripple_triggerZ_pupil(validid,:);
    rippleZ_id_pupil = ripple_id_pupil(validid);
    ripple_triggerZ_pupil_mean = nanmean(ripple_triggerZ_pupil(:,401:end),2);

    rippletimes = ripples.timestamps(rippleZ_id_pupil,:);
    %% sorted ripple by pupil size
    pupil_tile = nan(size(ripple_triggerZ_pupil_mean));
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
        pupil_tile(currentid) = i;
    end
    %% get ripple spikes and rank order for each pupil sextiles
    corr_all = [];
    for tile = 1:6
        ripIDs = find(pupil_tile == tile);
        riplist_all = [riplist_PRE;rippletimes(ripIDs,:)];
        riplist_labels = [zeros(size(riplist_PRE(:,1)));ones(size(ripIDs))];
        spikes_ripples = getRipSpikes(spikes,riplist_all,'basepath',animaldir,'saveMat',false);
        temp = RankOrder_noshuffle('basepath',animaldir,'spkEventTimes',spikes_ripples,'eventIDs',riplist_labels);
        rankcorr{tile} = temp.corrEvents;
        corr_all = [corr_all, temp.corrEvents];
    end

    for tile = 1:6
        % zscore correlation coefficients for each session
        rankcorr{tile} = (rankcorr{tile} - nanmean(corr_all))./nanstd(corr_all);
        if session_list == 1
            rankcorr_all{tile} = rankcorr{tile};
        else
            rankcorr_all{tile} = [rankcorr_all{tile},rankcorr{tile}];
        end
    end
end
    




