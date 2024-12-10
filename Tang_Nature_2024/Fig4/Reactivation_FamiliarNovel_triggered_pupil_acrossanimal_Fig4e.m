% this is the main script to generate Fig. 4e
clc
clear all
close all
%% define data
dir = '/Volumes/Pupil_data/';

% recent  = novel; remote = familiar
% data list, [animal, day, POST epoch, familiar-novel session pairs]
animal_info = [{'PPP7'},{12},{6},{[3,5]};...
    {'PPP7'},{14},{5},{[2,4]};...
    {'PPP8'},{12},{5},{[2,4]};...
    {'PPP8'},{15},{5},{[4,2]};...
    {'PPP8'},{14},{5},{[2,4]};...
    {'PPP8'},{17},{5},{[4,2]};...
    {'PPP15'},{8},{4},{[2,3]};...
    {'PPP13'},{9},{4},{[2,3]};...
    {'PPP13'},{12},{4},{[2,3]};...
    {'PVR4'},{1},{4},{[2,3]}];
%% set parameters
win = -30:0.1:30; % time window in s
%% gather data
ripple_trigger_pupil = [];
reactivationNovel_trigger_pupil = [];
reactivationFamiliar_trigger_pupil = [];
for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    ep = animal_info{session_list,3};
    eps_RUN = animal_info{session_list,4};
    
    if  strcmp(animalname,'PVR4')
        daystring = num2str(day);
        prefix = ['pupil',daystring];
        animaldir = [dir,animalname,'/pupil',daystring,'/'];
    else
        daystring = num2str(day);
        animaldir = [dir,animalname,'/day',daystring,'/'];
        prefix = ['day',daystring];
    end
    disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])

    % load data
    load(fullfile(animaldir,[prefix,'.MergePoints.events.mat'])) % Session info
    epochtimes = MergePoints.timestamps(ep,:);

    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % [time, normalized pupil diameter]

    load(fullfile(animaldir,[prefix,'.SleepState.states.mat'])) % sleep states
    Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM
    sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
    Sleep_state_ep = Sleep_state(sleepss_id,:);
    SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
    SWSid = find(Sleep_state_ep(:,2) == 3); % NREM
    SWS_vec_ep(SWSid) = 1;
    SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));
    SWSdur = SWSlist(:,2) - SWSlist(:,1);
    SWSlist = SWSlist(find(SWSdur > 10),:);

    % load ripples
    if strcmp(prefix,'day28') && strcmp(animalname,'PPP4')
        load(fullfile(animaldir,[prefix,'.ripples_task.events.mat'])) % ripple
    elseif strcmp(prefix,'day8') && strcmp(animalname,'PPP7')
        load(fullfile(animaldir,[prefix,'.dorsalripples.events.mat'])) % ripple
    else
        load(fullfile(animaldir,[prefix,'.ripples.events.mat'])) % ripple
    end
    
    riplist_all = [ripples.timestamps(:,1),ripples.timestamps(:,1) + ripples.duration,ripples.amplitude]; % [start, end, amplitude]
    
    % load reactivation strength (normalized by RS in PRE sleep)
    load(fullfile(animaldir,[prefix,'ReactivationStrength_diffTemp_PREnorm.mat']))
    strength_familiar = strength_all{eps_RUN(1)};
    strength_novel = strength_all{eps_RUN(2)};
    %% restrict to NREM
    if ~isempty(SWSlist)
        pupil_trace_id = find(pupil_trace(:,1)>= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
        pupil_trace = pupil_trace(pupil_trace_id,:);

        validid = find(riplist_all(:,1) >= pupil_trace(1,1) & riplist_all(:,2) <= pupil_trace(end,1));  % get time stamps
        riplist_ep = riplist_all(validid,1:2);
        riplist_RS_recent = [];
        riplist_RS_remote = [];

        riplist = [];
        for rip = 1:length(validid)
            ripid1 = find(strength_novel(:,1) >= riplist_ep(rip,1) & strength_novel(:,1) <= riplist_ep(rip,2));
            ripid2 = find(strength_familiar(:,1) >= riplist_ep(rip,1) & strength_familiar(:,1) <= riplist_ep(rip,2));
            if ~isempty(ripid1) && ~isempty(ripid2)
                riplist_RS_recent = [riplist_RS_recent;nanmean(strength_novel(ripid1,2))];
                riplist_RS_remote = [riplist_RS_remote;nanmean(strength_familiar(ripid2,2))];
                riplist = [riplist;riplist_ep(rip,:)];
            end
        end
        
        % restrict ripples to NREM
        riplist_ep = [];
        riplist_RS_recent_ep = [];
        riplist_RS_remote_ep = [];
        for swsseg = 1:length(SWSlist(:,1))
            validid = find(riplist(:,1) >= SWSlist(swsseg,1) & riplist(:,2) <= SWSlist(swsseg,2));  % get time stamps
            riplist_ep = [riplist_ep;riplist(validid,:)];
            riplist_RS_recent_ep = [riplist_RS_recent_ep;riplist_RS_recent(validid)];
            riplist_RS_remote_ep = [riplist_RS_remote_ep;riplist_RS_remote(validid)];
        end
        
        % zscore RS
        riplist_RSZ_recent_ep = (riplist_RS_recent_ep - nanmean(riplist_RS_recent_ep))./nanstd(riplist_RS_recent_ep);
        riplist_RSZ_remote_ep = (riplist_RS_remote_ep - nanmean(riplist_RS_remote_ep))./nanstd(riplist_RS_remote_ep);
        %% get RS triggered pupil info
        for i = 1:length(riplist_ep)
            rippletime =  riplist_ep(i,:); 
            tempid  = find(pupil_trace(:,1) >= rippletime(1)-40 & pupil_trace(:,1) <= rippletime(1)+40);
            pupil_d = pupil_trace(tempid,2);
            pupil_time = pupil_trace(tempid,1)-rippletime(1);
            if length(find(~isnan(pupil_d))) > 10
                pupil_d_interp = interp1(pupil_time(~isnan(pupil_d)),pupil_d(~isnan(pupil_d)),win);
                ripple_trigger_pupil = [ripple_trigger_pupil;pupil_d_interp];
                reactivationNovel_trigger_pupil = [reactivationNovel_trigger_pupil;riplist_RSZ_recent_ep(i)];
                reactivationFamiliar_trigger_pupil = [reactivationFamiliar_trigger_pupil;riplist_RSZ_remote_ep(i)];
            end
        end
    end
end
%% normalize by the baseline
baseline_mean = nanmean(ripple_trigger_pupil(:,1:300),2);
baseline_std = nanstd(ripple_trigger_pupil(:,1:300)')';
ripple_triggerZ_pupil = (ripple_trigger_pupil - repmat(baseline_mean,1,length(win)))./repmat(baseline_std,1,length(win));
ripple_triggerZ_pupil_sum = sum(ripple_triggerZ_pupil,2);
validid = find(~isnan(ripple_triggerZ_pupil_sum));

ripple_triggerZ_pupil = ripple_triggerZ_pupil(validid,:);
reactivationNovel_triggerZ_pupil = reactivationNovel_trigger_pupil(validid);
reactivationFamiliar_triggerZ_pupil = reactivationFamiliar_trigger_pupil(validid);
% zscore
reactivationNovel_triggerZ_pupil = (reactivationNovel_triggerZ_pupil - nanmean(reactivationNovel_triggerZ_pupil))./nanstd(reactivationNovel_triggerZ_pupil);
reactivationFamiliar_triggerZ_pupil = (reactivationFamiliar_triggerZ_pupil - nanmean(reactivationFamiliar_triggerZ_pupil))./nanstd(reactivationFamiliar_triggerZ_pupil);
ripple_triggerZ_pupil_mean = nanmean(ripple_triggerZ_pupil(:,301:end),2);

reactivation_Index_pupil = ( reactivationFamiliar_trigger_pupil(validid) -  reactivationNovel_trigger_pupil(validid));
reactivation_Index_pupil = (reactivation_Index_pupil-nanmean(reactivation_Index_pupil))./nanstd(reactivation_Index_pupil); % zscore
%% sort reactivation index by pupil size
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
    reactivation_Index_tile{i} = reactivation_Index_pupil(currentid);
end