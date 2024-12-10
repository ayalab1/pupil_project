% this is the preprocessing step to generate Fig. 2c and 2e
clc
clear all
close all
%% define data
% directory saved all the replay decoding results
% see replay decoding subfolder for more info: .../src/ReplayDecoding/
replaydir = '/Users/wt248/Replay_decoding_20230526/';
dir = '/Volumes/Pupil_data/';

% list of all data,[{animal},{day},{epoch}]
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
velfilter = 1; 
SWS_mindur = 10; % 10s
win = -30:0.1:30; % window around SWRs, -30s to 30s
%% session loop
ripple_info = [];
pupil_info = [];
for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    ep = animal_info{session_list,3};
    daystring = num2str(day);
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];
    animalname_prefix = [animalname,'-',prefix];
    
    disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])

    load(fullfile(animaldir,[prefix,'.MergePoints.events.mat'])) % Session info
    epochtimes = MergePoints.timestamps(ep,:);

    % load pupil info
    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % [time, normalized pupil diameter]

    % load sleep states
    load(fullfile(animaldir,[prefix,'.SleepState.states.mat'])) % sleep states
    Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM
    sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
    Sleep_state_ep = Sleep_state(sleepss_id,:);
    SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
    SWSid = find(Sleep_state_ep(:,2) == 3); %NREM
    SWS_vec_ep(SWSid) = 1;
    SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));
    SWSdur = SWSlist(:,2) - SWSlist(:,1);
    SWSlist = SWSlist(find(SWSdur > SWS_mindur),:); % last for a min duration
    
    % load ripple information
    if strcmp(prefix,'day28') && strcmp(animalname,'PPP4')
        load(fullfile(animaldir,[prefix,'.ripples_task.events.mat'])) % ripple
    elseif strcmp(prefix,'day8') && strcmp(animalname,'PPP7')
        load(fullfile(animaldir,[prefix,'.dorsalripples.events.mat'])) % ripple
    else
        load(fullfile(animaldir,[prefix,'.ripples.events.mat'])) % ripple
    end
    
    riplist_all = [ripples.timestamps(:,1),ripples.timestamps(:,1) + ripples.duration,ripples.duration, ripples.amplitude]; % [start, end, duration, amplitude]
    
    % load ripple-crossed HSE ripples
    load(sprintf('%s/%srippleHSEs.events.mat', animaldir, prefix));% rippleHSE info
    riptimes = [rippleHSEs{day}{ep}.starttime',rippleHSEs{day}{ep}.endtime'];  

    dur = 1000*(riptimes(:,2) - riptimes(:,1));
    keepidx = find(dur >= 50);%at least 5 bins, 50 ms for 10ms bins
    ripHSElist = riptimes(keepidx,:);
    
    % load replay information
    load(sprintf('%s%sreplayHSEseqdecode_CA1_%02d.mat', replaydir,animalname_prefix,ep));
    eventidx = replaytrajactory{ep}.eventidx;
    candidate_times = ripHSElist(eventidx,:);
    replay_pval = min(replaytrajactory{ep}.pvalue_r,[],2); %p_shuf
    
    HSElist_all = [candidate_times,replay_pval];
    
    % load reactivation strength
    % see .../src/ReactivationStrength/ for more info on calculating RS
    load(fullfile(animaldir,[prefix,'ReactivationStrength.mat']))  
    
    % load MUA firing rates
    filename = [animaldir,'MUA', '0',num2str(day),'-0',num2str(ep),'.mat'];
    load(filename)
    MUA_FR = [MUA.time',MUA.data(:,1)];
    
    % match replay HSEs and reactivation strength to ripples
    for event = 1:length(riplist_all(:,1))
        temp_RS = find(strength(:,1) <= riplist_all(event,2) & strength(:,1) >= riplist_all(event,1));

        if ~isempty(temp_RS)
            riplist_all(event,5) = nanmean(strength(temp_RS,2));
        else
            riplist_all(event,5) = 0;
        end
        
        temp_MUA = find(MUA_FR(:,1) <= riplist_all(event,2) & MUA_FR(:,1) >= riplist_all(event,1));        
        if ~isempty(temp_MUA)
            riplist_all(event,6) = nanmean(MUA_FR(temp_MUA,2));
        else
            riplist_all(event,6) = 0;
        end
            
        temp = find(HSElist_all(:,1) <= riplist_all(event,1) & HSElist_all(:,2) >= riplist_all(event,1));
        if ~isempty(temp)
            riplist_all(event,7) = HSElist_all(temp,3); % [start, end, dur, amplitude, RS, MUA FR, replay pval]
        else
            riplist_all(event,7) = 1; %pval == 1, max
        end
    end 
          
    %% restrict to NREM
    if ~isempty(SWSlist)  % epoch with NREM episodes
        pupil_trace_id = find(pupil_trace(:,1)>= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
        pupil_trace = pupil_trace(pupil_trace_id,:);
        
        validid = find(riplist_all(:,1) >= pupil_trace(1,1) & riplist_all(:,2) <= pupil_trace(end,1));  % get time stamps
        riplist = riplist_all(validid,:);
       
        % restrict ripples to NREM
        riplist_ep = [];
        for swsseg = 1:length(SWSlist(:,1))
            validid = find(riplist(:,1) >= SWSlist(swsseg,1) & riplist(:,2) <= SWSlist(swsseg,2));  % get time stamps
            riplist_ep = [riplist_ep;riplist(validid,:)];
        end
        
        % zscore ripple duration, amplitude, and RS
        riplist_durZ_ep = (riplist_ep(:,3) - nanmean(riplist_ep(:,3)))./nanstd(riplist_ep(:,3));
        riplist_ep(:,3) = riplist_durZ_ep;
        
        riplist_ampZ_ep = (riplist_ep(:,4) - nanmean(riplist_ep(:,4)))./nanstd(riplist_ep(:,4));
        riplist_ep(:,4) = riplist_ampZ_ep;
        
        riplist_RSZ_ep = (riplist_ep(:,5) - nanmean(riplist_ep(:,5)))./nanstd(riplist_ep(:,5));
        riplist_ep(:,5) = riplist_RSZ_ep;
        
        riplist_MUAZ_ep = (riplist_ep(:,6) - nanmean(riplist_ep(:,6)))./nanstd(riplist_ep(:,6));
        riplist_ep(:,6) = riplist_MUAZ_ep;

        %% calculated pupil size around ripples
        for i = 1:length(riplist_ep)
            rippletime =  riplist_ep(i,:); 
            tempid  = find(pupil_trace(:,1) >= rippletime(1)-40 & pupil_trace(:,1) <= rippletime(1)+40);
            pupil_d = pupil_trace(tempid,2);
            pupil_time = pupil_trace(tempid,1)-rippletime(1);
            if length(find(~isnan(pupil_d))) > 10
                pupil_d_interp = interp1(pupil_time(~isnan(pupil_d)),pupil_d(~isnan(pupil_d)),win);
                pupil_info = [pupil_info;pupil_d_interp];
                ripple_info = [ripple_info;riplist_ep(i,3:end)];
            end
        end
    end
end
%% normalization
baseline_mean = nanmean(pupil_info(:,1:300),2);
baseline_std = nanstd(pupil_info(:,1:300)')';
pupil_infoZ = (pupil_info - repmat(baseline_mean,1,length(win)))./repmat(baseline_std,1,length(win));
pupil_infoZ_sum = sum(pupil_infoZ,2);
validid = find(~isnan(pupil_infoZ_sum));

pupil_info = pupil_infoZ(validid,:);
ripple_info = ripple_info(validid,:);
%% save data
save('GLMgain_SWRproperties_pupil_all.mat','pupil_info','ripple_info');
