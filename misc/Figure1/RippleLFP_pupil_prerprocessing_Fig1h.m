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
win = -30:0.1:30;
velfilter = 1; 
%% session loop
for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    ep = animal_info{session_list,3};
    daystring = num2str(day);
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];
    disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])
    %% reset
    ripple_list = [];
    ripple_trigger_pupil = [];
    rippleamp_trigger_pupil = [];
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

    % ripples
    if strcmp(prefix,'day28') && strcmp(animalname,'PPP4')
        load(fullfile(animaldir,[prefix,'.ripples_task.events.mat'])) % ripple
    elseif strcmp(prefix,'day8') && strcmp(animalname,'PPP7')
        load(fullfile(animaldir,[prefix,'.dorsalripples.events.mat'])) % ripple
    else
        load(fullfile(animaldir,[prefix,'.ripples.events.mat'])) % ripple
    end
    riplist_all = [ripples.timestamps(:,1),ripples.timestamps(:,1) + ripples.duration,ripples.amplitude,ripples.peaks]; % [start, end, amplitude, peak time]

    %% restrict to NREM
    if ~isempty(SWSlist)
        pupil_trace_id = find(pupil_trace(:,1)>= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
        pupil_trace = pupil_trace(pupil_trace_id,:);
    
        validid = find(riplist_all(:,1) >= pupil_trace(1,1) & riplist_all(:,2) <= pupil_trace(end,1));  % get time stamps
        riplist_ep = riplist_all(validid,[1:2,4]);
        
        riplist = riplist_ep;
        riplist_amp = riplist(:,3);
            
        % restrict ripples to sws
        riplist_ep = [];
        riplist_amp_ep = [];
        riplist_dur_ep = [];
        for swsseg = 1:length(SWSlist(:,1))
            validid = find(riplist(:,1) >= SWSlist(swsseg,1) & riplist(:,2) <= SWSlist(swsseg,2));  % get time stamps
            riplist_ep = [riplist_ep;riplist(validid,:)];
            riplist_amp_ep = [riplist_amp_ep;riplist_amp(validid)];
        end
            
        % zscore ripple amp
        riplist_ampZ_ep = (riplist_amp_ep - nanmean(riplist_amp_ep))./nanstd(riplist_amp_ep);
        %%
        for i = 1:length(riplist_ep)
            rippletime =  riplist_ep(i,:); 
            tempid  = find(pupil_trace(:,1) >= rippletime(1)-40 & pupil_trace(:,1) <= rippletime(1)+40);
            pupil_d = pupil_trace(tempid,2);
            pupil_time = pupil_trace(tempid,1)-rippletime(1);
            if length(find(~isnan(pupil_d))) > 10
                ripple_list = [ripple_list;rippletime];
                pupil_d_interp = interp1(pupil_time(~isnan(pupil_d)),pupil_d(~isnan(pupil_d)),win);
                ripple_trigger_pupil = [ripple_trigger_pupil;pupil_d_interp];
                rippleamp_trigger_pupil = [rippleamp_trigger_pupil;riplist_ampZ_ep(i)];
            end
        end
    end
    %% normalization
    baseline_mean = nanmean(ripple_trigger_pupil(:,1:300),2);
    baseline_std = nanstd(ripple_trigger_pupil(:,1:300)')';
    ripple_triggerZ_pupil = (ripple_trigger_pupil - repmat(baseline_mean,1,length(win)))./repmat(baseline_std,1,length(win));
    ripple_triggerZ_pupil_sum = sum(ripple_triggerZ_pupil,2);
    validid = find(~isnan(ripple_triggerZ_pupil_sum));


    ripple_list = ripple_list(validid,:);
    ripple_triggerZ_pupil_mean = nanmean(ripple_triggerZ_pupil(validid,301:end),2);
    ripple_triggerZ_pupil = ripple_triggerZ_pupil(validid,:);
    rippleamp_triggerZ_pupil = rippleamp_trigger_pupil(validid);
    %% load ripple LFP 
    % load(fullfile(animaldir,'swrCh.mat'))
    % RipCH = swrCh.Ripple_Channel;
    if strcmp(animalname,'HYC2')
        RipCH = ripples.detectorinfo.detectionparms.channel; % for HYC2
    else
        RipCH = ripples.detectorinfo.detectionparms.Channels(1);
    end
    % Read in xml of file parameters
    par        = LoadXml( [animaldir, prefix, '.xml'] );
    % pull parameters from .xml file
    SR         = par.lfpSampleRate; % lfp sampling rate
    
    Ripple_lfp = getLFP([RipCH], 'basepath', animaldir,'basename', prefix); 
   
    % filter LFP to ripple band
    Ripple_lfp_flit = bz_Filter(Ripple_lfp.data,'passband',[80 inf],'filter','fir1','nyquist',SR/2);
    
    Ripple_lfp.data = Ripple_lfp_flit;

    [lfp_ripples,twin] = eventLFP (Ripple_lfp, ripple_list(:,3),'twin',[0.025,0.025]);
    lfp_ripples = squeeze(lfp_ripples);
    lfp_ripples = lfp_ripples';
    %% save data
    save(fullfile(animaldir,[prefix,'.RippleLFPfeatures-NoRS-','EP',num2str(ep),'.LFP.mat']),'lfp_ripples','rippleamp_triggerZ_pupil','ripple_triggerZ_pupil_mean');
end
