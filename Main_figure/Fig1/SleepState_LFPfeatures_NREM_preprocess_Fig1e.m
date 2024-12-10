% this is the preprocessing step to generate Fig. 1e - gathering LFP
% features during NREM sleep for further calculating SI
clc
clear all
close all
%% define data
dir = '/Volumes/Pupil_data/'; % data directory
% list of all data, [{animal},{day},{[epochs]}]
animal_info = [{'HYC2'},{2},{[2,3]};...
    {'HYC3'},{8},{1};...
    {'HYC3'},{9},{1};...
    {'PPP4'},{8},{3};...
    {'PPP4'},{9},{5};...
    {'PPP4'},{10},{5};...
    {'PPP4'},{11},{[1,3,5]};...
    {'PPP4'},{18},{5};...
    {'PPP7'},{8},{[3,5]};...
    {'PPP7'},{12},{[1,4]};...
    {'PPP7'},{14},{1};...
    {'PPP7'},{23},{[3,6]};...
    {'PPP8'},{7},{[1,3,5]};...
    {'PPP8'},{8},{[1,3]}];
%% define parameters
savedata = 1;
mindur = 10; % min NREM duration
window = 2; % LFP smooth window
noverlap = window-1; %1s dt
smoothfact = 15;
% Narrowband theta and delta
f_all = [1 20];
f_theta = [5 10];
f_delta = [1 4];
IRASA = true;
ThIRASA = true;
%% animal-day loop
for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    eps = animal_info{session_list,3};
    daystring = num2str(day);
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];

    load(fullfile(animaldir,[prefix,'.MergePoints.events.mat'])) % Session info

    % pupil info
    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % time, normalized diameter

    % sleep state
    load(fullfile(animaldir,[prefix,'.SleepState.states.mat'])) % sleep states
    Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM

    % LFP features
    load(fullfile(animaldir,[prefix,'.EMGFromLFP.LFP.mat'])) % EMG from LFP
    EMG_LFP = [EMGFromLFP.timestamps,EMGFromLFP.data];

    load(fullfile(animaldir,[prefix,'.SleepScoreLFP.LFP.mat'])) % sleep score LFP, theta and SW (spectrum slope)
    Th_LFP = double(SleepScoreLFP.thLFP);
    Sw_LFP = double(SleepScoreLFP.swLFP);
    %% calculate SW (spectrum slope)
    %Put the LFP in the right structure format
    lfp.data = SleepScoreLFP.swLFP;
    lfp.timestamps = SleepScoreLFP.t;
    lfp.samplingRate = SleepScoreLFP.sf;
    %Calculate PSS, bz_PowerSpectrumSlope.m from AYALab neurocode
    [specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,'frange',[4 90],'IRASA',IRASA);
    broadbandSlowWave = -specslope.data; %So NREM is higher as opposed to lower
    t_clus = specslope.timestamps;
    specdt = 1./specslope.samplingRate;

    %Remove transients before calculating SW histogram
    zFFTspec = NormToInt(spec.amp,'modZ');  % NormToInt from AYALab neurocode
    totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
    badtimes = find(totz>3);

    %Smooth and 0-1 normalize
    broadbandSlowWave(badtimes) = nan;
    broadbandSlowWave = smooth(broadbandSlowWave,smoothfact./specdt);
    %% calculate theta and delta
    %Put the LFP in the right structure format
    lfp.data = SleepScoreLFP.thLFP;
    %Calculate PSS
    [specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
        'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
    t_thclu = specslope.timestamps;
    specdt = 1./specslope.samplingRate;
    thFFTspec = specslope.resid';
    thFFTspec(thFFTspec<0)=0;

    % Remove transients
    zFFTspec = NormToInt(spec.amp,'modZ');
    totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
    badtimes_TH = find(totz>3);

    % theta
    thFFTfreqs = specslope.freqs';
    thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
    thratio = max((thFFTspec(thfreqs,:)),[],1);

    thratio(badtimes_TH) = nan;     %Remove transients
    thratio = smooth(thratio,smoothfact./specdt);

    % delta
    dlfreqs = (thFFTfreqs>=f_delta(1) & thFFTfreqs<=f_delta(2));
    dlratio = max((thFFTspec(dlfreqs,:)),[],1);

    dlratio(badtimes_TH) = nan;     %Remove transients
    dlratio = smooth(dlratio,smoothfact./specdt);
    %% get ripple LFPs
    if strcmp(prefix,'day28') && strcmp(animalname,'PPP4')
        load(fullfile(animaldir,[prefix,'.ripples_task.events.mat'])) % ripple
    elseif strcmp(prefix,'day8') && strcmp(animalname,'PPP7')
        load(fullfile(animaldir,[prefix,'.dorsalripples.events.mat'])) % ripple
    else
        load(fullfile(animaldir,[prefix,'.ripples.events.mat'])) % ripple
    end
    % find ripple channel
    if strcmp(ripples.detectorinfo.detectorname,'FindRipples') % detected without SW
        RipCH = ripples.detectorinfo.detectionparms.channel;
        RipBp = ripples.detectorinfo.detectionparms.passband;
    else % detected with SW
        RipCH = ripples.detectorinfo.detectionparms.Channels(1);
        RipBp = ripples.detectorinfo.detectionparms.ripBP;
    end
    Ripple_lfp = getLFP([RipCH], 'basepath', animaldir,'basename', prefix);

    [specslope,spec] = bz_PowerSpectrumSlope(Ripple_lfp,window,window-noverlap,...
        'nfreqs',200,'frange',[50,400],'IRASA',IRASA);
    specdt = 1./specslope.samplingRate;
    RipFFTspec = specslope.resid';
    RipFFTspec(RipFFTspec<0)=0;

    % Remove transients
    zFFTspec = NormToInt(spec.amp,'modZ');
    totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
    badtimes_RIP = find(totz>3);
    % ripple
    ripFFTfreqs = specslope.freqs';
    ripfreqs = (ripFFTfreqs>=RipBp(1) & ripFFTfreqs<=RipBp(2));
    ripratio = max((RipFFTspec(ripfreqs,:)),[],1);

    ripratio(badtimes_RIP) = nan;     %Remove transients
    ripratio = smooth(ripratio,smoothfact./specdt);
    %% make sure that all the LFP features have the same length
    EMG_LFP_interp = interp1(EMG_LFP(:,1), EMG_LFP(:,2),specslope.timestamps);
    EMG_LFP = EMG_LFP_interp;
    %% epoch loop
    for ep = eps
        disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])
        epochtimes = MergePoints.timestamps(ep,:);

        valididx = find(pupil_trace(:,1) >= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
        pupil_trace_ep = pupil_trace(valididx,:);

        % get NREM episodes
        sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
        Sleep_state_ep = Sleep_state(sleepss_id,:);
        SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
        SWSid = find(Sleep_state_ep(:,2) == 3);
        SWS_vec_ep(SWSid) = 1;
        SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));
        SWSdur = SWSlist(:,2) - SWSlist(:,1);
        SWSlist = SWSlist(find(SWSdur > mindur),:);

        %% restrict LFP features to the current epoch
        valididx = find(specslope.timestamps >= epochtimes(1) & specslope.timestamps <= epochtimes(2));
        LFP_features = [specslope.timestamps(valididx),dlratio(valididx),thratio(valididx),ripratio(valididx),...
            broadbandSlowWave(valididx),EMG_LFP(valididx)]; % time, delta, theta, ripple, SW (or slope), EMG
        LFP_features = double(LFP_features);
        
        % make sure the pupil has the same length as the LFP
        pupil_trace_interp = interp1(pupil_trace_ep(:,1),pupil_trace_ep(:,2),LFP_features(:,1));
        pupil_trace_ep = [LFP_features(:,1),pupil_trace_interp];
        %% restrcit LFP features to NREM periods
        pupil_trace_NREM = [];
        LFP_features_NREM = [];
        if ~isempty(SWSlist)
            for swsseg = 1:length(SWSlist(:,1))
                validid = find(LFP_features(:,1) >= SWSlist(swsseg,1) & LFP_features(:,1) <= SWSlist(swsseg,2));  % get time stamps
                LFP_features_NREM = [LFP_features_NREM;LFP_features(validid,:)];
                pupil_trace_NREM = [pupil_trace_NREM;pupil_trace_ep(validid,:)];
            end
        end
        % zscore pupil
        pupil_label = (pupil_trace_NREM(:,2) - nanmean(pupil_trace_NREM(:,2)))./nanstd(pupil_trace_NREM(:,2));

       if savedata
           save(fullfile(animaldir,[prefix,'.SleepStateLFPfeatures-','EP',num2str(ep),'.LFP.mat']),'LFP_features_NREM','pupil_label');
       end
    end
end
