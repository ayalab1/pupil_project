% this is the main script to generate Fig. 2b
clc
clear all
close all
%% define data
% directory save all the replay decoding results
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
win = -30:0.1:30; % window around SWRs, -30s to 30s
%% session loop
ripple_trigger_pupil = []; % reset
replayPval_trigger_pupil = [];
MUA_trigger_pupil = [];
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

    % pupil info
    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % [time, normalized pupil diameter]

    % sleep states
    load(fullfile(animaldir,[prefix,'.SleepState.states.mat'])) % sleep states
    Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM
    sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
    Sleep_state_ep = Sleep_state(sleepss_id,:);
    SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
    SWSid = find(Sleep_state_ep(:,2) == 3);
    SWS_vec_ep(SWSid) = 1;
    SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));
    SWSdur = SWSlist(:,2) - SWSlist(:,1);
    SWSlist = SWSlist(find(SWSdur > 10),:);

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
    replay_pval = 1 - min(replaytrajactory{ep}.pvalue_r,[],2); % replay pval = 1-p_shuf
    
    % load MUA firing rates
    filename = [animaldir,'MUA', '0',num2str(day),'-0',num2str(ep),'.mat'];
    load(filename)
    MUA_FR = [MUA.time',MUA.data(:,1)];
    
    %% restrict ripple to NREM
    if ~isempty(SWSlist)
        pupil_trace_id = find(pupil_trace(:,1)>= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
        pupil_trace = pupil_trace(pupil_trace_id,:);
        
        validid = find(candidate_times(:,1) >= pupil_trace(1,1) & candidate_times(:,2) <= pupil_trace(end,1));  % get time stamps
        candidate_ep = candidate_times(validid,1:2);
        candidate_pval = replay_pval(validid);
        
        candidate_MUA = [];
        candidate_list = [];
        candidate_pvals = [];
        for rip = 1:length(validid)
            ripid = find(MUA_FR(:,1) >= candidate_ep(rip,1) & MUA_FR(:,1) <= candidate_ep(rip,2));
            if ~isempty(ripid)
                candidate_MUA = [candidate_MUA;nanmean(MUA_FR(ripid,2))];
                candidate_list = [candidate_list;candidate_ep(rip,:)];
                candidate_pvals = [candidate_pvals;candidate_pval(rip)];
            end
        end
        
        % restrict ripples to NREM
        riplist_ep = [];
        riplist_MUA_ep = [];
        riplist_Pval_ep = [];
        for swsseg = 1:length(SWSlist(:,1))
            validid = find(candidate_list(:,1) >= SWSlist(swsseg,1) & candidate_list(:,2) <= SWSlist(swsseg,2));  % get time stamps
            riplist_ep = [riplist_ep;candidate_list(validid,:)];
            riplist_Pval_ep = [riplist_Pval_ep;candidate_pvals(validid)];
            riplist_MUA_ep = [riplist_MUA_ep;candidate_MUA(validid)];
        end
        
        % zscore MUA FR
        riplist_PvalZ_ep = riplist_Pval_ep; % replay pval is not zscored 
        riplist_MUAZ_ep = (riplist_MUA_ep - nanmean(riplist_MUA_ep))./nanstd(riplist_MUA_ep);
        %%
        for i = 1:length(riplist_ep)
            rippletime =  riplist_ep(i,:); 
            tempid  = find(pupil_trace(:,1) >= rippletime(1)-40 & pupil_trace(:,1) <= rippletime(1)+40);
            pupil_d = pupil_trace(tempid,2);
            pupil_time = pupil_trace(tempid,1)-rippletime(1);
            if length(find(~isnan(pupil_d))) > 10
                pupil_d_interp = interp1(pupil_time(~isnan(pupil_d)),pupil_d(~isnan(pupil_d)),win);
                ripple_trigger_pupil = [ripple_trigger_pupil;pupil_d_interp];
                replayPval_trigger_pupil = [replayPval_trigger_pupil;riplist_PvalZ_ep(i)];
                MUA_trigger_pupil = [MUA_trigger_pupil;riplist_MUAZ_ep(i)];
            end
        end
    end
end
%% normalized to baseline
replayPval_triggerZ_pupil = (replayPval_trigger_pupil - nanmean(replayPval_trigger_pupil))./nanstd(replayPval_trigger_pupil);
baseline_mean = nanmean(ripple_trigger_pupil(:,1:300),2);
baseline_std = nanstd(ripple_trigger_pupil(:,1:300)')';
ripple_triggerZ_pupil = (ripple_trigger_pupil - repmat(baseline_mean,1,length(win)))./repmat(baseline_std,1,length(win));
ripple_triggerZ_pupil_sum = sum(ripple_triggerZ_pupil,2);
validid = find(~isnan(ripple_triggerZ_pupil_sum));

ripple_triggerZ_pupil = ripple_triggerZ_pupil(validid,:);
replayPval_triggerZ_pupil = replayPval_triggerZ_pupil(validid);
MUA_trigger_pupil = MUA_trigger_pupil(validid);
ripple_triggerZ_pupil_mean = nanmean(ripple_triggerZ_pupil(:,301:end),2);
%% sorted replay by pupil size
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
    replayPval_tile{i} = replayPval_triggerZ_pupil(currentid);
    MUA_tile{i} = MUA_trigger_pupil(currentid);
end
%% resample to match the MUA during replay
samplenum = 1320; % randomly take half of the data
MUA_trigger_pupil_template = MUA_tile{6};
% downsample for matching
randid = randperm(length(MUA_trigger_pupil_template));
MUA_trigger_pupil_template = MUA_trigger_pupil_template(randid(1:samplenum));
for i = 1:5
    [indices, ~] = resample_distribution(MUA_trigger_pupil_template,MUA_tile{i});
    replayPval_tile_resample{i} = replayPval_tile{i}(indices);
end
replayPval_tile_resample{6} = replayPval_tile{6}(randid(1:samplenum));
