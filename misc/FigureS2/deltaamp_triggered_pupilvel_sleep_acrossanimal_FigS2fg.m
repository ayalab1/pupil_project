clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
%% list of data
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
    {'PPP8'},{7},{1};...
    {'PPP8'},{7},{3};...
    {'PPP8'},{7},{5};...
    {'PPP8'},{8},{1};...
    {'PPP8'},{8},{3}];
%% set parameters
velfilter = 1; 
win = -30:0.1:30;
%% gather data
delta_trigger_pupil = [];
deltaamp_trigger_pupil = [];
deltadur_trigger_pupil = [];

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
    epochtimes = MergePoints.timestamps(ep,:);

    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % time, normalized diameter

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

    load(fullfile(animaldir,[prefix,'.deltaWaves.events.mat'])) % delta wave
    deltalist_all = [deltaWaves.timestamps,deltaWaves.peakNormedPower - deltaWaves.troughValue]; % start, end, amplitude
    %% restrict to NREM
    if ~isempty(SWSlist)
        pupil_trace_id = find(pupil_trace(:,1)>= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
        pupil_trace = pupil_trace(pupil_trace_id,:);

        validid = find(deltalist_all(:,1) >= pupil_trace(1,1) & deltalist_all(:,2) <= pupil_trace(end,1));  % get time stamps
        deltalist = deltalist_all(validid,1:2);
        deltalist_amp = deltalist_all(validid,3);
        
        % confine delta waves to sws
        deltalist_ep = [];
        deltalist_amp_ep = [];       

        for swsseg = 1:length(SWSlist(:,1))
            validid = find(deltalist(:,1) >= SWSlist(swsseg,1) & deltalist(:,2) <= SWSlist(swsseg,2));  % get time stamps
            deltalist_ep = [deltalist_ep;deltalist(validid,:)];
            deltalist_amp_ep = [deltalist_amp_ep;deltalist_amp(validid)];
        end
        
        % zscore delta amplitude
        deltalist_ampZ_ep = (deltalist_amp_ep - nanmean(deltalist_amp_ep))./nanstd(deltalist_amp_ep);
        %%
        for i = 1:length(deltalist_ep)
            deltatime =  deltalist_ep(i,:); 
            tempid  = find(pupil_trace(:,1) >= deltatime(1)-40 & pupil_trace(:,1) <= deltatime(1)+40);
            pupil_d = pupil_trace(tempid,2);
            pupil_time = pupil_trace(tempid,1)-deltatime(1);
            if length(find(~isnan(pupil_d))) > 10
                pupil_d_interp = interp1(pupil_time(~isnan(pupil_d)),pupil_d(~isnan(pupil_d)),win);
                delta_trigger_pupil = [delta_trigger_pupil;pupil_d_interp];
                deltaamp_trigger_pupil = [deltaamp_trigger_pupil;deltalist_ampZ_ep(i)];
            end
        end
    end
end
%%
baseline_mean = nanmean(delta_trigger_pupil(:,1:300),2);
baseline_std = nanstd(delta_trigger_pupil(:,1:300)')';
delta_triggerZ_pupil = (delta_trigger_pupil - repmat(baseline_mean,1,length(win)))./repmat(baseline_std,1,length(win));
delta_triggerZ_pupil_sum = sum(delta_triggerZ_pupil,2);
validid = find(~isnan(delta_triggerZ_pupil_sum));

delta_triggerZ_pupil = delta_triggerZ_pupil(validid,:);
deltaamp_triggerZ_pupil = deltaamp_trigger_pupil(validid);
delta_triggerZ_pupil_mean = nanmean(delta_triggerZ_pupil(:,301:end),2);
%% correlation 
r = corr(delta_triggerZ_pupil_mean,deltaamp_triggerZ_pupil,'type','spearman');
%% sorted pupil by delta amplitude
figure,
Y = prctile(deltaamp_triggerZ_pupil,1/6*100);
tile1 = [min(deltaamp_triggerZ_pupil),Y];
tile2 = [Y,prctile(deltaamp_triggerZ_pupil,2/6*100)];
tile3 = [prctile(deltaamp_triggerZ_pupil,2/6*100),prctile(deltaamp_triggerZ_pupil,3/6*100)];
tile4 = [prctile(deltaamp_triggerZ_pupil,3/6*100),prctile(deltaamp_triggerZ_pupil,4/6*100)];
tile5 = [prctile(deltaamp_triggerZ_pupil,4/6*100),prctile(deltaamp_triggerZ_pupil,5/6*100)];
tile6 = [prctile(deltaamp_triggerZ_pupil,5/6*100),prctile(deltaamp_triggerZ_pupil,6/6*100)];
sixtiles = [tile1;tile2;tile3;tile4;tile5;tile6];

for i = 1:6
    currentid = find(deltaamp_triggerZ_pupil >= sixtiles(i,1) & deltaamp_triggerZ_pupil <= sixtiles(i,2));
    delta_triggerZ_pupil_tile = delta_triggerZ_pupil(currentid,:);
    delta_trigger_pupil_stats{i} = [win',nanmean(delta_triggerZ_pupil_tile)',nanstd(delta_triggerZ_pupil_tile)',sum(double((~isnan(delta_triggerZ_pupil_tile))))'];
    plot(win',nanmean(delta_triggerZ_pupil_tile)')
    hold on
end
%% sorted delta amplitude by pupil
Y = prctile(delta_triggerZ_pupil_mean,1/6*100);
tile1 = [min(delta_triggerZ_pupil_mean),Y];
tile2 = [Y,prctile(delta_triggerZ_pupil_mean,2/6*100)];
tile3 = [prctile(delta_triggerZ_pupil_mean,2/6*100),prctile(delta_triggerZ_pupil_mean,3/6*100)];
tile4 = [prctile(delta_triggerZ_pupil_mean,3/6*100),prctile(delta_triggerZ_pupil_mean,4/6*100)];
tile5 = [prctile(delta_triggerZ_pupil_mean,4/6*100),prctile(delta_triggerZ_pupil_mean,5/6*100)];
tile6 = [prctile(delta_triggerZ_pupil_mean,5/6*100),prctile(delta_triggerZ_pupil_mean,6/6*100)];
sixtiles = [tile1;tile2;tile3;tile4;tile5;tile6];

for i = 1:6
    currentid = find(delta_triggerZ_pupil_mean >= sixtiles(i,1) & delta_triggerZ_pupil_mean <= sixtiles(i,2));
    deltaamp_tile{i} = deltaamp_triggerZ_pupil(currentid);
end
