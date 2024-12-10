% this is the main script to generate Fig. 5f (opto generated ripples)
clc
clear all
close all
%% define data
dir = '/Volumes/Pupil_data/';

% data list, [{animal},{day},{epoch}]
animal_info = [{'PPP10'},{23},{1};...
    {'PPP10'},{26},{1};...
    {'PPP10'},{27},{1};...
    {'PPP11'},{1},{1};...
    {'PPP13'},{9},{5};...
    {'PPP13'},{12},{4};...
    {'PPP14'},{1},{1};...
    {'PPP12'},{1},{1};...
    {'PPP15'},{2},{1}];
%% set parameters
win = -30:0.1:30; % window around ripple, in sec
%% gather results
ripple_trigger_pupil = []; % reset
rippleamp_trigger_pupil = [];

for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    ep = animal_info{session_list,3};
    daystring = num2str(day);
    
    if  strcmp(animalname,'PPP13') && day == 1
        prefix = ['day',daystring,'_rippleg'];
        animaldir = [dir,animalname,'/day',daystring,'/',prefix,'/'];
    elseif strcmp(animalname,'PPP14') && day == 1
        prefix = ['day',daystring,'_ripple_g'];
        animaldir = [dir,animalname,'/day',daystring,'/',prefix,'/'];
    elseif strcmp(animalname,'PPP10') && day == 27
        prefix = ['day',daystring,'_g'];
        animaldir = [dir,animalname,'/day',daystring,'/',prefix,'/'];
    elseif strcmp(animalname,'PPP12') && day == 1
        prefix = 'ripple_g';
        animaldir = [dir,animalname,'/day',daystring,'/',prefix,'/'];
    else
        prefix = ['day',daystring];
        animaldir = [dir,animalname,'/day',daystring,'/'];
    end
    
    disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])

    % load data
    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % [time, normalized pupil diameter]
    epochtimes = [behavior.epochs{ep}.startTime,behavior.epochs{ep}.stopTime];  

    % ripple
    load(fullfile(animaldir,[prefix,'.optoripples.events.mat'])) % optogenetically generated ripples

    if ~isempty(freq_threshold)
        rip_freq = ripples.frequency;
        validid = find(rip_freq > freq_threshold); % remove some pulses that failed to generate events at ripple frequency
        riplist_all = [ripples.timestamps(validid,1),ripples.timestamps(validid,1) + ripples.duration(validid),ripples.amplitude(validid)]; % [start, end, amplitude]
    else
        riplist_all = [ripples.timestamps(:,1),ripples.timestamps(:,1) + ripples.duration,ripples.amplitude]; % [start, end, amplitude]
    end

    validid = find(riplist_all(:,1) >= epochtimes(1) & riplist_all(:,2) <= epochtimes(2));  % get time stamps
    riplist = riplist_all(validid,1:2);
    riplist_amp = riplist_all(validid,3);
    
    % zscore ripple amplitude
    riplist_ampZ_ep = (riplist_amp - nanmean(riplist_amp))./nanstd(riplist_amp);
    
    pupil_trace_id = find(pupil_trace(:,1)>= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
    pupil_trace = pupil_trace(pupil_trace_id,:);
    
    for i = 1:length(riplist)
        rippletime =  riplist(i,:); 
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
%% normalization
baseline_mean = nanmean(ripple_trigger_pupil(:,1:300),2);
baseline_std = nanstd(ripple_trigger_pupil(:,1:300)')';
ripple_triggerZ_pupil = (ripple_trigger_pupil - repmat(baseline_mean,1,length(win)))./repmat(baseline_std,1,length(win));
ripple_triggerZ_pupil_sum = sum(ripple_triggerZ_pupil,2);
validid = find(~isnan(ripple_triggerZ_pupil_sum));

ripple_triggerZ_pupil = ripple_triggerZ_pupil(validid,:);
rippleamp_triggerZ_pupil = rippleamp_trigger_pupil(validid);
ripple_triggerZ_pupil_mean = nanmean(ripple_triggerZ_pupil(:,301:end),2);
%% correlation
[rval,pval] = corr(rippleamp_triggerZ_pupil,ripple_triggerZ_pupil_mean,'type','spearman');
%% sorted pupil by ripple amplitude
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