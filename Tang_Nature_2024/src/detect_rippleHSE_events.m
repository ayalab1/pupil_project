%----------------------------------------------------------%
%  This is the main script for detecting SWR-associated    %     
%  high-synchrony events (HSEs)                            %
%  -- Wenbo Tang (Sep 13, 2023)                            %
%----------------------------------------------------------%
clc
close all
clear all;
%% define data
dir = '/Volumes/Pupil_data/';

animalprefix = 'PPP8'; % animal prefix
day = 7; % day
daystring = num2str(day);
epochs = 1:5; % epochs
dir = [dir,animalprefix,'/day',daystring,'/'];
prefix = ['day',daystring];
%% define parameters
savedata = 1;
binsize = 1; % ms
min_overlap_duration = 20; %ms
%%
load(fullfile(dir,[prefix,'.MergePoints.events.mat'])) % Session info

% ripple info
if strcmp(prefix,'day28') && strcmp(animalprefix,'PPP4')
    load(fullfile(dir,[prefix,'.ripples_task.events.mat'])) % ripple
elseif strcmp(prefix,'day8') && strcmp(animalprefix,'PPP7')
    load(fullfile(dir,[prefix,'.dorsalripples.events.mat'])) % ripple
else
    load(fullfile(dir,[prefix,'.ripples.events.mat'])) % ripple
end

load(sprintf('%s/%sHSEs.events.mat', dir, prefix));% HSE info
%% epoch loop
for ep = epochs
    event_count = 0;
    epochtimes = MergePoints.timestamps(ep,:); % epoch start and end time
    validid = find(ripples.timestamps(:,1) >= epochtimes(1) & ripples.timestamps(:,1) <= epochtimes(2));
    riptimes = ripples.timestamps(validid,:);
    timevec = epochtimes(1):binsize/1000:epochtimes(2);
    ripvec = list2vec(riptimes,timevec);
    ripvec = double(ripvec);
    HSEtimes = [HSEs{day}{ep}.starttime,HSEs{day}{ep}.endtime];
    % check if a HSE is associated with a SWR
    for event = 1:length(HSEtimes(:,1))
        HSE_event = HSEtimes(event,:);
        evt_id = find(timevec >= HSE_event(1) & timevec <= HSE_event(2));
        overlap_dur = sum(ripvec(evt_id));
        if overlap_dur >= min_overlap_duration
            event_count = event_count +1;
            ripHSE.starttime(event_count) = HSEs{day}{ep}.starttime(event);
            ripHSE.endtime(event_count) = HSEs{day}{ep}.endtime(event);
            ripHSE.activecellnum(event_count) = HSEs{day}{ep}.activecellnum(event);
        end
    end
    
    if event_count == 0
        ripHSE.starttime = [];
        ripHSE.endtime = [];
        ripHSE.activecellnum = [];
    end
    
    ripHSE.samprate = HSEs{day}{ep}.samprate;
    ripHSE.threshold = HSEs{day}{ep}.threshold;
    ripHSE.baseline = HSEs{day}{ep}.baseline;
    ripHSE.std = HSEs{day}{ep}.std;
    ripHSE.minimum_duration = HSEs{day}{ep}.minimum_duration;
    ripHSE.minimum_overlap = min_overlap_duration;
    ripHSE.cellthresh = HSEs{day}{ep}.cellthresh;
    
    rippleHSEs{day}{ep} = ripHSE;
    clear ripHSE;
end
%% save data
if savedata
    save(sprintf('%s/%srippleHSEs.events.mat', dir, prefix), 'rippleHSEs');
end


            
            
        
        

