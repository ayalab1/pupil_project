%--------------------------------------------------------------------%
%  This is the main script for detecting high-synchrony events (HSEs)%
%  -- Wenbo Tang (Sep 13, 2023)                                      %
%--------------------------------------------------------------------%
clc
close all
clear all
%% define data
dir = '/Volumes/Pupil_data/';
animalprefix = 'PPP8'; % animal prefix
day = 12; % day
daystring = num2str(day);
epochs = 1:5; % epochs
dir = [dir,animalprefix,'/day',daystring,'/'];
prefix = ['day',daystring];
%% define parameters
savedata = 1;
binsize = 1; % ms
thresh_nstd = 3; % threshold, 3 stds
min_suprathresh_duration = 0.030; %s, typically 30ms
velfilter = 2; % speed filter, cm/s
nstd = 15; % gaussian smooth window
activecellthresh = 0.15; % fraction of cells active required
%% load data
load(fullfile(dir,[prefix,'.animal.behavior.mat'])) % behavior
vel  = behavior.speed_smooth(:,2); % running speed
load(fullfile(dir,[prefix,'.MergePoints.events.mat'])) % Session info
% load spikes
if strcmp(animalprefix,'PPP4') ||  (strcmp(animalprefix,'PPP8') && day == 8)
    spikes = importSpikes('basepath',dir,'CellType',"Pyramidal Cell");
elseif strcmp(animalprefix,'PPP7')
    spikes = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"dCA1");
else
    spikes = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"CA1");
end
cellnum = length(spikes.UID);
%% epoch loop
for ep = epochs
    epochtimes = MergePoints.timestamps(ep,:); % epoch start and end time
    epoch = ep;
    totaltime = epochtimes(2)- epochtimes(1);
    timevec = epochtimes(1)*1000:binsize:epochtimes(2)*1000;
    % Get position
    % --------------------
    valid_id = find(behavior.timestamps >= epochtimes(1) & behavior.timestamps <= epochtimes(2));
    if isempty(valid_id) % sleep epochs don't have position and speed info
        sleepepoch = 1;
    else
        sleepepoch = 0;
        velocity = vel(valid_id);
        veltimes = behavior.timestamps(valid_id);
    end
    %%
    MUA_spikes_mat = [];
    % Get the spike times
    % --------------------
    cellcount = 0;
    for cell = 1:cellnum
        if ~isempty(spikes.times{cell})
            spikeu = spikes.times{cell}*1000;  % in ms;
            valid_id = find(spikeu >= epochtimes(1)*1000 & spikeu<= epochtimes(2)*1000);
            spikeu = spikeu(valid_id);
            if length(spikeu) >= 100
                cellcount = cellcount +1;
                % Get firing rate of neuron
                % -------------------------
                histspks_all = histc(spikeu,timevec);
                MUA_spikes_mat = [MUA_spikes_mat,histspks_all];
            end
        end
    end
    %%
    MUA_spikes = sum(MUA_spikes_mat,2);
    id = find(MUA_spikes_mat > 0);
    MUA_spikes_cellid = zeros(size(MUA_spikes_mat));
    MUA_spikes_cellid(id) = 1;
    MUA_spikes_cellnum = sum(MUA_spikes_cellid,2);
    
    % gaussian smooth
    g1 = gaussian(nstd, 3*nstd+1);
    MUA_spikes = smoothvect(MUA_spikes, g1);
    samprate = 1/(binsize/1000);
    
    cellthresh = activecellthresh * cellcount;    
    
    MUA.data(:,1) = MUA_spikes; % total number of spikes
    MUA.data(:,2) = MUA_spikes_cellnum; % total number of active cells
    MUA.time = timevec/1000;
    MUA.filtersamprate = 1/(binsize/1000);
    MUA.starttime = epochtimes(1);
    MUA.endtime = epochtimes(2);
    %%
    % save MUA trace
    if savedata
        filename = [dir,'MUA', '0',num2str(day),'-0',num2str(ep),'.mat'];
        save(filename,'MUA')
    end
    clear MUA
    %%
    % Detecting HSEs
    % --------------------
    % speed threshold? only for awake epochs
    if ~sleepepoch
        if ~isempty(velfilter)
            immobile = vec2list(velocity < velfilter,veltimes); % generate [start end] list of immobile epochs
        end
        [imtime,im_vec] = wb_list2vec(immobile,timevec/1000);
        im_index = find(im_vec);
        MUA_spikes_im = MUA_spikes(im_index);
    else
        MUA_spikes_im = MUA_spikes;
    end
    
    
    baseline = mean(MUA_spikes_im);
    stdev = std(MUA_spikes_im);
    thresh = baseline + thresh_nstd * stdev;
    mindur = round(min_suprathresh_duration * samprate);
    
    if (thresh > 0) && any(find(MUA_spikes_im<baseline))
        tmpHSE = extractevents(MUA_spikes, thresh, baseline, 0, mindur, 0)';
        [HSEtime,HSEvec] = wb_list2vec(tmpHSE(:,1:2)/samprate + epochtimes(1),timevec/1000);
        if ~sleepepoch
            if ~isempty(velfilter)     
                HSEvec = HSEvec & im_vec; 
            end
        end
        HSElist = vec2list(HSEvec,timevec/1000);

        tmpHSE = HSElist;
            
        HSE.starttime = tmpHSE(:,1);
        HSE.endtime = tmpHSE(:,2);
        
        % apply cell threshold and duration threshold
        valid_timeid = [];
        valid_activecellnum = [];
        for event = 1:length(tmpHSE(:,1))
            currenttime_id = find(timevec/1000 > tmpHSE(event,1)& timevec/1000 < tmpHSE(event,2));
            current_activecell = MUA_spikes_cellid(currenttime_id,:);
            current_activecell = sum(current_activecell);
            activecellnum = length(find(current_activecell > 0));
                
            current_dur = tmpHSE(event,2)-tmpHSE(event,1);
            if (activecellnum > cellthresh) && (current_dur < 0.7)
                valid_timeid = [valid_timeid;event];
                valid_activecellnum = [valid_activecellnum;activecellnum];
            end
        end
        HSE.starttime = tmpHSE(valid_timeid,1);
        HSE.endtime = tmpHSE(valid_timeid,2);
        HSE.activecellnum = valid_activecellnum;
            
     else
        HSE.starttime = [];
        HSE.endtime = [];
     end
     HSE.samprate = samprate;
     HSE.threshold = thresh;
     HSE.baseline = baseline;
     HSE.std = stdev;
     HSE.minimum_duration = min_suprathresh_duration;
     HSE.cellthresh = cellthresh;

     HSEs{day}{ep} = HSE;
     HSEduration = round(sum(HSE.endtime-HSE.starttime));
     disp(sprintf('d %d e %d: %d s of HSEs',day,ep,HSEduration))
     clear HSE;
end
%%
if savedata
    save(sprintf('%s/%sHSEs.events.mat', dir, prefix), 'HSEs');
end
