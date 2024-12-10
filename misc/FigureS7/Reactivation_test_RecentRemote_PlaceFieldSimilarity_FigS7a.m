clc
close all
clear all;
%% add codes to path
% addpath(genpath('/Users/wenbotang/Documents/Remote_HMM_files/AYALab_Code/'))
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
% list of data, [animal, day, POST sleep, familiar-novel session pairs
% animal_info = [{'PPP7'},{12},{6},{[3,5]};...
%     {'PPP7'},{14},{5},{[2,4]};...
%     {'PPP8'},{12},{5},{[2,4]};...
%     {'PPP8'},{15},{5},{[4,2]};...
%     {'PPP8'},{14},{5},{[2,4]};...
%     {'PPP8'},{17},{5},{[4,2]};...
%     {'PPP15'},{8},{4},{[2,3]};...
%     {'PPP13'},{9},{4},{[2,3]};...
%     {'PPP13'},{12},{4},{[2,3]};...
%     {'PVR4'},{1},{4},{[2,3]}];

dir = '/Volumes/Extreme Pro/Pupil_data/';
animalprefix = 'PPP7';
day = 12;
eps_RUN = [3,5];

if  strcmp(animalprefix,'SNCG09')|| strcmp(animalprefix,'PVR4')|| strcmp(animalprefix,'PVR5')
    daystring = num2str(day);
    prefix = ['pupil',daystring];
    dir = [dir,animalprefix,'/pupil',daystring,'/'];
else
    daystring = num2str(day);
    dir = [dir,animalprefix,'/day',daystring,'/'];
    prefix = ['day',daystring];
end
%% set paramters
cell_excluded = [];
speedthresh = 2; % cm/s
%% load spikes
if strcmp(animalprefix,'PPP4') ||  strcmp(animalprefix,'PPP15') || strcmp(animalprefix,'PPP13') ||... 
        (strcmp(animalprefix,'PPP8') && day == 8) || (strcmp(animalprefix,'PPP8') && (day == 15 || day == 14 || day == 17 || day == 19))||...
        strcmp(animalprefix,'PVR4') ||  strcmp(animalprefix,'PVR5') 
    spikeStructure = importSpikes('basepath',dir,'CellType',"Pyramidal Cell");
elseif strcmp(animalprefix,'PPP7') % dorsal CA1 cells only
    if day == 18 || day == 17
        spikeStructure = importSpikes('basepath',dir,'CellType',"Pyramidal Cell");
    else
        spikeStructure = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"dCA1");
    end
else
    spikeStructure = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"CA1");
end

spikes = [];
cell_count = 0;
for i = 1:length(spikeStructure.times)
    if ~ismember(i,cell_excluded)
        cell_count = cell_count + 1;
        timestamps = spikeStructure.times{i};
        id = ones(size(timestamps))*i; % the id for that unit
        spikes = [spikes; timestamps id]; % add [timestamps id] for this unit to the matrix
    end
end
spikes = sortrows(spikes); % sort spikes accoring to their time    
%% load session info and sleep states
load(fullfile(dir,[prefix,'.animal.behavior.mat'])) % behavior
vel  = behavior.speed_smooth; % running speed
load(fullfile(dir,[prefix,'.MergePoints.events.mat'])) % Session info
%% calculate the template for recent and remote memories separately, and compute the PF similarity
load(fullfile(dir,  [prefix,'.firingMapsAvg_2D.cellinfo.mat']));% load firing rate maps
PF_corr_all = [];

for ep = eps_RUN
    PF_corr_ep = [];
    epochtimes = MergePoints.timestamps(ep,:);
    vel_epID = find(vel(:,1) >= epochtimes(1) & vel(:,1) <= epochtimes(2));
    vel_ep = vel(vel_epID,:);
    
    % caculate template activity
    RUN = vec2list(vel_ep(:,2) > speedthresh,vel_ep(:,1)); % generate [start end] list of immobile epochs
    validid = find((RUN(:,2)-RUN(:,1)) > 0.05);
    RUN = RUN(validid,:);
    taskIntervals = SplitIntervals(RUN,'pieceSize',0.05); % 50 ms bins in behavior
    [templates,~,weights,~] = ActivityTemplates(spikes,'bins',taskIntervals,'mode','ica');
    assemble_num = length(weights(1,:));

    for i = 1:assemble_num
        current_weights = weights(:,i);
        assemble_mean = nanmean(abs(current_weights));
        assemble_std = nanstd(abs(current_weights));
        threshold =  assemble_mean + 3*assemble_std; % van de Ven 2016
        memeber_IDs = find(current_weights > threshold);
        memeber_num = length(memeber_IDs);

        if memeber_num > 1
            pairind = combnk(1:memeber_num,2);
            % calculate place-field similarity
            for pair = 1:length(pairind(:,1))
                cid1 = memeber_IDs(pairind(pair,1));
                cid2 = memeber_IDs(pairind(pair,2));

                rm1 = firingMaps_2D{ep}.rateMaps{cid1};
                rm2 = firingMaps_2D{ep}.rateMaps{cid2};

                rm1(find(rm1 < 0)) = NaN;
                rm2(find(rm2 < 0)) = NaN;

                PF_corr = corr(rm1(:),rm2(:),'rows','complete');

                ep_alt = setdiff(eps_RUN,ep);
                rm1 = firingMaps_2D{ep_alt}.rateMaps{cid1};
                rm2 = firingMaps_2D{ep_alt}.rateMaps{cid2};

                rm1(find(rm1 < 0)) = NaN;
                rm2(find(rm2 < 0)) = NaN;

                PF_corr_alt = corr(rm1(:),rm2(:),'rows','complete');

                PF_corr_ep = [PF_corr_ep;PF_corr,PF_corr_alt];
            end
        end
    end

    PF_corr_all = [PF_corr_all;nanmean(PF_corr_ep)];

end