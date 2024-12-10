clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wenbotang/Documents/Remote_HMM_files/AYALab_Code/'))
%% setting for plotting
% defaultGraphicsSetttings
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';
animalname = 'PPP7';
day = 8;
epoch = 3;
daystring = num2str(day);
dir = [dir,animalname,'/day',daystring,'/'];
animalprefix = ['day',daystring];
EI_UIDs = [258,83]; % PYR, INT UID
%% set parameters (from cell excplorer ce_MonoSynConvClick parameters)
binSize = 4/10000; 
binSize_ccg = 0.001; 
duration = 0.12;
%% load data
spikes = GetAyaSpikes(dir); % all spikes

% separate INT and PYR spikes
if  strcmp(animalname,'PPP7')
    spikes_pyr = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"dCA1");
    spikes_int = importSpikes('basepath',dir,'CellType',["Narrow Interneuron", "Wide Interneuron"],'brainRegion',"dCA1");
elseif strcmp(animalname,'PPP4') || (strcmp(animalname,'PPP8') && day == 8) || (strcmp(animalname,'PPP8') && day == 15)
    spikes_pyr = importSpikes('basepath',dir,'CellType',"Pyramidal Cell");
    spikes_int = importSpikes('basepath',dir,'CellType',["Narrow Interneuron", "Wide Interneuron"]);
else
    spikes_pyr = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"CA1");
    spikes_int = importSpikes('basepath',dir,'CellType',["Narrow Interneuron", "Wide Interneuron"],'brainRegion',"CA1");
end

load(fullfile(dir,[animalprefix,'.MergePoints.events.mat'])) % Session info
epochtimes = MergePoints.timestamps(epoch,:);

load(fullfile(dir,[animalprefix,'.cell_metrics.cellinfo.mat'])) % cell info

load([dir,animalprefix,'.PupilTile_NREM','EP',num2str(epoch),'.mat']) % pupil info

% sleep states
load(fullfile(dir,[animalprefix,'.SleepState.states.mat'])) % sleep states
Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  
sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
Sleep_state_ep = Sleep_state(sleepss_id,:);
%1 awake, 3 NREM, 5 REM
SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
SWSid = find(Sleep_state_ep(:,2) == 3);
SWS_vec_ep(SWSid) = 1;
SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));
%% calculate CCG
for tile = 1:3 % small, middle, large
    interval = [pupil_tile{tile*2-1};pupil_tile{tile*2}];
    restricted = Restrict(spikes,interval);
    [ccg,t] = CCG(restricted(:,1),restricted(:,2),'binSize',binSize,'duration',duration,'normtype','rate');
    ccgEI{tile} = ccg(:,EI_UIDs(1),EI_UIDs(2));
end
%% plot CCG
figure
plot(t,ccgEI{3})
hold on
plot(t,ccgEI{2})
plot(t,ccgEI{1})
plot([0,0],[0,500],'k--')
xlim([-0.05,0.05])
%% plot ACG from cell metrices
figure
acg_time =  -100:100; % time axis (in ms)
plot(acg_time,cell_metrics.acg.narrow(:,EI_UIDs(2)))
hold on
plot(acg_time,cell_metrics.acg.narrow(:,EI_UIDs(1)))
