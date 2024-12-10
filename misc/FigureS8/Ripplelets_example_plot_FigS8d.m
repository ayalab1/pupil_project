clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
filedir = '/Volumes/Extreme Pro/Pupil_data/';

animalname = 'PPP4';
day = 8;
ep = 3;
eventnum = [13186,13189]; % long trace

CA1ripCH = 59; % ripple channel
CA1SwCH = 123; % Sharp wave channel 

daystring = num2str(day);
dir = [filedir,animalname,'/day',daystring,'/'];
animalprefix = ['day',daystring];
animalname = [animalname,'-',animalprefix];

% get ripplelets time
load(sprintf('%s/%s.ripplelets.events.mat', dir, animalprefix));% rippleHSE info
single_eventtime = [ripplelets.timestamps(eventnum(1),1),ripplelets.timestamps(eventnum(end),2)];

load(fullfile(dir,[animalprefix,'.animal.behavior.mat'])) % behavior
pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]);
%% pupil info
load(fullfile(dir,[animalprefix,'.MergePoints.events.mat'])) % Session info
epochtimes = MergePoints.timestamps(ep,:); % epoch start and end time
%% get ripple LFPs around the event
%-----load eeg using CA1 ripple channel-----%
swlfp_all = getLFP([CA1SwCH], 'basepath', dir,'basename', animalprefix); 
lfp_all = getLFP([CA1ripCH], 'basepath', dir,'basename', animalprefix); 
%%
eegtid = find(lfp_all.timestamps >= single_eventtime(1)-0.8  & lfp_all.timestamps <= single_eventtime(2)+0.1);
eegtimes_single = lfp_all.timestamps(eegtid);
eeg_single = lfp_all.data(eegtid);
sweeg_single = swlfp_all.data(eegtid);

figure,plot(eegtimes_single- single_eventtime(1),eeg_single)
hold on
plot(eegtimes_single- single_eventtime(1),sweeg_single-2000)
