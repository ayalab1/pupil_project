clc
clear all
close all
%%
% directory where replay decoding results are saved
replaydir = '/Users/wenbotang/Documents/Remote_HMM_files/AYALab_Code/wb_scripts/Hongyu_scripts/Replay_decoding_20230526/';
% add codes to path
addpath('/Users/wenbotang/Documents/MATLAB/colormapformulae/')
addpath(genpath('/Users/wenbotang/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
filedir = '/Volumes/Extreme Pro/Pupil_data/';

% example events in Fig2a are PPP7 day8 epoch5 event 427 (decoded traj 1)
% and 430 (decoded traj 2)

animalname = 'PPP7';
day = 8;
ep = 5;
epRUN = 4;
exampleevent_idx = 430; % 430 or 427
decoded_traj = 2; % 1 or 2; T-maze has 2 different trajectories
track = [1,2];%possible trajectories
CA1ripCH = 17; % ripple channel
CA1SwCH = 27; % Sharp wave channel 

daystring = num2str(day);
dir = [filedir,animalname,'/day',daystring,'/'];
animalprefix = ['day',daystring];
animalname = [animalname,'-',animalprefix];
%% parameters used by replay decoding
tBinSz = 10; % step size; default temporal bin in ms; O'Neill 2017 and Farooq, 2019
tBinSz_sm = 20; % window size; in ms
%% load data
load(fullfile(dir,[animalprefix,'.MergePoints.events.mat'])) % Session info
epochtimes = MergePoints.timestamps(ep,:); % epoch start and end time
load(fullfile(dir,[animalprefix,'.animal.behavior.mat'])) % behavior

% get ripple time
load(sprintf('%s/%srippleHSEs.events.mat', dir, animalprefix));% rippleHSE info
riptimes = [rippleHSEs{day}{ep}.starttime',rippleHSEs{day}{ep}.endtime'];  

dur = 1000*(riptimes(:,2) - riptimes(:,1));
keepidx = find(dur >= 50);%at least 5 bins, 50 ms for 10ms bins
riptimes = riptimes(keepidx,:);

% get replay decoding info
load(sprintf('%s%sreplayHSEseqdecode_CA1_%02d.mat', replaydir,animalname,ep), 'replaytrajactory');
eventidx = replaytrajactory{ep}.eventidx;
eventtimes = riptimes(eventidx,:);

% get event info
exampleevent_times = eventtimes(exampleevent_idx,:);

pMat = replaytrajactory{ep}.pMat{exampleevent_idx}{decoded_traj}.pMat;
bstline = replaytrajactory{ep}.pMat{exampleevent_idx}{decoded_traj}.bstline;
bstline_time = replaytrajactory{ep}.pMat{exampleevent_idx}{decoded_traj}.timevec;

pvalues = replaytrajactory{ep}.pvalue(exampleevent_idx,track);
rvalues = replaytrajactory{ep}.rvalues(exampleevent_idx,track);
pvalues_r = replaytrajactory{ep}.pvalue_r(exampleevent_idx,track);
%% get the fitted line
eventtime(1) = exampleevent_times(1);
eventtime(2) = exampleevent_times(2);

timevec = exampleevent_times(1):tBinSz/1000:exampleevent_times(2);

timevec_raw = bstline_time-bstline_time(1);
timevec_raw = timevec_raw./max(timevec_raw);
timevec_interp = timevec - timevec(1)+tBinSz_sm/2000;
timevec_interp = timevec_interp./max(timevec_interp);
bestline_interp = interp1(timevec_raw,bstline,timevec_interp,'PCHIP');
%% plot posterior prob matrix
figure('Position',[100,100,170,170]),
imagesc(pMat)
colormap(colormapformulae([34,35,36])); % it is same as 'black-red-yellow-white'
caxis([0,0.3]);
hold on
plot(1:length(timevec),bestline_interp,'color','cyan','linewidth',2)
colorbar
%% get ripple LFPs around the event
%-----load eeg using CA1 ripple channel-----%
lfp_all = getLFP([CA1ripCH], 'basepath', dir,'basename', animalprefix); 

eegtid = find(lfp_all.timestamps >= eventtime(1)-0.2  & lfp_all.timestamps <= eventtime(2)+0.2);
eegtimes = lfp_all.timestamps(eegtid);
eeg_rip = lfp_all.data(eegtid);
figure,plot(eegtimes- eventtime(1),eeg_rip)
%% get Sharp-wave LFPs around the event
%-----load eeg using CA1 sw channel-----%
lfp_all = getLFP([CA1SwCH], 'basepath', dir,'basename', animalprefix); 

eegtid = find(lfp_all.timestamps >= eventtime(1)-0.2  & lfp_all.timestamps <= eventtime(2)+0.2 );
eegtimes = lfp_all.timestamps(eegtid);
eeg = lfp_all.data(eegtid);
figure,plot(eegtimes- eventtime(1),eeg)
%% get MUA firing around the event
%-----load MUA activity-----%
load(sprintf('%sMUA%02d-%02d.mat', dir,day, ep));
MUAtid = find(MUA.time >= eventtime(1)-0.2  & MUA.time <= eventtime(2)+0.2);
MUAtimes = MUA.time(MUAtid);
MUAseg = MUA.data(MUAtid,1);
figure,plot(MUAtimes- eventtime(1),MUAseg)