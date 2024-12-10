clc
clear all
close all
%% add code to path
addpath(genpath('/Users/wenbotang/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/Pupil_OptoBehav/';
animalprefix = 'PPP13';
day = 8;
offset_sessID = 2; % time offset for spike data, from data concatenation 
daystring = num2str(day);

OptoCtrl = 0; % Opto session = 1; Ctrl Session = 0
prefix = ['day',daystring,'_post'];
spikeprefix = ['day',daystring];

if OptoCtrl
    daystring = num2str(day);
    animaldir = [dir,'/Opto/',animalprefix,'/day',daystring,'/day',daystring,'_post','/'];
    spikedir = [dir,'/Opto/',animalprefix,'/day',daystring,'/','all','/'];

else
    daystring = num2str(day);
    animaldir = [dir,'/Ctrl/',animalprefix,'/day',daystring,'/day',daystring,'_post','/'];
    spikedir = [dir,'/Ctrl/',animalprefix,'/day',daystring,'/','all','/'];
end
cd(animaldir)
%% load data
load(fullfile(animaldir,[prefix,'.pulses.events.mat'])) % opto pulse
load([animaldir,prefix,'.animal.behavior.mat']) % behavior file
pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % time, normalized diameter

load(fullfile(animaldir,[prefix,'.ripples.events.mat'])) % ripple
riplist = [ripples.timestamps, ripples.timestamps + ripples.duration];
%% pupil trace
for i = 1:length(behavior.onlinetracking.CurrentTime)
    if isempty(behavior.onlinetracking.Area{i})
        pupil_area(i,1) = behavior.onlinetracking.CurrentTime(i);
        pupil_area(i,2) = nan(1);
    else
        pupil_area(i,:) = [behavior.onlinetracking.CurrentTime(i),str2num( behavior.onlinetracking.Area{i})];
    end
end
pupil_d_fit = locsmooth(pupil_area(:,2),30,6);
%%  get ripple LFP
% ripple channels
if strcmp(animalprefix,'PPP10')
    RipCH = 9; % PPP10
elseif strcmp(animalprefix,'PPP11')
    RipCH = 26; % PPP11
elseif strcmp(animalprefix,'PPP12')
    RipCH = 45; % PPP12
elseif strcmp(animalprefix,'PPP13')
    RipCH = 45; % PPP13
elseif strcmp(animalprefix,'PPP14')
    RipCH = 33; % PPP14
elseif strcmp(animalprefix,'PPP15')
    RipCH = 45; % PPP15
elseif strcmp(animalprefix,'PPP16')
    RipCH = 3; % PPP16
end

% load LFP
Ripple_lfp = getLFP([RipCH], 'basepath', animaldir,'basename', prefix); 

par        = LoadXml( [animaldir, prefix, '.xml'] );
% pull parameters from .xml file
SR         = par.lfpSampleRate; % lfp sampling rate
% filter to ripple band, see also Sebastian et al., Nat Neu 2023
Ripple_lfp_flit = bz_Filter(Ripple_lfp.data,'passband',[80 250],'filter','fir1','nyquist',SR/2);
%% plot the example trace
figure,
subplot(311)
plot(pupil_trace(:,1),pupil_trace(:,2))
hold on
plot(pupil_trace(1:end-1,1),pupil_d_fit)
plot(pupil_trace(1:end-1,1),behavior.onlinetracking.pupil_sws)
xlim([2292,2352]);

subplot(312)
plot(Ripple_lfp.timestamps,Ripple_lfp_flit)
xlim([2292,2352]);
ylim([-4000,4000])

subplot(313)
plot(riplist(:,1),ones(size(riplist(:,1))),'o')
xlim([2292,2352]);


subplot(313)
hold on
plot(pulses.timestamps(:,1),ones(size(pulses.timestamps(:,1))),'o')
xlim([2292,2352]);

