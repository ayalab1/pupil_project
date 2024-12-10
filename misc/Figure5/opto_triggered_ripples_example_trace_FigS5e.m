clc
clear all
close all
%% add code to path
addpath(genpath('/Users/wenbotang/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';
animalprefix = 'PPP15';
day = 2;

daystring = num2str(day);
prefix = ['day',daystring];
dir = [dir,animalprefix,'/day',daystring,'/'];
cd(dir)
%% define parameters
plotRipples = 0;
offset_time = 0;
RipCH = 45; % PPP15
example_timerange = [1509,1549];
example_pulse = [1512.8,1513.8;...
                1532,1533;...
                1547,1548];
%% load ripples and stim pulses
% get ripple LFP
Ripple_lfp = getLFP([RipCH], 'basepath', dir,'basename', prefix); 

par        = LoadXml( [dir, prefix, '.xml'] );
% pull parameters from .xml file
SR         = par.lfpSampleRate; % lfp sampling rate

% filter to ripple band, see also Sebastian et al., Nat Neu 2023
Ripple_lfp_flit = bz_Filter(Ripple_lfp.data,'passband',[80 inf],'filter','fir1','nyquist',SR/2);
%% load pupil info 
load(fullfile(dir,[prefix,'.animal.behavior.mat'])) % behavior file
pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % time, normalized diameter
% smooth out jumping point
pupil_d = pupil_trace(:,2);
pupil_time = pupil_trace(:,1);
pupil_d_interp = interp1(pupil_time(~isnan(pupil_d)),pupil_d(~isnan(pupil_d)),pupil_time);
pupil_SR = 1/nanmedian(diff(pupil_trace(:,1)));
pupil_d_fit = locsmooth(pupil_d_interp,pupil_SR,5);
%% pick out example traces
validid = find(pupil_time >= example_timerange(1) & pupil_time <= example_timerange(2));
pupil_d_ex = [pupil_time(validid),pupil_d_fit(validid)'];
%% example LFP segments
for n = 1:3
    validid = find(Ripple_lfp.timestamps>= example_pulse(n,1) & Ripple_lfp.timestamps <= example_pulse(n,2));
    pulse_lfp{n} = [Ripple_lfp.timestamps(validid),double(Ripple_lfp.data(validid)),double(Ripple_lfp_flit(validid))];
end
%% load actual pulse time
load(fullfile(dir,[prefix,'.pulses.events.mat'])) % opto pulse

pulses_dur = pulses.duration;
validid = find(pulses_dur > 0.09);
pulselist = pulses.timestamps(validid,:) + offset_time;
%% spectrograms
defaultGraphicsSetttings
for n = 1:3
    validid = find(Ripple_lfp.timestamps>= example_pulse(n,1) & Ripple_lfp.timestamps <= example_pulse(n,2));
    temp_time = Ripple_lfp.timestamps(validid);
    wavespec = WaveSpec(Ripple_lfp,'intervals',example_pulse(n,:),'frange',[50 250],...
        'roundfreqs',false,'nfreqs',100,'ncyc',5,...
        'fvector',[],'space','log');
    figure
    contourf(temp_time,wavespec.freqs,abs(wavespec.data'),30,'LineColor','none');
    hold on;
    set(gca,'YScale','log');
    ylim([wavespec.freqs(1) wavespec.freqs(end)]);
    colormap parula;
    xlabel('time (ms)'); ylabel('frequency (Hz)');
    clim([0,4000])
    pulse_spec{n} = wavespec;
end




