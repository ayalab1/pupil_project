% this is the main script to generate Fig. 1c, left
clc
clear all
close all
%% define data
% example data comes from PPP4, day 8, epoch 3
dir = '/Volumes/Pupil_data/';
animalname = 'PPP4'; 
day = 8;
ep = 3; %epoch
daystring = num2str(day);
animaldir = [dir,animalname,'/day',daystring,'/'];
prefix = ['day',daystring];
%% define parameters
mindur = 600; % in s
winsize = 300; % in s, for acg
%% load data
load(fullfile(animaldir,[prefix,'.MergePoints.events.mat'])) % Session info
epochtimes = MergePoints.timestamps(ep,:);

% pupil info
load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % [time, normalized pupil diameter]
valididx = find(pupil_trace(:,1) >= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
pupil_trace = pupil_trace(valididx,:);

% smooth out jumping point
pupil_d = pupil_trace(:,2);
pupil_time = pupil_trace(:,1);
pupil_d_interp = interp1(pupil_time(~isnan(pupil_d)),pupil_d(~isnan(pupil_d)),pupil_time);
pupil_SR = 1/nanmedian(diff(pupil_trace(:,1)));
pupil_d_fit = locsmooth(pupil_d_interp,pupil_SR,2); % from chronux toolbox
%%  load sleep states (1 WAKE, 3 NREM, 5 REM)
load(fullfile(animaldir,[prefix,'.SleepState.states.mat'])) % sleep states
Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM
sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
Sleep_state_ep = Sleep_state(sleepss_id,:);
SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
SWSid = find(Sleep_state_ep(:,2) == 3); % restrict to NREM
SWS_vec_ep(SWSid) = 1;
SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));
SWSdur = SWSlist(:,2) - SWSlist(:,1);
SWSlist = SWSlist(find(SWSdur > mindur),:); % last for a minimal duration
%% ACG and power spectrum of ACG for each NREM episode
for SWSseg = 1:length(SWSlist(:,1))
    current_SWS = SWSlist(SWSseg,:);
    tempid = find(pupil_time >= current_SWS(1) & pupil_time <= current_SWS(2));
    pupil_d_SWS = pupil_d_fit(tempid);

    [pupil_ACG,lags] = autocorr(pupil_d_SWS,'NumLags',round(winsize*pupil_SR)); %ACG

    lags = lags/pupil_SR;
    lags = [-lags(end:-1:2),lags];
    pupil_ACG = [pupil_ACG(:,end:-1:2),pupil_ACG];
    pupil_ACG_corrected_tmp(SWSseg,:) = pupil_ACG./(1-abs(lags/pupil_SR)/(SWSlist(SWSseg,2)-SWSlist(SWSseg,1)));% corrected, Mizuseki et al. (2009)

    % spectrogram
    pupil_acg.timestamps = lags';
    pupil_acg.data = pupil_ACG_corrected_tmp(SWSseg,:)';
    pupil_acg.samplingRate = pupil_SR;
    wavespec = WaveSpec(pupil_acg,'frange',[0.001 1],...
        'roundfreqs',false,'nfreqs',100,'ncyc',5,...
        'fvector',[],'space','log');

    wavespec_amp = abs(wavespec.data);
    avgSpect_tmp =  nanmean(wavespec_amp);
    avgSpect(SWSseg,:) = avgSpect_tmp;
end
%% plot results
pupil_ACG_corrected_tmp_full = pupil_ACG_corrected_tmp;
figure,
subplot(211)
plot(lags,pupil_ACG_corrected_tmp_full) % ACG
subplot(212)
semilogx(wavespec.freqs,avgSpect) % spectrum

