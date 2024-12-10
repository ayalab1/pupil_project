clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
addpath(genpath('/Users/wt248/Documents/MATLAB/matplotlib'))
%% setting for the plots
defaultGraphicsSetttings
%% define data
% example data comes from PPP7, day 12, epoch 4
dir = '/Volumes/Extreme Pro/Pupil_data/';
animalname = 'PPP7'; 
day = 12;
ep = 4;
daystring = num2str(day);
animaldir = [dir,animalname,'/day',daystring,'/'];
prefix = ['day',daystring];
%% define parameters
no_xbins = 20;
no_ybins = 20;
total_bins = no_xbins * no_ybins;
mindur = 10;
window = 2;
noverlap = window-1; %1s dt
smoothfact = 15;

%for NarrowbandTheta, see also sleepscoreMaster.m
f_all = [1 20]; % broad band, for PSD slope, in Hz
f_theta = [5 10];  % theta band,in Hz
f_delta = [1 4]; % delta band, in Hz
IRASA = true;
ThIRASA = true;
%% load data
load(fullfile(animaldir,[prefix,'.MergePoints.events.mat'])) % Session info
epochtimes = MergePoints.timestamps(ep,:);

% pupil info
load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % time, normalized diameter
valididx = find(pupil_trace(:,1) >= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
pupil_trace = pupil_trace(valididx,:);

% sleep state (1 WAKE, 3 NREM, 5 REM)
state_ids = [1,3,5]; % WAKE, NREM, REM
load(fullfile(animaldir,[prefix,'.SleepState.states.mat'])) % sleep states
Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM
sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
Sleep_state_ep = Sleep_state(sleepss_id,:);

for state = 1:3
    State_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
    Stateid = find(Sleep_state_ep(:,2) == state_ids(state));
    State_vec_ep(Stateid) = 1;
    Statelist = vec2list(State_vec_ep, Sleep_state_ep(:,1));
    Statedur = Statelist(:,2) - Statelist(:,1);
    Statelist = Statelist(find(Statedur > mindur),:);
    Statelist_all{state} = Statelist;
end

% LFP features, see also sleepscoreMaster.m
load(fullfile(animaldir,[prefix,'.EMGFromLFP.LFP.mat'])) % EMG from LFP
EMG_LFP = [EMGFromLFP.timestamps,EMGFromLFP.data];

load(fullfile(animaldir,[prefix,'.SleepScoreLFP.LFP.mat'])) % NarrowbandTheta and PSD slope
Th_LFP = double(SleepScoreLFP.thLFP);
Sw_LFP = double(SleepScoreLFP.swLFP);
%% calculate SW (spectrum slope)
% Put the LFP in the right structure format
lfp.data = SleepScoreLFP.swLFP;
lfp.timestamps = SleepScoreLFP.t;
lfp.samplingRate = SleepScoreLFP.sf;

% Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,'frange',[4 90],'IRASA',IRASA);
broadbandSlowWave = -specslope.data; %So NREM is higher as opposed to lower
t_clus = specslope.timestamps;
specdt = 1./specslope.samplingRate;
    
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes = find(totz>3);

% Smooth and 0-1 normalize
broadbandSlowWave(badtimes) = nan;
broadbandSlowWave = smooth(broadbandSlowWave,smoothfact./specdt);
%% calculate theta 
% Put the LFP in the right structure format
lfp.data = SleepScoreLFP.thLFP;
% Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
        'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
   
% Remove transients
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);

% theta
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);

thratio(badtimes_TH) = nan;     %Remove transients
thratio = smooth(thratio,smoothfact./specdt);    

% delta
dlfreqs = (thFFTfreqs>=f_delta(1) & thFFTfreqs<=f_delta(2));
dlratio = max((thFFTspec(dlfreqs,:)),[],1);

dlratio(badtimes_TH) = nan;     %Remove transients
dlratio = smooth(dlratio,smoothfact./specdt);
%% get ripple LFPs
if strcmp(prefix,'day28') && strcmp(animalname,'PPP4')
    load(fullfile(animaldir,[prefix,'.ripples_task.events.mat'])) % ripple
elseif strcmp(prefix,'day8') && strcmp(animalname,'PPP7')
    load(fullfile(animaldir,[prefix,'.dorsalripples.events.mat'])) % ripple
else
    load(fullfile(animaldir,[prefix,'.ripples.events.mat'])) % ripple
end
RipCH = ripples.detectorinfo.detectionparms.Channels(1);
RipBp = ripples.detectorinfo.detectionparms.ripBP;
Ripple_lfp = getLFP([RipCH], 'basepath', animaldir,'basename', prefix); 

[specslope,spec] = bz_PowerSpectrumSlope(Ripple_lfp,window,window-noverlap,...
        'nfreqs',200,'frange',[50,400],'IRASA',IRASA);
specdt = 1./specslope.samplingRate;
RipFFTspec = specslope.resid';
RipFFTspec(RipFFTspec<0)=0;

% Remove transients
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_RIP = find(totz>3);
    
ripFFTfreqs = specslope.freqs';
ripfreqs = (ripFFTfreqs>=RipBp(1) & ripFFTfreqs<=RipBp(2));
ripratio = max((RipFFTspec(ripfreqs,:)),[],1);

ripratio(badtimes_RIP) = nan;     %Remove transients
ripratio = smooth(ripratio,smoothfact./specdt); 
%% make sure that all the LFP features have the same length
EMG_LFP_interp = interp1(EMG_LFP(:,1), EMG_LFP(:,2),specslope.timestamps);
EMG_LFP = EMG_LFP_interp;
%% restrict LFP features to the current epoch
valididx = find(specslope.timestamps >= epochtimes(1) & specslope.timestamps <= epochtimes(2));
LFP_features = [specslope.timestamps(valididx),dlratio(valididx),thratio(valididx),ripratio(valididx),...
    broadbandSlowWave(valididx),EMG_LFP(valididx)]; % time, delta, theta, ripple, SW (or slope), EMG
LFP_features = double(LFP_features);

% make pupil trace the same length as the LFP
pupil_trace_interp = interp1(pupil_trace(:,1),pupil_trace(:,2),LFP_features(:,1));
pupil_trace = [LFP_features(:,1),pupil_trace_interp];
%% restrict pupil trace to each sleep state for zscoring
for state = 1:3
    pupil_trace_State = [];
    LFP_features_State = [];
    Statelist = Statelist_all{state};
    validid_gather = [];
    if ~isempty(Statelist)
        for seg = 1:length(Statelist(:,1))
            validid = find(pupil_trace(:,1) >= Statelist(seg,1) & pupil_trace(:,1) <= Statelist(seg,2));  % get time stamps
            pupil_trace_State = [pupil_trace_State;pupil_trace(validid,:)];
            validid_gather = [validid_gather;validid];
        end
    end
    % zscore pupil
    pupil_traceZ_State = (pupil_trace_State(:,2) - nanmean(pupil_trace_State(:,2)))./nanstd(pupil_trace_State(:,2));
    pupil_traceZ{state} = [pupil_trace(validid_gather,1),validid_gather,pupil_traceZ_State]; % [time, time ID, pupil size zscored]
end
%% run UMAP across all sleep states
addpath(genpath('/Users/wt248/Documents/MATLAB/umapAndEppFileExchange (4.1)'))
% when running the first time, save template:
% LFP_features_umap = run_umap(LFP_features(:,2:end),'n_components',
% 2,'min_dist',0.6,'n_neighbors',20,'metric','cosine','save_template_file', 'SleepStates_UMAP_PPP7Day12Ep4_template3.mat'); 

% use the template for replicability
LFP_features_umap = run_umap(LFP_features(:,2:end),'n_components', 2,'min_dist',0.6,'n_neighbors',20,'metric','cosine','template_file', 'SleepStates_UMAP_PPP7Day12Ep4_template3.mat'); 
%% color each state with its own zscored pupil labels
state_indices_gather =[];
state_label_gather = [];
pupil_label_gather = [];
for state = 1:3
    state_indices_gather = [state_indices_gather;pupil_traceZ{state}(:,2)];
    state_label_gather = [state_label_gather;ones(size(pupil_traceZ{state}(:,2)))*state];
    pupil_label_gather= [pupil_label_gather;pupil_traceZ{state}(:,3)];
end
state_indices_gather = state_indices_gather(end:-1:1);
pupil_label_gather = pupil_label_gather(end:-1:1);
state_label_gather = state_label_gather(end:-1:1);
%% rearrange state ids for visualization purpose
state_label_gather_flip = state_label_gather;
id1 = find(state_label_gather == 1);
id2 = find(state_label_gather == 2);
id3 = find(state_label_gather == 3);
state_label_gather_flip(id1) = 3;
state_label_gather_flip(id2) = 2;
state_label_gather_flip(id3) = 1;
%% plot results
colormap_states = plasma(100);
colormap_rip = viridis(100);
figure,
scatter(LFP_features_umap(state_indices_gather,1),LFP_features_umap(state_indices_gather,2),30*ones(size(state_label_gather)),state_label_gather_flip,'filled',...
    'MarkerEdgeAlpha',0.7,'MarkerFaceAlpha',0.7)
colormap(colormap_states)

figure,
scatter(LFP_features_umap(state_indices_gather,1),LFP_features_umap(state_indices_gather,2),30*ones(size(state_label_gather)),pupil_trace(state_indices_gather,2),'filled',...
    'MarkerEdgeAlpha',0.7,'MarkerFaceAlpha',0.7)
colormap(colormap_rip)
%%
figure,
color_gray = [134, 138, 145 ] ./ 255;
NREM_indices = pupil_traceZ{2}(end:-1:1,2);
pupil_label_NREM = pupil_traceZ{2}(end:-1:1,3);

scatter(LFP_features_umap(state_indices_gather,1),LFP_features_umap(state_indices_gather,2),20,color_gray,'filled')
hold on
scatter(LFP_features_umap(NREM_indices,1),LFP_features_umap(NREM_indices,2),30*ones(size(pupil_label_NREM)),pupil_label_NREM,'filled',...
    'MarkerEdgeAlpha',0.7,'MarkerFaceAlpha',0.7)
colormap(colormap_rip)
caxis([-2,2])
%%
figure,
colormap_states = plasma(100);
color_gray = [134, 138, 145 ] ./ 255;
NREM_indices = pupil_traceZ{1}(end:-1:1,2);
pupil_label_NREM = pupil_traceZ{1}(end:-1:1,3);

scatter(LFP_features_umap(state_indices_gather,1),LFP_features_umap(state_indices_gather,2),20,color_gray,'filled')
hold on
scatter(LFP_features_umap(NREM_indices,1),LFP_features_umap(NREM_indices,2),30*ones(size(pupil_label_NREM)),pupil_label_NREM,'filled',...
    'MarkerEdgeAlpha',0.7,'MarkerFaceAlpha',0.7)
colormap(colormap_rip)
caxis([-2,2])

%% for supplemental figure, raw LFP trace
figure,
scatter3(LFP_features(state_indices_gather,3),LFP_features(state_indices_gather,5),LFP_features(state_indices_gather,6),50*ones(size(state_indices_gather)),pupil_trace(state_indices_gather,2),'filled',...
    'MarkerEdgeAlpha',0.7,'MarkerFaceAlpha',0.7)
colormap(colormap_rip)
xlabel('theta power')
ylabel('PSD slope')
zlabel('EMG from LFP')
caxis([0.14,0.82])










