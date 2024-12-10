% this is the main script to generate Fig. 2c
clc
clear all
close all
%% load data
load('GLMgain_SWRproperties_pupil_all_results.mat') % preprocess results
ripple_triggerZ_pupil = pupil_info;
ripple_triggerZ_pupil_mean = nanmean(ripple_triggerZ_pupil(:,301:end),2);
repplayPval_triggerZ_pupil = ripple_info(:,5);
%% define parameters
pval_threshold = 0.025; % replay threshold
samplenum = 4000; % randomly take 4000 events
nRun = 100; % run 100 times
%% seperate replay by pupil
Y = prctile(ripple_triggerZ_pupil_mean,1/6*100);
tile1 = [min(ripple_triggerZ_pupil_mean),Y];
tile2 = [Y,prctile(ripple_triggerZ_pupil_mean,2/6*100)];
tile3 = [prctile(ripple_triggerZ_pupil_mean,2/6*100),prctile(ripple_triggerZ_pupil_mean,3/6*100)];
tile4 = [prctile(ripple_triggerZ_pupil_mean,3/6*100),prctile(ripple_triggerZ_pupil_mean,4/6*100)];
tile5 = [prctile(ripple_triggerZ_pupil_mean,4/6*100),prctile(ripple_triggerZ_pupil_mean,5/6*100)];
tile6 = [prctile(ripple_triggerZ_pupil_mean,5/6*100),prctile(ripple_triggerZ_pupil_mean,6/6*100)];
sixtiles = [tile1;tile2;tile3;tile4;tile5;tile6];
for i = 1:6
    currentid = find(ripple_triggerZ_pupil_mean >= sixtiles(i,1) & ripple_triggerZ_pupil_mean <= sixtiles(i,2));
    replayPval_tile{i} = repplayPval_triggerZ_pupil(currentid);
end
%% shuffled data
replayPval_triggerZ_pupil_shuf = repplayPval_triggerZ_pupil(randperm(length(reactivation_triggerZ_pupil)));

for i = 1:6
    currentid = find(ripple_triggerZ_pupil_mean >= sixtiles(i,1) & ripple_triggerZ_pupil_mean <= sixtiles(i,2));
    replayPval_tile_shuf{i} = replayPval_triggerZ_pupil_shuf(currentid);
end
%% estimate replay percentage for each of pupil sextile
for run = 1:nRun
    for i = 1:6
        randid = randperm(length(reactivation_tile{i}));
        temp = reactivation_tile{i}(randid(1:samplenum));
        replaynum(i) = length(find(temp < pval_threshold));
        nonreplaynum(i) = length(find(temp >= pval_threshold));
        replayprc(run,i) = replaynum(i)./samplenum;
    end
end
%% shuffled percentage
for run = 1:nRun
    for i = 1:6
        randid = randperm(length(reactivation_tile_shuf{i}));
        temp = reactivation_tile_shuf{i}(randid(1:samplenum));
        replaynum_shuf(i) = length(find(temp < pval_threshold));
        replayprc_shuf(run,i) = replaynum_shuf(i)./samplenum;
    end
end