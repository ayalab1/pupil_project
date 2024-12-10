clc
clear all
close all 
%% add codes to path
addpath(genpath('/Users/wenbotang/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data

%-----the sessions used are listed below:-----%
% "PPP13 day3"
% "PPP13 day5"
% "PPP13 day7" (example trace in Fig. 1a)
% "PPP13 day8"
% "PPP10 day13"
% "PPP10 day14"
% "PPP14 day2"
% "PPP14 day4"
% "PPP14 day6"
% "PPP15 day3"
% * only one epoch in each day
%----------------------------------------------%

dir = '/Volumes/Extreme Pro/Pupil_data/Pupil_OnlineOffline/';
animalprefix = 'PPP13';
day = 7;
ep = 1;
daystring = num2str(day);
figopt = 1; % plot the interim result? 
savedata = 1; % save the result
%%
if strcmp(animalprefix,'PPP13') || strcmp(animalprefix,'PPP14')
    dir = [dir,animalprefix,'/day',daystring,'_post','/'];
    prefix = ['day',daystring,'_post'];
elseif strcmp(animalprefix,'PPP10')
    dir = [dir,animalprefix,'/day',daystring,'/'];
    prefix = ['day',daystring];
elseif strcmp(animalprefix,'PPP15')
    dir = [dir,animalprefix,'/day',daystring,'/post_sleep'];
    prefix = 'post_sleep';
end
%% load data
load(fullfile(dir,[prefix,'.animal.behavior.mat'])) % behavior file
%% get pupil info from DLC offline tracking
epochtimes = [behavior.epochs{ep}.startTime,behavior.epochs{ep}.stopTime];
epoch_timevecid = find(behavior.pupil.timestamps >= epochtimes(1) & behavior.pupil.timestamps <= epochtimes(2));
epoch_timevec = behavior.pupil.timestamps(epoch_timevecid)';

[xlength,ylength] = size(behavior.pupil.eyenasal_x_point);
if xlength > ylength
    eyentxy = [behavior.pupil.eyenasal_x_point behavior.pupil.eyenasal_y_point behavior.pupil.eyetemporal_x_point behavior.pupil.eyetemporal_y_point];
    pupilxy = [behavior.pupil.pupilnasal_x_point behavior.pupil.pupilnasal_y_point behavior.pupil.pupiltemporal_x_point behavior.pupil.pupiltemporal_y_point];
else
    eyentxy = [behavior.pupil.eyenasal_x_point' behavior.pupil.eyenasal_y_point' behavior.pupil.eyetemporal_x_point' behavior.pupil.eyetemporal_y_point'];
    pupilxy = [behavior.pupil.pupilnasal_x_point' behavior.pupil.pupilnasal_y_point' behavior.pupil.pupiltemporal_x_point' behavior.pupil.pupiltemporal_y_point'];
end

eyentxy_ep = eyentxy(epoch_timevecid,:);
pupilxy_ep = pupilxy(epoch_timevecid,:);
%% pupil info from online tracking
if strcmp(animalprefix,'PPP10')
    if day == 14
        CurrentTime = behavior.onlinetracking.timestamps(1:length(behavior.onlinetracking.Area));
        for i = 1:length(CurrentTime)
            if isempty(behavior.onlinetracking.Area{i})
                pupil_area(i,1) = CurrentTime(i);
                pupil_area(i,2) = nan(1);
            else
                pupil_area(i,:) = [CurrentTime(i),str2num( behavior.onlinetracking.Area{i})];
            end
        end
    else
        pupil_area = [behavior.onlinetracking.timestamps(1:length(behavior.onlinetracking.Area)),behavior.onlinetracking.Area];
    end
elseif (strcmp(animalprefix,'PPP13') && (day == 7 || day == 8)) || (strcmp(animalprefix,'PPP14') && (day == 6))
    for i = 1:length(behavior.onlinetracking.CurrentTime)
        if isempty(behavior.onlinetracking.Area{i})
            pupil_area(i,1) = behavior.onlinetracking.CurrentTime(i);
            pupil_area(i,2) = nan(1);
        else
            pupil_area(i,:) = [behavior.onlinetracking.CurrentTime(i),str2num( behavior.onlinetracking.Area{i})];
        end
    end
else
    pupil_area = [behavior.onlinetracking.CurrentTime,behavior.onlinetracking.Area];
end
%% interpolate the missing values
% eye traces
for i = 1:4
    temp = eyentxy_ep(:,i);
    validid = find(~isnan(temp));
    invalidid = find(isnan(temp));
    temp = temp(validid);
    %take care of the edge issue
    if validid(end) < length(eyentxy_ep(:,i)) % the last pos value is missing
        validid(end+1) = length(eyentxy_ep(:,i));
        temp(end+1) = temp(end);
    end
    
    if validid(1) > 1 % the first pos value is missing
        validid = [1;validid];
        temp = [temp(1);temp];
    end
    temp_interp = interp1(validid,temp,(1:length(eyentxy_ep(:,i)))','pchip');
    
    if ~isempty(invalidid)
        [lo,hi]= findcontiguous(invalidid);  %find contiguous NaNs
        long_nan_period_id =  find(hi-lo > 50);
        for ii = long_nan_period_id' % don't pad long periods with NaNs
            temp_interp(lo(ii):hi(ii)) = nan;
        end
    end
    eyentxy_ep(:,i) = temp_interp;
end
% pupil traces
for i = 1:4
    temp = pupilxy_ep(:,i);
    validid = find(~isnan(temp));
    invalidid = find(isnan(temp));
    
    temp = temp(validid);
    %take care of the edge issue
    if validid(end) < length(pupilxy_ep(:,i)) % the last pos value is missing
        validid(end+1) = length(pupilxy_ep(:,i));
        temp(end+1) = temp(end);
    end
    
    if validid(1) > 1 % the first pos value is missing
        validid = [1;validid];
        temp = [temp(1);temp];
    end
    temp_interp = interp1(validid,temp,(1:length(pupilxy_ep(:,i)))','pchip');
    if ~isempty(invalidid)
        [lo,hi]= findcontiguous(invalidid);  %find contiguous NaNs
        long_nan_period_id =  find(hi-lo > 50);
        for ii = long_nan_period_id' % don't pad long periods with NaNs
            temp_interp(lo(ii):hi(ii)) = nan;
        end
    end
    pupilxy_ep(:,i) = temp_interp;
end

% online pupil traces
temp = pupil_area(:,2);
validid = find(~isnan(temp));
invalidid = find(isnan(temp));
temp = temp(validid);
%take care of the edge issue
if validid(end) < length(pupil_area(:,2)) % the last pos value is missing
    validid(end+1) = length(pupil_area(:,2));
    temp(end+1) = temp(end);
end

if validid(1) > 1 % the first pos value is missing
    validid = [1;validid];
    temp = [temp(1);temp];
end
temp_interp = interp1(validid,temp,(1:length(pupil_area(:,2)))','pchip');

if ~isempty(invalidid)
    [lo,hi]= findcontiguous(invalidid);  %find contiguous NaNs
    long_nan_period_id =  find(hi-lo > 50);
    for ii = long_nan_period_id' % don't pad long periods with NaNs
        temp_interp(lo(ii):hi(ii)) = nan;
    end
end
pupil_area(:,2) = temp_interp;

pupil_area_interp = interp1(pupil_area(:,1),pupil_area(:,2),epoch_timevec);
pupil_area = [epoch_timevec,pupil_area_interp];
%% check the tracking positions
if figopt
    figure,
    plot(eyentxy_ep(:,1),eyentxy_ep(:,2),'.')
    hold on
    plot(eyentxy_ep(:,3),eyentxy_ep(:,4),'.')

    plot(pupilxy_ep(:,1),pupilxy_ep(:,2),'.')
    plot(pupilxy_ep(:,3),pupilxy_ep(:,4),'.')
    legend('eye nasal','eye temporal','pupil nasal','pupil temporal')
end
%% detect outliers
%     eye_outlier = isoutlier(eyentxy_ep,"grubbs",'ThresholdFactor',1e-5);
%     pupil_outlier = isoutlier(pupilxy_ep,"grubbs",'ThresholdFactor',1e-4);
eye_outlier = isoutlier(eyentxy_ep,"movmean",5000); % use moving window 
pupil_outlier = isoutlier(pupilxy_ep,"movmean",5000);

pupil_area_outlier = isoutlier(pupil_area(:,2),"movmean",5000);
%% remove outliers
outlier_ids = any(pupil_outlier,2) | any(eye_outlier,2) | pupil_area_outlier;
eyentxy_ep_valid = eyentxy_ep(~outlier_ids,:);
pupilxy_ep_valid = pupilxy_ep(~outlier_ids,:);
epoch_timevec_valid = epoch_timevec(~outlier_ids);
pupil_area_valid = pupil_area(~outlier_ids,:);
%% plot results after removing outliers
if figopt
    figure,
    plot(eyentxy_ep_valid(:,1),eyentxy_ep_valid(:,2),'.')
    hold on
    plot(eyentxy_ep_valid(:,3),eyentxy_ep_valid(:,4),'.')

    plot(pupilxy_ep_valid(:,1),pupilxy_ep_valid(:,2),'.')
    plot(pupilxy_ep_valid(:,3),pupilxy_ep_valid(:,4),'.')
end
%% calcualte diameters
eye_d = sqrt( (eyentxy_ep_valid(:,1)-eyentxy_ep_valid(:,3)).^2 + (eyentxy_ep_valid(:,2)-eyentxy_ep_valid(:,4)).^2);
pupil_d = sqrt( (pupilxy_ep_valid(:,1)-pupilxy_ep_valid(:,3)).^2 + (pupilxy_ep_valid(:,2)-pupilxy_ep_valid(:,4)).^2);
pupil_d_norm = pupil_d./eye_d;
pupil_d = medfilt1(pupil_d,300);

pupil_d_online = sqrt(pupil_area_valid(:,2)./pi) * 2;
pupil_d_z = (pupil_d -nanmean(pupil_d))./nanstd(pupil_d);
pupil_d_online_z = (pupil_d_online -nanmean(pupil_d_online))./nanstd(pupil_d_online);
%% calculate the xcorr to fix the timestamp difference
validid = find(~isnan(pupil_d_z) & ~isnan(pupil_d_online_z));
[pupil_corr,lag] = xcorr(pupil_d_z(validid),pupil_d_online_z(validid),10000);

figure
plot(lag,pupil_corr)
[~,maxid] = max(pupil_corr);
maxlag = lag(maxid);
if maxlag > 0
   temp = [nan(maxlag,1);pupil_d_online(1:end-maxlag)];
   pupil_d_online = temp;
else
   maxlag = -maxlag;
   temp = [pupil_d_online(maxlag+1:end);nan(maxlag,1)];
   pupil_d_online = temp;
end
%% zscore for comparsion
pupil_d_z = (pupil_d -nanmean(pupil_d))./nanstd(pupil_d);
pupil_d_online_z = (pupil_d_online -nanmean(pupil_d_online))./nanstd(pupil_d_online);

% plot
figure('position',[100,100,800,300]),
plot(epoch_timevec_valid - epoch_timevec(1),pupil_d_z,'.')
title('pupil diameter')

hold on
plot(epoch_timevec_valid - epoch_timevec(1),pupil_d_online_z,'.')
%% calculate tracking errors
mse = nanmean((pupil_d_z-pupil_d_online_z).^2); % mean squared error, GLM no offset
rss = nanmean((pupil_d_z-nanmean(pupil_d_z)).^2); % squared error of spike train
tracking_perf = 1-mse/rss;
%% structure the data
pupil_d_offline = nan(size(epoch_timevec));
pupil_d_offline(~outlier_ids) = pupil_d;
pupil_d_online_all = nan(size(epoch_timevec));
pupil_d_online_all(~outlier_ids) = pupil_d_online;

pupil_info = [epoch_timevec,pupil_d_offline,pupil_d_online_all];

behavior.pupil.pupil_info = pupil_info;
behavior.pupil.pupil_info_fields = ...
            'time  diameter-DLC  diameter-Online';
behavior.pupil.tracking_mse = mse;
behavior.pupil.tracking_rss = rss;
behavior.pupil.tracking_perf = tracking_perf;
%% save the updated behavior file
if savedata
    save([dir,prefix,'.animal.behavior.mat'],'behavior') % update behavior file
end
