clc
clear all
close all
%% add codes to path 
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/Pupil_OptoBehav/';
OptoCtrl = 1; % 1 Opto, 0 Ctrl group
animalname = 'PPP12';
day = 10;

daystring = num2str(day);
if OptoCtrl
    animaldir = [dir,'/Opto/',animalname,'/day',daystring,'/'];
else
    animaldir = [dir,'/Ctrl/',animalname,'/day',daystring,'/'];
end
prefix = ['day',daystring];        
%% define parameters
posfilt = gaussian(5*0.5, 10); % gaussian smoothing kernal
max_posdiff = 5;
bin = 0.02; % 20ms
%% load behavior info
load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
%% training trials
fig = figure;
plot_number = 0;

traintrial_num = length(behavior.train(:,1));
pathlength = zeros(traintrial_num,1);
for tr = 1:traintrial_num
    trialtime = behavior.train(tr,:);
    trialid = find(behavior.timestamps >= trialtime(1) & behavior.timestamps <= trialtime(2));
    trial_pos_raw = [behavior.position.x(trialid)',behavior.position.y(trialid)'];
    trial_pos_time_raw = behavior.timestamps(trialid)';
    trial_pos_timebin = trial_pos_time_raw(1):bin:trial_pos_time_raw(end);
   
    % detect NaNs
    validid = find(~isnan(trial_pos_raw(:,1) + trial_pos_raw(:,2)));
    trial_pos_time = trial_pos_time_raw(validid);
    trial_pos = trial_pos_raw(validid,:);

    % replace jumping point
    % smooth
    posx_smooth = filtfilt(posfilt,1,trial_pos(:,1));
    diffposx = diff(posx_smooth); % get running direction on x-axis

    posy_smooth = filtfilt(posfilt,1,trial_pos(:,2));
    diffposy = diff(posy_smooth); % get running direction on y-axis

    posx_diff = abs(posx_smooth - trial_pos(:,1));
    posy_diff = abs(posy_smooth - trial_pos(:,1));

    invalid_id  = find(posx_diff > max_posdiff | posy_diff > max_posdiff);

    trial_pos(invalid_id,1) = posx_smooth(invalid_id);
    trial_pos(invalid_id,2) = posy_smooth(invalid_id);

    % interpolation for the same bin size
    trial_posx_interp = interp1(trial_pos_time,trial_pos(:,1),trial_pos_timebin);
    trial_posy_interp = interp1(trial_pos_time,trial_pos(:,2),trial_pos_timebin);

    % plot
    plot_number = plot_number +1;
    supersubplot(fig, 3, 3, plot_number);
    plot(trial_pos_raw(:,1),trial_pos_raw(:,2))
    hold on
    plot(trial_posx_interp,trial_posy_interp)

    
    % calculate path length
    posx_diff = diff(trial_posx_interp).^2;
    posy_diff = diff(trial_posx_interp).^2;
    pos_diff = sqrt(posx_diff + posy_diff);
    pathlength(tr) = nansum(pos_diff);

    % update behavior structure
    behavior.train_pos{tr} = [trial_pos_timebin',trial_posx_interp',trial_posy_interp'];
    behavior.train_pathlength = pathlength;
end

%% test trials
fig = figure;
plot_number = 0;

testtrial_num = length(behavior.test(:,1));
pathlength = zeros(testtrial_num,1);
for tr = 1:testtrial_num
    trialtime = behavior.test(tr,:);
    trialid = find(behavior.timestamps >= trialtime(1) & behavior.timestamps <= trialtime(2));
    trial_pos_raw = [behavior.position.x(trialid)',behavior.position.y(trialid)'];
    trial_pos_time_raw = behavior.timestamps(trialid)';
    trial_pos_timebin = trial_pos_time_raw(1):bin:trial_pos_time_raw(end);

  
    % detect NaNs
    validid = find(~isnan(trial_pos_raw(:,1) + trial_pos_raw(:,2)));
    trial_pos_time = trial_pos_time_raw(validid);
    trial_pos = trial_pos_raw(validid,:);

    % replace jumping point
    % smooth
    posx_smooth = filtfilt(posfilt,1,trial_pos(:,1));
    diffposx = diff(posx_smooth); % get running direction on x-axis

    posy_smooth = filtfilt(posfilt,1,trial_pos(:,2));
    diffposy = diff(posy_smooth); % get running direction on y-axis

    posx_diff = abs(posx_smooth - trial_pos(:,1));
    posy_diff = abs(posy_smooth - trial_pos(:,1));

    invalid_id  = find(posx_diff > max_posdiff | posy_diff > max_posdiff);

    trial_pos(invalid_id,1) = posx_smooth(invalid_id);
    trial_pos(invalid_id,2) = posy_smooth(invalid_id);

    % interpolation for the same bin size
    trial_posx_interp = interp1(trial_pos_time,trial_pos(:,1),trial_pos_timebin);
    trial_posy_interp = interp1(trial_pos_time,trial_pos(:,2),trial_pos_timebin);

    % plot
    plot_number = plot_number +1;
    supersubplot(fig, 3, 3, plot_number);
    plot(trial_pos_raw(:,1),trial_pos_raw(:,2))
    hold on
    plot(trial_posx_interp,trial_posy_interp)

    % calculate path length
    posx_diff = diff(trial_posx_interp).^2;
    posy_diff = diff(trial_posx_interp).^2;
    pos_diff = sqrt(posx_diff + posy_diff);
    pathlength(tr) = nansum(pos_diff);

    % update behavior structure
    behavior.test_pos{tr} = [trial_pos_timebin',trial_posx_interp',trial_posy_interp'];
    behavior.test_pathlength = pathlength;
end
%%
save(fullfile(animaldir,[prefix,'.animal.behavior.mat']),'behavior')