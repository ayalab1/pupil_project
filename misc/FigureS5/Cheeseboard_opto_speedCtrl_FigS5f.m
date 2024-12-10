clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/Pupil_OptoBehav/';
% small pupil data set, [animal, day]
Opto_animal_info = [{'PPP13'},{[7]};...
    {'PPP13'},{[7]};
    {'PPP14'},{[6]};...
    {'PPP14'},{[6]};...
    {'PPP15'},{[3]};...
    {'PPP15'},{[5]};...
    {'PPP16'},{[5]};...
    {'PPP12'},{[10]}];

% large pupil data set
% PPP13 day8 missed video tracking, so it is excluded
Ctrl_animal_info = [{'PPP13'},{[4]};...
    {'PPP13'},{[6]};...
    {'PPP14'},{[3]};...
    {'PPP14'},{[7]};...
    {'PPP15'},{[6]};...
    {'PPP16'},{[16]};...
    {'PPP12'},{[11]}];
%% set parameters
do_smooth = 1; % smooth the positions
order = 2; %the order for Kalman filter for running velocity
%% Opto group (small pupil)
opto_vel = [];
for session_list = 1:length(Opto_animal_info)
    animalname = Opto_animal_info{session_list,1};
    day = Opto_animal_info{session_list,2};
    daystring = num2str(day);
    animaldir = [dir,'/Opto/',animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];        
    %% load behavior info
    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    %% calculate running speed
    opto_pos{1} = behavior.train_pos{end};
    opto_pos{2} = behavior.test_pos{1};
    opto_pos{3} = behavior.test_pos{2};
    for i = 1:3
        positions = opto_pos{i};
        %find valid post
        post = positions(:,1);
        posx = positions(:,2);
        posy = positions(:,3);
        validid = find(~isnan(posx) & ~isnan(posy));
        post_val = post(validid);
        posx_val = posx(validid);
        posy_val = posy(validid);

        %take care of the edge issue
        if post_val(end) < post(end) % the last pos value is missing
            post_val(end+1) = post(end);
            posx_val(end+1) = posx_val(end);
            posy_val(end+1) = posy_val(end);
        end

        if post_val(1) > post(1) % the first pos value is missing
            post_val = [post(1);post_val];
            posx_val = [posx_val(1);posx_val];
            posy_val = [posy_val(1);posy_val];
        end

        posx_pad = interp1(post_val,posx_val,post,'pchip');
        posy_pad = interp1(post_val,posy_val,post,'pchip');

        if do_smooth
            %smooth
            posfilt = gaussian(10*0.5, 20); % gaussian smoothing for velocity filter
            posx_pad = filtfilt(posfilt,1,posx_pad);
            posy_pad = filtfilt(posfilt,1,posy_pad);
        end
       
        [~,~,~,vx,vy,~,~] = KalmanVel(posx_pad,posy_pad,post,order);%calculate speed with Kalman filters
        vel_temp = sqrt(vx.^2+vy.^2);

        vel(i) = nanmean(vel_temp);
    end
    opto_vel = [opto_vel;vel];
end

%% Ctrl group (large pupil)
ctrl_vel = [];
for session_list = 1:length(Ctrl_animal_info)
    animalname = Ctrl_animal_info{session_list,1};
    day = Ctrl_animal_info{session_list,2};
    daystring = num2str(day);
    animaldir = [dir,'/Ctrl/',animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];       
    %% load behavior info
    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    %% calculate running speed
    opto_pos{1} = behavior.train_pos{end};
    opto_pos{2} = behavior.test_pos{1};
    opto_pos{3} = behavior.test_pos{2};
    for i = 1:3
        positions = opto_pos{i};
        %find valid post
        post = positions(:,1);
        posx = positions(:,2);
        posy = positions(:,3);
        validid = find(~isnan(posx) & ~isnan(posy));
        post_val = post(validid);
        posx_val = posx(validid);
        posy_val = posy(validid);

        %take care of the edge issue
        if post_val(end) < post(end) % the last pos value is missing
            post_val(end+1) = post(end);
            posx_val(end+1) = posx_val(end);
            posy_val(end+1) = posy_val(end);
        end

        if post_val(1) > post(1) % the first pos value is missing
            post_val = [post(1);post_val];
            posx_val = [posx_val(1);posx_val];
            posy_val = [posy_val(1);posy_val];
        end

        posx_pad = interp1(post_val,posx_val,post,'pchip');
        posy_pad = interp1(post_val,posy_val,post,'pchip');

        if do_smooth
            %smooth
            posfilt = gaussian(10*0.5, 20); % gaussian smoothing for velocity filter
            posx_pad = filtfilt(posfilt,1,posx_pad);
            posy_pad = filtfilt(posfilt,1,posy_pad);
        end

        [~,~,~,vx,vy,~,~] = KalmanVel(posx_pad,posy_pad,post,order);%calculate speed with Kalman filters
        vel_temp = sqrt(vx.^2+vy.^2);

        vel(i) = nanmean(vel_temp);
    end
    ctrl_vel = [ctrl_vel;vel];
end







