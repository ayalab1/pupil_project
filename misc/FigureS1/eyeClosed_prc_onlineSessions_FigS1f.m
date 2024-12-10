clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/Pupil_OnlineOffline/';
% list of all data
animal_info = [{'PPP13'},{3},{1};...
    {'PPP13'},{5},{1};...
    {'PPP13'},{7},{1};...
    {'PPP13'},{8},{1};...
    {'PPP10'},{13},{1};...
    {'PPP10'},{14},{1};...
    {'PPP14'},{2},{1};...
    {'PPP14'},{4},{1};...
    {'PPP14'},{6},{1};...
    {'PPP15'},{3},{1}];
%% set parameters
mindur = 5; %5s
%% gather data
% reset
eyeclose_prc = [];
for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    ep = animal_info{session_list,3};
    daystring = num2str(day);

    if strcmp(animalname,'PPP13') || strcmp(animalname,'PPP14')
        animaldir = [dir,animalname,'/day',daystring,'_post','/'];
        prefix = ['day',daystring,'_post'];
    elseif strcmp(animalname,'PPP10')
        animaldir = [dir,animalname,'/day',daystring,'/'];
        prefix = ['day',daystring];
    elseif strcmp(animalname,'PPP15')
        animaldir = [dir,animalname,'/day',daystring,'/post_sleep'];
        prefix = 'post_sleep';
    end
    
    disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])

    % pupil info
    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,2]); % time, normalized diameter
    epochtimes = [behavior.epochs{ep}.startTime,behavior.epochs{ep}.stopTime];  
    session_dur = epochtimes(2) - epochtimes(1);
    
    tempid  = find(pupil_trace(:,1) >= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
    pupil_trace = pupil_trace(tempid,:);
    
    Pupil_SR = nanmedian(diff(pupil_trace(:,1)));
    a = find(isnan(pupil_trace(:,2)));
    [lo,hi]= findcontiguous(a);  %find contiguous NaNs
    eyeclose_time = 0;
    for i = 1:length(lo)-1
        if (hi(i)-lo(i)) > mindur/Pupil_SR
            eyeclose_time = eyeclose_time + (hi(i)-lo(i))*Pupil_SR;
        end
    end
    eyeclose_prc = [eyeclose_prc;double(eyeclose_time/session_dur)];
end
    
    
    
                

    
  
   