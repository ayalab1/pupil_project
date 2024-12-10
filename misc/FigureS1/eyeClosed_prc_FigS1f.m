clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';
% list of all data
animal_info = [{'HYC2'},{2},{2};...
    {'HYC2'},{2},{3};...
    {'HYC3'},{8},{1};...
    {'HYC3'},{9},{1};...
    {'PPP4'},{8},{3};...
    {'PPP4'},{9},{5};...
    {'PPP4'},{10},{5};...
    {'PPP4'},{11},{1};...
    {'PPP4'},{11},{3};...
    {'PPP4'},{11},{5};...
    {'PPP4'},{18},{5};...
    {'PPP7'},{8},{3};...
    {'PPP7'},{8},{5};...
    {'PPP7'},{12},{1};...
    {'PPP7'},{12},{4};...
    {'PPP7'},{14},{1};...
    {'PPP7'},{23},{3};...
    {'PPP7'},{23},{6};...
    {'PPP8'},{7},{1};...
    {'PPP8'},{7},{3};...
    {'PPP8'},{7},{5};...
    {'PPP8'},{8},{1};...
    {'PPP8'},{8},{3}];
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
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];
    
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
    
    
    
                

    
  
   