clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';
% list of all data, in [animal, day, epochs];
animal_info = [{'HYC2'},{2},{[2,3]};...
    {'HYC3'},{8},{1};...
    {'HYC3'},{9},{1};...
    {'PPP4'},{8},{3};...
    {'PPP4'},{9},{5};...
    {'PPP4'},{10},{5};...
    {'PPP4'},{11},{[1,3,5]};...
    {'PPP4'},{18},{5};...
    {'PPP7'},{8},{[3,5]};...
    {'PPP7'},{12},{[1,4]};...
    {'PPP7'},{14},{1};...
    {'PPP7'},{23},{[3,6]};...
    {'PPP8'},{7},{[1,3,5]};...
    {'PPP8'},{8},{[1,3]}];

state_IDs = [1,3,5];%1 awake, 3 NREM, 5 REM
%% animal-day loop
for state = 1:3
    pupil_state{state} = [];
end
for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    eps = animal_info{session_list,3};
    daystring = num2str(day);
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];

    load(fullfile(animaldir,[prefix,'.MergePoints.events.mat'])) % Session info

    % pupil info
    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % time, normalized diameter

    % sleep state
    load(fullfile(animaldir,[prefix,'.SleepState.states.mat'])) % sleep states
    Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM
    %% epoch loop
    for ep = eps
        disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])
        epochtimes = MergePoints.timestamps(ep,:);

        valididx = find(pupil_trace(:,1) >= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
        pupil_trace_ep = pupil_trace(valididx,:);
        
        % zscore
        pupil_trace_ep(:,2) = (pupil_trace_ep(:,2) - nanmean(pupil_trace_ep(:,2)))./nanstd(pupil_trace_ep(:,2));

        % get NREM episodes
        sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
        Sleep_state_ep = Sleep_state(sleepss_id,:);

        for state = 1:3 %1 awake, 3 NREM, 5 REM
            current_state = state_IDs(state);
            State_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
            Stateid = find(Sleep_state_ep(:,2) == current_state);
            State_vec_ep(Stateid) = 1;
            Statelist = vec2list(State_vec_ep, Sleep_state_ep(:,1));
            %% restrcit pupil to the current state
            pupil_trace_temp = [];
            if ~isempty(Statelist)
                for seg = 1:length(Statelist(:,1))
                    validid = find(pupil_trace_ep(:,1) >= Statelist(seg,1) & pupil_trace_ep(:,1) <= Statelist(seg,2));  % get time stamps
                    pupil_trace_temp = [pupil_trace_temp;pupil_trace_ep(validid,2)];
                end
            else
                pupil_trace_temp = nan;
            end
            pupil_state{state} = [pupil_state{state};pupil_trace_temp];
        end
    end
end
%% histogram
edges = -5:0.1:5;
for state = 1:3
    temp = histcounts(pupil_state{state},edges);
    pupil_hist{state}= temp'./length(~isnan(pupil_state{state}));
end
bin = edges(1:end-1)' + (edges(2)-edges(1))/2;
