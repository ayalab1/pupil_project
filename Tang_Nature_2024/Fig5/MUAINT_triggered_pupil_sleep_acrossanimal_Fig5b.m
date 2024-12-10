% this is the code to generate Fig. 5b and EDFig. 10g
clc
clear all
close all
%% define data
dir = '/Volumes/Pupil_data/';

% list of all data, [{animal},{day},{epoch}]
animal_info = [{'HYC3'},{8},{1};...
    {'HYC3'},{9},{1};...
    {'PPP4'},{8},{3};...
    {'PPP4'},{10},{5};...
    {'PPP4'},{11},{1};...
    {'PPP4'},{11},{3};...
    {'PPP4'},{11},{5};...
    {'PPP7'},{8},{3};...
    {'PPP7'},{8},{5};...
    {'PPP7'},{12},{1};...
    {'PPP7'},{12},{4};...
    {'PPP7'},{14},{1};...
    {'PPP8'},{7},{1};...
    {'PPP8'},{7},{3};...
    {'PPP8'},{7},{5};...
    {'PPP8'},{8},{1};...
    {'PPP8'},{8},{3}];
%% define parameters
win = -30:0.1:30; % window around SWRs, in sec
%% gather data
ripple_trigger_pupil = []; % reset
rippleamp_trigger_pupil = [];% ripple amplitude
MUA_trigger_pupil = []; % PYR
MUAINT_trigger_pupil = []; %INT

for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    ep = animal_info{session_list,3};
    daystring = num2str(day);
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];
    
    disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])

    load(fullfile(animaldir,[prefix,'.MergePoints.events.mat'])) % Session info
    epochtimes = MergePoints.timestamps(ep,:);

    % pupil info
    load(fullfile(animaldir,[prefix,'.animal.behavior.mat'])) % behavior
    pupil_trace = behavior.pupil.pupil_sleep_info(:,[1,3]); % [time, normalized pupil diameter]

    % sleep states
    load(fullfile(animaldir,[prefix,'.SleepState.states.mat'])) % sleep states
    Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  %1 awake, 3 NREM, 5 REM
    sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
    Sleep_state_ep = Sleep_state(sleepss_id,:);
    SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
    SWSid = find(Sleep_state_ep(:,2) == 3); %NREM
    SWS_vec_ep(SWSid) = 1;
    SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));
    SWSdur = SWSlist(:,2) - SWSlist(:,1);
    SWSlist = SWSlist(find(SWSdur > 10),:);

    % ripples 
    if strcmp(prefix,'day28') && strcmp(animalname,'PPP4')
        load(fullfile(animaldir,[prefix,'.ripples_task.events.mat'])) % ripple
    elseif strcmp(prefix,'day8') && strcmp(animalname,'PPP7')
        load(fullfile(animaldir,[prefix,'.dorsalripples.events.mat'])) % ripple
    else
        load(fullfile(animaldir,[prefix,'.ripples.events.mat'])) % ripple
    end
    
    riplist_all = [ripples.timestamps(:,1),ripples.timestamps(:,1) + ripples.duration,ripples.amplitude]; % [start, end, amplitude]
    
    % PYR MUA firing rates
    filename = [animaldir,'MUA', '0',num2str(day),'-0',num2str(ep),'.mat'];
    load(filename)
    MUA_FR = [MUA.time',MUA.data(:,1)];
    
    % INT MUA firing rates
    filename = [animaldir,'MUA_INT', '0',num2str(day),'-0',num2str(ep),'.mat'];
    load(filename)
    MUAINT_FR = [MUA.time',MUA.data(:,1)];
    %% restrict to NREM
    if ~isempty(SWSlist)
        pupil_trace_id = find(pupil_trace(:,1)>= epochtimes(1) & pupil_trace(:,1) <= epochtimes(2));
        pupil_trace = pupil_trace(pupil_trace_id,:);

        validid = find(riplist_all(:,1) >= pupil_trace(1,1) & riplist_all(:,2) <= pupil_trace(end,1));  % get time stamps
        riplist_ep = riplist_all(validid,1:2);
        riplist_amp_ep = riplist_all(validid,3);

        riplist = [];
        riplist_amp = [];
        riplist_MUA = [];
        riplist_MUAINT = [];
        for rip = 1:length(validid)
            ripid = find(MUA_FR(:,1) >= riplist_ep(rip,1) & MUA_FR(:,1) <= riplist_ep(rip,2));
            if ~isempty(ripid)
                riplist = [riplist;riplist_ep(rip,:)];
                riplist_amp = [riplist_amp;riplist_amp_ep(rip)];
                riplist_MUA = [riplist_MUA;nanmean(MUA_FR(ripid,2))];
                riplist_MUAINT = [riplist_MUAINT;nanmean(MUAINT_FR(ripid,2))];
            end
        end
                    
        % restrict ripples to NREM
        riplist_ep = [];
        riplist_amp_ep = [];
        riplist_MUA_ep = [];
        riplist_MUAINT_ep = [];
        for swsseg = 1:length(SWSlist(:,1))
            validid = find(riplist(:,1) >= SWSlist(swsseg,1) & riplist(:,2) <= SWSlist(swsseg,2));  % get time stamps
            riplist_ep = [riplist_ep;riplist(validid,:)];
            riplist_amp_ep = [riplist_amp_ep;riplist_amp(validid)];
            riplist_MUA_ep = [riplist_MUA_ep;riplist_MUA(validid)];            
            riplist_MUAINT_ep = [riplist_MUAINT_ep;riplist_MUAINT(validid)];
        end
        
        % zscore firing rate, because each session has different cells
        % recorded     
        riplist_ampZ_ep = (riplist_amp_ep - nanmean(riplist_amp_ep))./nanstd(riplist_amp_ep);
        riplist_MUAZ_ep = (riplist_MUA_ep - nanmean(riplist_MUA_ep))./nanstd(riplist_MUA_ep);
        riplist_MUAINTZ_ep = (riplist_MUAINT_ep - nanmean(riplist_MUAINT_ep))./nanstd(riplist_MUAINT_ep);
        %% get pupil info for each ripple
        for i = 1:length(riplist_ep)
            rippletime =  riplist_ep(i,:); 
            tempid  = find(pupil_trace(:,1) >= rippletime(1)-40 & pupil_trace(:,1) <= rippletime(1)+40);
            pupil_d = pupil_trace(tempid,2);
            pupil_time = pupil_trace(tempid,1)-rippletime(1);
            if length(find(~isnan(pupil_d))) > 10
                pupil_d_interp = interp1(pupil_time(~isnan(pupil_d)),pupil_d(~isnan(pupil_d)),win);
                ripple_trigger_pupil = [ripple_trigger_pupil;pupil_d_interp];
                rippleamp_trigger_pupil = [rippleamp_trigger_pupil;riplist_ampZ_ep(i)];
                MUA_trigger_pupil = [MUA_trigger_pupil;riplist_MUAZ_ep(i)];
                MUAINT_trigger_pupil = [MUAINT_trigger_pupil;riplist_MUAINTZ_ep(i)];
            end
        end
    end
end
%% normalized to baseline
MUA_ratio = MUAINT_trigger_pupil -MUA_trigger_pupil; % INT - PYR
baseline_mean = nanmean(ripple_trigger_pupil(:,1:300),2);
baseline_std = nanstd(ripple_trigger_pupil(:,1:300)')';
ripple_triggerZ_pupil = (ripple_trigger_pupil - repmat(baseline_mean,1,length(win)))./repmat(baseline_std,1,length(win));
ripple_triggerZ_pupil_sum = sum(ripple_triggerZ_pupil,2);
validid = find(~isnan(ripple_triggerZ_pupil_sum));

ripple_triggerZ_pupil = ripple_triggerZ_pupil(validid,:);
MUA_ratio = MUA_ratio(validid);
rippleamp_triggerZ_pupil = rippleamp_trigger_pupil(validid);
ripple_triggerZ_pupil_mean = nanmean(ripple_triggerZ_pupil(:,301:end),2);
MUAINT_trigger_pupil = MUAINT_trigger_pupil(validid);
MUA_trigger_pupil  = MUA_trigger_pupil(validid);
%% sorted pupil by INT-PYR MUA ratio
Y = prctile(MUA_ratio,1/6*100);
tile1 = [min(MUA_ratio),Y];
tile2 = [Y,prctile(MUA_ratio,2/6*100)];
tile3 = [prctile(MUA_ratio,2/6*100),prctile(MUA_ratio,3/6*100)];
tile4 = [prctile(MUA_ratio,3/6*100),prctile(MUA_ratio,4/6*100)];
tile5 = [prctile(MUA_ratio,4/6*100),prctile(MUA_ratio,5/6*100)];
tile6 = [prctile(MUA_ratio,5/6*100),prctile(MUA_ratio,6/6*100)];
sixtiles = [tile1;tile2;tile3;tile4;tile5;tile6];
figure
ripple_trigger_pupil_win_stats = [];
for i = 1:6
    currentid = find(MUA_ratio >= sixtiles(i,1) & MUA_ratio <= sixtiles(i,2));
    ripple_triggerZ_pupil_tile = ripple_triggerZ_pupil(currentid,:);
    temp = ripple_triggerZ_pupil_tile(:,250:350);
    ripple_trigger_pupil_win_stats = [ripple_trigger_pupil_win_stats;nanmean(temp(:)),nanstd(temp(:))',sum(double((~isnan(temp(:)))))'];
    ripple_trigger_pupil_stats{i} = [win',nanmean(ripple_triggerZ_pupil_tile)',nanstd(ripple_triggerZ_pupil_tile)',sum(double((~isnan(ripple_triggerZ_pupil_tile))))'];
    plot(win',nanmean(ripple_triggerZ_pupil_tile)')
    hold on
end
%% sorted MUA rate by pupil size
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
    ripple_amp_tile{i} = rippleamp_triggerZ_pupil(currentid);
    MUA_tile{i} = MUA_ratio(currentid);
    MUAINT_tile{i} = MUAINT_trigger_pupil(currentid);
    MUAPYR_tile{i} = MUA_trigger_pupil(currentid);

end
%% sorted MUA rate by ripple amplitudes
Y = prctile(rippleamp_triggerZ_pupil,1/6*100);
tile1 = [min(rippleamp_triggerZ_pupil),Y];
tile2 = [Y,prctile(rippleamp_triggerZ_pupil,2/6*100)];
tile3 = [prctile(rippleamp_triggerZ_pupil,2/6*100),prctile(rippleamp_triggerZ_pupil,3/6*100)];
tile4 = [prctile(rippleamp_triggerZ_pupil,3/6*100),prctile(rippleamp_triggerZ_pupil,4/6*100)];
tile5 = [prctile(rippleamp_triggerZ_pupil,4/6*100),prctile(rippleamp_triggerZ_pupil,5/6*100)];
tile6 = [prctile(rippleamp_triggerZ_pupil,5/6*100),prctile(rippleamp_triggerZ_pupil,6/6*100)];
sixtiles = [tile1;tile2;tile3;tile4;tile5;tile6];

for i = 1:6
    currentid = find(rippleamp_triggerZ_pupil >= sixtiles(i,1) & rippleamp_triggerZ_pupil <= sixtiles(i,2));
    MUA_tile_amp{i} = MUA_ratio(currentid);
    MUAINT_tile_amp{i} = MUAINT_trigger_pupil(currentid);
    MUAPYR_tile_amp{i} = MUA_trigger_pupil(currentid);
end
%% calculate boostrapped distribution of MUA-pupil correlation, for EDFig. 10
[r_all_pupil,p_all_pupil] = corr(MUA_ratio,ripple_triggerZ_pupil_mean,'rows','complete','type','spearman');
Bootcorrfun = @(a,b)(corr(a,b,'rows','complete','type','spearman'));
bootstats_pupil = bootstrp(1000,Bootcorrfun,MUA_ratio,ripple_triggerZ_pupil_mean);
se_corr_pupil = std(bootstats_pupil);
mean_corr_pupil = mean(bootstats_pupil);
%% calculate boostrapped distribution of MUA-ripple amplitude correlation, for EDFig. 10
[r_all_ripple,p_all_ripple] = corr(rippleamp_triggerZ_pupil,ripple_triggerZ_pupil_mean,'rows','complete','type','spearman');
Bootcorrfun = @(a,b)(corr(a,b,'rows','complete','type','spearman'));
bootstats_amp = bootstrp(1000,Bootcorrfun,rippleamp_triggerZ_pupil,ripple_triggerZ_pupil_mean);
se_corr_amp = std(bootstats_amp);
mean_corr_amp = mean(bootstats_amp);

