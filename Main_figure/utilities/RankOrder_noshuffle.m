function [rankStats] = RankOrder_noshuffle(varargin)
% [rankStats] = RankOrder_noshuffle(varargin)
%
%  Get rank order of spikes inside events from previously calculated 
% bz_getRipSpikes structure. The 'corrEvents' output field structure shows
% rank correlation of each event with the selected 'templateType', and the
% 'pvalEvents' field states how significantly different from chance is that
% correlation, so presumably those events with low 'pvalEvents' will have a
% rank order that repeats over time. In the case that there is more than 
% one rank sequence, the output field structure 'rankClusters' will number
% the events equally if they have similar sequences.
%
% INPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'basepath'	    full path where session is located (default pwd)
%                       E.g. /mnt/Data/buddy140_060813_reo/buddy140_060813_reo
%     'spkEventTimes'	Field to be used is:    
%                   		.EventRel: 2xN cell matrix. In the first row, relative 
%                   		times of spikes for that particular event across all units.
%                   		In the second row, the UID associated to the above spike
%     'templateType'    String that can be:
%                       - 'Pairwise': (default) it computes the rank correlation
%                           of each event pairwise with each other event.
%                       - 'MeanRank': it computes the rank correlation 
%                           of each event against the mean rank
%                           of the rest of the events.
%                       - 'Peak': searchs for the bz_findPlaceFieldsTemplate
%                           output, the X.placeFieldTemplate.mat, loads it
%                           and takes the 'Peak' field, a 1 x (# conditions) 
%                           cell array. Within each cell there is a (# units) x 3
%                           matrix, which has bins corresponding to the firing 
%                           map peak for each unit (NaN for units without place
%                           field); second column has the unit ID; third column
%                           has the position of  the unit in the UID vector. The
%                           rest of the cells contain this similar matrix for the
%                           other conditions.
%                       - 'CenterofMass': searchs for the bz_findPlaceFieldsTemplate
%                           output, the X.placeFieldTemplate.mat, loads it
%                           and takes the 'CenterofMass', a 1 x (# conditions) 
%                           cell array. Within each cell there is a (# units) x 3
%                           matrix, which has bins corresponding to the firing 
%                           map peak for each unit (NaN for units without place 
%                           field); second column has the unit ID; third column
%                           has the position of  the unit in the UID vector. 
%                           The rest of the cells contain this similar matrix 
%                           for the other conditions.
%     'eventIDs'        A (#events)-length array of 0s and 1s that indicates
%                       which event belongs to one of two different groups 
%                       of events. By default all events will belong to a
%                       single group, and the rank correlation of each event
%                       will be performed against the rest of the events.
%                       However, if an 'eventIDs' array establishes two
%                       different groups of events (e.g. spontaneous vs induced
%                       ripples), then the rank correlation of events in 
%                       one group will be performed against the rank of the
%                       other group.
%     'timeSpike'       A string, to determine what time reference of spike:
%                       1. 'first': (default) takes into account the first time
%                           the unit fires 
%                       2. 'mean': takes the mean of the spike time
%     'minUnits'        Minimum number of units parcitipating in each event.
%                       Events with less than 'minUnits' are not considered
%                       for correlation.
%     'normalize'       Work with normalized ranks, logical (default: true)
%                       So instead of ranking 1, 4, 5, 2, 3, it will be 0,
%                       0.8, 1.0, 0.2, 0.4.
%     'doPlot'	        Make plots, logical (default: true) 
%     'saveMat'	        Saves file, logical (default: true) 
%
%    =========================================================================
%
% OUTPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'rankStats'	    structure with the following subfields:
%         .corrMean          (#conditions x 1) Mean of rank correlations 
%                            of all events. If there are several conditions
%                            each row is the mean correlation of events
%                            against the template given by each condition.
%         .corrStd           (#conditions x 1) Standard deviation of rank 
%                            of all events. If there are several conditions
%                            each row is the mean correlation of events
%                            against the template given by each condition.
%         .corrEvents        (#conditions x #events matrix) Rank correlation
%                            of each event. Same as with 'corrMean' happens
%                            for each row.
%         .rankUnits         (#units x #events matrix) Rank of units for
%                            each event. This is the matrix that has been
%                            used to compute the correlation.
%                            E.g.:  A  [ 0.0  0.4  0.6  nan       0.6 
%                                   B    0.2  0.2  nan  nan       nan
%                                   C    0.4  0.0  nan  nan       nan
%                                   D    0.6  0.6  1.0  nan  ...  1.0
%                                   E    0.8  0.8  0.3  nan       0.3
%                                   F    1.0  1.0  nan  nan       nan ]
%  See also bz_getRipSpikes, bz_findPlaceFieldsTemplate
%

%  RankOrder.m by Andrea Navas-Olive, 2019. Antonio FR, 2020. 
%  Modified by Wenbo Tang, 2023, removed shuffling procedures from original 
%  RankOrder.m to speed up the calculation.

%% Parse inputs 
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'spkEventTimes',{},@isstruct);
addParameter(p,'templateType','MeanRank',@isstr);
addParameter(p,'eventIDs',1,@isnumeric);
addParameter(p,'timeSpike','first',@isstr);
addParameter(p,'minUnits',5,@isnumeric);
addParameter(p,'normalize',true,@islogical);
addParameter(p,'doPlot', true, @islogical);
addParameter(p,'saveMat', true, @islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
spkEventTimes = p.Results.spkEventTimes;
templateType = p.Results.templateType;
eventIDs = p.Results.eventIDs;
normalize = p.Results.normalize;
timeSpike = p.Results.timeSpike;
minUnits = p.Results.minUnits;
doPlot = p.Results.doPlot;
saveMat = p.Results.saveMat;

% Get session info
basename = basenameFromBasepath(basepath);
% load([basepath filesep basename '.sessionInfo.mat']);
% Load default spkEventTimes
if isempty(spkEventTimes)
    spkEventTimes = load([basepath filesep basename '.spkEventTimes.mat']);
    spkEventTimes = spkEventTimes.spkEventTimes;
end
% Relative times of spikes for each particular event across all units.
evtTimes = spkEventTimes.EventRel;

% If template type is 'Peak' or 'CenterofMass' load 'placeFieldTemplate'
if strcmp(templateType,'Peak') || strcmp(templateType,'CenterofMass')
    % External template
    templateExt = load([basepath filesep basename '.placeFieldTemplate.mat']);
    templateExt = templateExt.placeFieldTemplate;
    templateExt = templateExt.(templateType);
else
    templateExt = [];
end


%% Rank matrix
% Create (#events) x (#units) matrix with position of spikes in event. It
% considers just the first spikes of each unit
rankUnits = nan*ones(size(spkEventTimes.UnitEventRel));
for event = 1:length(evtTimes)
    % Take into account just first spike
    if strcmp(timeSpike,'first')
        units = unique(evtTimes{event}(2,:),'stable');
    elseif strcmp(timeSpike,'mean')
        means = [];
        for jj = unique(evtTimes{event}(2,:))
            means  = [means [mean(evtTimes{event}(1,evtTimes{event}(2,:)==jj)); jj]];
        end
        if ~isempty(means)
            means = sortrows(means')';
            units = means(2,:);
        else
            units = [];
        end
    else
        warning('The variable "timeSpike" is invalid');
    end
    nUnits = length(units);
    % Set event as nan if it has no enough units
    if nUnits < minUnits
        rankUnits(:,event) = nan;
    % Rank units
    else
        rankUnits(units,event) = 1:nUnits;
        % If normalize, set ranks between 0 and 1        
        if normalize
            rankUnits(units,event) = rankUnits(units,event) / nUnits;
        end
    end
end

%% Compute correlation

% Compute mean and standard deviation of rank correlation, and correlation
% event by event
[corrMean, corrStd, corrEvents] = compute_rank_correlation(templateType, rankUnits, evtTimes, templateExt, eventIDs);

%% Save and plot

% Save statistics

    rankStats.corrMean = corrMean;
    rankStats.corrStd = corrStd;
    rankStats.corrEvents = corrEvents;
    rankStats.rankUnits = rankUnits;

end




% RANK CORRELATION...
function [corrMean, corrStd, corrEvents] = compute_rank_correlation(templateType, rankUnits, evtTimes, templateExt, eventIDs)


    % ... WITHOUT EXTERNAL TEMPLATE
    
    % Template method: a template for each specific ripple or theta event is
    %   constructed based on the averaged rank of all units over all other events
    %   (i.e., excluding that specific event). Then the rank correlation between
    %   each specific event and its template was computed, and averaged over all
    %   events
    if strcmp(templateType,'MeanRank')
        corrEvents = zeros(size(evtTimes));
        % - If there are no two groups of events
        if length(unique(eventIDs))==1
            parfor event = 1:length(evtTimes)
                % Order of units in this event
                rankThisEvent = rankUnits(:,event);
                % Mean order of units along the rest of the events (template)
                rankTemplate = nanmean(rankUnits(:,[1:event-1,event+1:end]),2);
                % Correlation between these previous variables
                corrEvents(event) = corr(rankThisEvent, rankTemplate, 'rows', 'complete');
            end
        % - If there are two groups of events, test this event against 
        %   the meank rank of the other group
        else
            parfor event = 1:length(evtTimes)
                % Order of units in this event
                rankThisEvent = rankUnits(:,event);
                % Mean order of units along the rest of the events (template)
                rankTemplate = nanmean(rankUnits(:,eventIDs==(1-eventIDs(event))),2);
                % Correlation between these previous variables
                corrEvents(event) = corr(rankThisEvent, rankTemplate, 'rows', 'complete');
            end
        end
        % Mean and standard deviation of rank correlation
        corrMean = nanmean(corrEvents);
        corrStd = nanstd(corrEvents);
        
    % Pairwise method:
    elseif strcmp(templateType,'Pairwise')
        corrEvents = corr(rankUnits, 'rows', 'pairwise');
        corrEvents(eye(size(corrEvents))==1) = nan;
        % Mean and standard deviation of rank correlation...
        % - If there are no two groups of events
        if length(unique(eventIDs))==1
            corrEvents = nanmean(corrEvents);
        % - If there are two groups of events, test this event against 
        %   the meank rank of the other group
        else
            corrEventsMean = zeros(size(evtTimes));
            parfor event = 1:length(evtTimes)
                corrEventsMean(event) = nanmean(corrEvents(event,eventIDs==(1-eventIDs(event))));
            end
            corrEvents = corrEventsMean;
        end
        % Mean and standard deviation of rank correlation
        corrMean = nanmean(corrEvents);
        corrStd = nanstd(corrEvents);
    end

    
    % ... WITH EXTERNAL TEMPLATE
    % Template method: external template
    if ismember(templateType,{'Peak','CenterofMass'})
        % Initialize
        corrEvents = zeros(length(templateExt),size(evtTimes,2));
        corrMean = zeros(length(templateExt),1);
        corrStd = zeros(length(templateExt),1);
        for iCond = 1:length(templateExt)
            parfor event = 1:size(evtTimes,2)
                % Order of units in this event
                rankThisEvent = rankUnits(:,event);
                % External template (we need to transform it to an array
                % similar to rankThisEvent: a (# units) x 1 array, where in
                % the position of each unit the rank is written 
                rankTemplate = nan*ones(size(rankThisEvent));
                rankTemplate(templateExt{iCond}(:,3)) = templateExt{iCond}(:,1);
                % Normalize
                if max(max(rankUnits))==1
                    rankTemplate = (rankTemplate - nanmin(templateExt{iCond}(:,1)))/ (nanmax(templateExt{iCond}(:,1))- nanmin(templateExt{iCond}(:,1)));
                end
                % Correlation between these previous variables
                corrEvents(iCond,event) = corr(rankThisEvent, rankTemplate, 'rows', 'complete');
            end
        end
        % Mean and standard deviation of rank correlation
        corrMean = nanmean(corrEvents,2);
        corrStd = nanstd(corrEvents,[],2);
    end
end

