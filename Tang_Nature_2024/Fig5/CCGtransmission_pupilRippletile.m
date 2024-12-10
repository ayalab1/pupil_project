function  pTransmission_tile = CCGtransmission_pupilRippletile(dir,animalname,animalprefix,day,epoch)
%---------------------------------------------------------------%
%  This is the function for calculating PYR-INT spike           %
%  transmission probability across pupil tiles.                 %
%  -- Wenbo Tang (Oct 16, 2023)                                 %
%---------------------------------------------------------------%
%
% INPUTS:
%
%    dir = data directory.
%    animalname = animal name.
%    animalprefix = animal prefix, usually the same as experimental day in string.
%    day = experimental day string.
%    epoch = epoch.
%
% OUTPUTS:
%
%    pTransmission_tile = [pTs for tiles 1-3, cell pair UIDs];
%
% Calls CCG.m and ce_GetTransProb.m from AYALab neurocode
%% define parameters
binSize = 4/10000; duration = 0.12; % cell excplorer ce_MonoSynConvClick parameters
%% load data
% load spikes
spikes = GetAyaSpikes(dir);
if  strcmp(animalname,'PPP7')
    spikes_pyr = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"dCA1");
    spikes_int = importSpikes('basepath',dir,'CellType',["Narrow Interneuron", "Wide Interneuron"],'brainRegion',"dCA1");
elseif strcmp(animalname,'PPP4') || (strcmp(animalname,'PPP8') && day == 8) || (strcmp(animalname,'PPP8') && day == 15)
    spikes_pyr = importSpikes('basepath',dir,'CellType',"Pyramidal Cell");
    spikes_int = importSpikes('basepath',dir,'CellType',["Narrow Interneuron", "Wide Interneuron"]);
else
    spikes_pyr = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"CA1");
    spikes_int = importSpikes('basepath',dir,'CellType',["Narrow Interneuron", "Wide Interneuron"],'brainRegion',"CA1");
end

load(fullfile(dir,[animalprefix,'.MergePoints.events.mat'])) % Session info
epochtimes = MergePoints.timestamps(epoch,:);

load(fullfile(dir,[animalprefix,'.cell_metrics.cellinfo.mat'])) % cell info

load([dir,animalprefix,'.PupilRipples','EP',num2str(epoch),'.events.mat']) % pupil info

% extract monosynaptic cell pairs
monosynaptic_list = cell_metrics.putativeConnections.excitatory;
EI_monosynaptic_list = [];
for pair = 1:length(monosynaptic_list(:,1))
    current_pair_UIDs = monosynaptic_list(pair,:);
    if ismember(current_pair_UIDs(1),spikes_pyr.UID) && ismember(current_pair_UIDs(2),spikes_int.UID)
        EI_monosynaptic_list = [EI_monosynaptic_list;current_pair_UIDs];
    end
end
%% calculate spike transmission probability
if ~isempty(EI_monosynaptic_list)
    for tile = 1:3 % small, middle, large pupil tile
        ripIDs = find(ripples.pupil_tile == (tile*2-1)| ripples.pupil_tile == tile*2);
        
        interval = [ripples.timestamps(ripIDs,1),...
            ripples.timestamps(ripIDs,1) + ripples.duration(ripIDs)];
        restricted = Restrict(spikes,interval);

        % CCG.m from AYALab neurocode
        [ccg,~] = CCG(restricted(:,1),restricted(:,2),'binSize',binSize,'duration',duration);

        for i = 1:size(EI_monosynaptic_list,1)
            try
                ccgEI = ccg(:,EI_monosynaptic_list(i,1),EI_monosynaptic_list(i,2));
                spikesI = sum(restricted(:,2) == EI_monosynaptic_list(i,1));
                spikesE = sum(restricted(:,2) == EI_monosynaptic_list(i,2));
                if spikesI > 100 % at least 100 spikes to calculate robust pT
                    % ce_GetTransProb.m from AYALab neurocode
                    [trans,~,~,~] = ce_GetTransProb(ccgEI,  spikesI,  binSize,  0.020);
                    pTransmission_tile(i,tile) = trans;
                else %otherwise, set to NaN
                    pTransmission_tile(i,tile)  = nan;
                end
                    
            catch
                pTransmission_tile(i,tile)  = nan;
            end
            
        end
    end
    pTransmission_tile(:,4:5) = EI_monosynaptic_list; % save UID info
else
    pTransmission_tile = [];
end
