function  pTransmission_tile = CCGtransmission_pupiltile_ver2(dir,animalname,animalprefix,day,epoch)

%  calculate CCG for PYR-INT pairs for each pupil sextile
%  processes
%  Usage: pTransmission_tile = CCGtransmission_pupiltile_ver2(dir,animalname,animalprefix,day,epoch)
%  Inputs:
%  dir   (data directory) 
%  animalname    (animalname)
%  animalprefix    (datapath prefix)
%  day   (recording day)
%  epoch   (sleep session number)
% 
% Output:
% pTransmission_tile   (transmission probablity, N*8 - [pTs in 1-6 pupil tile, PYR UID, INT UID]).

binSize = 4/10000; duration = 0.12; % cell excplorer ce_MonoSynConvClick parameters

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

load([dir,animalprefix,'.PupilTile_NREM','EP',num2str(epoch),'.mat']) % pupil info

% sleep states
load(fullfile(dir,[animalprefix,'.SleepState.states.mat'])) % sleep states
Sleep_state = [SleepState.idx.timestamps,SleepState.idx.states];  
sleepss_id = find(Sleep_state(:,1) <= epochtimes(2) & Sleep_state(:,1) >= epochtimes(1));
Sleep_state_ep = Sleep_state(sleepss_id,:);
%1 awake, 3 NREM, 5 REM
SWS_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
SWSid = find(Sleep_state_ep(:,2) == 3);
SWS_vec_ep(SWSid) = 1;
SWSlist = vec2list(SWS_vec_ep, Sleep_state_ep(:,1));

WAKE_vec_ep = zeros(length(Sleep_state_ep(:,1)),1);
WAKEid = find(Sleep_state_ep(:,2) == 1);
WAKE_vec_ep(WAKEid) = 1;
WAKElist = vec2list(WAKE_vec_ep, Sleep_state_ep(:,1));

intervals{1} = SWSlist;
intervals{2} = WAKElist;
%% extract PYR-INT monosynaptic cell pairs
monosynaptic_list = cell_metrics.putativeConnections.excitatory;
EI_monosynaptic_list = [];
for pair = 1:length(monosynaptic_list(:,1)) 
    current_pair_UIDs = monosynaptic_list(pair,:);
    if ismember(current_pair_UIDs(1),spikes_pyr.UID) && ismember(current_pair_UIDs(2),spikes_int.UID)
        EI_monosynaptic_list = [EI_monosynaptic_list;current_pair_UIDs];
    end
end
%% calculate transmission prob during waking and NREM
if ~isempty(EI_monosynaptic_list)
    for period = 1:2
        interval = intervals{period};
        restricted = Restrict(spikes,interval);
        [ccg,~] = CCG(restricted(:,1),restricted(:,2),'binSize',binSize,'duration',duration);

        for i = 1:size(EI_monosynaptic_list,1)
            try
                ccgEI = ccg(:,EI_monosynaptic_list(i,1),EI_monosynaptic_list(i,2));
                spikesI = sum(restricted(:,2) == EI_monosynaptic_list(i,1));
                if spikesI > 100 % minimal number of spikes needed 
                    [trans,~,~,~] = ce_GetTransProb(ccgEI,  spikesI,  binSize,  0.020);
                    pTransmission(i,period) = trans;
                else
                    pTransmission(i,period)  = nan;
                end
            catch

                pTransmission(i,period)  = nan;
            end
        end
    end
    pTransmission(:,3:4) = EI_monosynaptic_list;
else
    pTransmission = [];
end
%% find the pair showing an decrease from WAKE to Sleep, and calculate their transmission probability for each pupil sextile
if ~isempty(pTransmission)
    dpT = pTransmission(:,2) - pTransmission(:,1); % wake - sleep
    valid_idx = find(dpT > 0);

    EI_monosynaptic_list_select = pTransmission(valid_idx,3:4);
    %%
    for tile = 1:6
        interval = pupil_tile{tile};
        restricted = Restrict(spikes,interval);
        [ccg,~] = CCG(restricted(:,1),restricted(:,2),'binSize',binSize,'duration',duration);

        for i = 1:size(EI_monosynaptic_list_select,1)
            try
                ccgEI = ccg(:,EI_monosynaptic_list_select(i,1),EI_monosynaptic_list_select(i,2));
                spikesI = sum(restricted(:,2) == EI_monosynaptic_list_select(i,1));
                [trans,~,~,~] = ce_GetTransProb(ccgEI,  spikesI,  binSize,  0.020);
                pTransmission_tile(i,tile) = trans;
            catch
                pTransmission_tile(i,tile)  = nan;
            end
            
        end
    end
    pTransmission_tile(:,7:8) = EI_monosynaptic_list_select; % add UID info to the output
else
    pTransmission_tile = [];
end

        



