clc
close all
clear all;
%% add codes to path
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
addpath(genpath('/Users/wt248/Documents/MATLAB/matplotlib'))
%% setting for the plots
% defaultGraphicsSetttings
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';

animalprefix = 'PPP7';
day = 8;
daystring = num2str(day);
dir = [dir,animalprefix,'/day',daystring,'/'];
prefix = ['day',daystring];
ep = 4;
%% set parameters
conditions = 2; % trajectory types, left vs. right
pos_binnum = 100;
pos_bin = 1:100;
%% gather data
for track = 1:conditions
    rm1 = []; % ratemap matrix
    peakrate1 = []; % peak rate matrix
    peakloc1 = []; % peak location matrix
    sparsity1 = [];
    % load spikes
    if (strcmp(animalprefix,'PPP8') && day == 8)
        spikes = importSpikes('basepath',dir,'CellType',"Pyramidal Cell");
    elseif strcmp(animalprefix,'PPP4') || strcmp(animalprefix,'PPP7')
        spikes = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"dCA");
    else
        spikes = importSpikes('basepath',dir,'CellType',"Pyramidal Cell",'brainRegion',"CA1");
    end

    cellnum = length(spikes.UID);

    %-----load the ratemaps-----%
    load(fullfile(dir,  [prefix,'.firingMapsAvg_linear.cellinfo.mat']));% load firing rate maps
    %% cell loop
    nn = 0;
    for i = 1:cellnum
        UID = spikes.UID(i);
        cind = find(firingMaps_linear{ep}.UID == UID);  
        rm_cell = [];
        peakrate_cell = [];
        peakloc_cell = [];
        if ~isempty(cind)
            linfield1 = firingMaps_linear{ep}.rateMaps{cind};  
        else
            linfield1=[];
        end

        if ~isempty(linfield1)
            linfield_hp = linfield1{track};
            a = find(linfield_hp < 0);

            %pad nan
            if ~isempty(a)
                [lo,hi]= findcontiguous(a);  %find contiguous NaNs
                for ii = 1:length(lo)
                    if lo(ii) > 1 & hi(ii) < length(linfield_hp)
                        fill = linspace(linfield_hp(lo(ii)-1), ...
                            linfield_hp(hi(ii)+1), hi(ii)-lo(ii)+1);
                        linfield_hp(lo(ii):hi(ii)) = fill;
                    end
                end
            end

        else
            linfield_hp = zeros(size(1:pos_binnum));
        end
        [rm1_peak,rm1_peakloc] = max(linfield_hp);
        rm1_peakloc = pos_bin(rm1_peakloc);
        linfield_hp(isnan(linfield_hp)) = 0;
        idx = find(linfield_hp > 0.25*rm1_peak);
        sparsity = length(idx)/length(linfield_hp);
        rm1 = [rm1;linfield_hp];
        peakrate1 = [peakrate1;rm1_peak];
        peakloc1 = [peakloc1;rm1_peakloc];
        sparsity1 = [sparsity1;sparsity];
    end
    rm{track}.ratemap = rm1;
    rm{track}.peakrate = peakrate1;
    rm{track}.peakloc = peakloc1;
    rm{track}.sparsity = sparsity1;
end
%%  sort rate maps by peak locations on the left trajectory
validid = find(rm{1}.peakrate > 3 | rm{2}.peakrate > 3);

figure,
[peakloc_sorted,sortedid] = sort(rm{1}.peakloc(validid));

% normalized by peak rate
ratemap2 = rm{2}.ratemap(validid(sortedid),:)./rm{2}.peakrate(validid(sortedid));
ratemap1 = rm{1}.ratemap(validid(sortedid),:)./rm{1}.peakrate(validid(sortedid));
ratemap1(isnan(ratemap1)) = 0;
ratemap2(isnan(ratemap2)) = 0;

zerocellid = find(peakloc_sorted == 0);
allrate = sum(ratemap1(zerocellid,:),2);
[~,indx2]  = sort(allrate);
ratemap2(zerocellid,:) = ratemap2(zerocellid(indx2),:);
ratemap1(zerocellid,:) = ratemap1(zerocellid(indx2),:);
%%
figure
m=100;
cm_inferno = inferno(m);
subplot(1,2,1)
imagesc(ratemap1)
colormap(cm_inferno)
caxis([0.2 1])
subplot(1,2,2)
imagesc(ratemap2)
colormap(cm_inferno)
caxis([0.2 1])

%% by peak locations on the right trajectory
figure,
[peakloc_sorted,sortedid] = sort(rm{2}.peakloc(validid));

% normalized by peak firing rate
ratemap2 = rm{2}.ratemap(validid(sortedid),:)./rm{2}.peakrate(validid(sortedid));
ratemap1 = rm{1}.ratemap(validid(sortedid),:)./rm{1}.peakrate(validid(sortedid));
ratemap1(isnan(ratemap1)) = 0;
ratemap2(isnan(ratemap2)) = 0;

zerocellid = find(peakloc_sorted == 0);
allrate = sum(ratemap2(zerocellid,:),2);
[~,indx2]  = sort(allrate);
ratemap2(zerocellid,:) = ratemap2(zerocellid(indx2),:);
ratemap1(zerocellid,:) = ratemap1(zerocellid(indx2),:);


subplot(1,2,1)
imagesc(ratemap1)
colormap(cm_inferno)
caxis([0.2 1])
subplot(1,2,2)
imagesc(ratemap2)
colormap(cm_inferno)
caxis([0.2 1])

