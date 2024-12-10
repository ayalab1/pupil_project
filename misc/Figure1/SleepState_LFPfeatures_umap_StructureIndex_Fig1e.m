clc
clear all
close all
%%
% python toolbox from the Prida Lab is needed, get it here: https://github.com/PridaLab/structure_index
addpath(genpath('/Users/wt248/Documents/Remote_HMM_files/AYALab_Code/'))
addpath(genpath('/Users/wt248/Documents/MATLAB/umapAndEppFileExchange (4.1)'))
%% define data
dir = '/Volumes/Extreme Pro/Pupil_data/';

% list of all data
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

%% define parameters
nshuffles = 500; % number of shuffles
savedata = 1; % save data

si_pymod = py.importlib.import_module('structure_index'); % import the python toolbox for structure index calculation
%%  animal-day loop
for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    eps = animal_info{session_list,3};
    daystring = num2str(day);
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];

    % load LFP features
    LFP_features_all = [];
    pupil_label_all = [];
    for ep = eps
        disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])
        load(fullfile(animaldir,[prefix,'.SleepStateLFPfeatures-','EP',num2str(ep),'.LFP.mat']));
        LFP_features_all = [LFP_features_all;LFP_features_NREM];
        pupil_label_all = [pupil_label_all;pupil_label];
    end
    
    % UMAP
    LFP_features_NREM_umap = run_umap(LFP_features_all(:,2:end),'n_components', 2,'min_dist',0.6,'n_neighbors',20,'metric','cosine');

    %% calculate structure index
    LFP_features_NREM_umap_trans = LFP_features_NREM_umap';
    features_py = py.numpy.array(LFP_features_NREM_umap_trans(:).');
    features_py  = py.numpy.reshape(features_py ,int64(size(LFP_features_NREM_umap_trans')));

    label_py = py.numpy.array(pupil_label_all(:).');

    result = si_pymod.compute_structure_index(pyargs(...
        'data',features_py,'label',label_py,...
        'num_shuffles',int16(nshuffles)));
    Structure_index = double(result(1));
    %% save results
    shuf_SI_vals = result{4};
    shuf_SI_vals = double(py.array.array('d',py.numpy.nditer(shuf_SI_vals)));

    if savedata
        save(fullfile(animaldir,[prefix,'.StructureIndex-','.mat']),'LFP_features_NREM_umap','pupil_label','Structure_index','shuf_SI_vals');
    end
    close all
end
