% this is the second step to generate Fig. 1h - calculating SI based on the
% ripple LFP within NREM
% python toolbox from the Prida Lab is needed, get it here: https://github.com/PridaLab/structure_index
clc
clear all
close all
%% define data
dir = '/Volumes/Pupil_data/';
% list of all data, [{animal},{day},{[epochs]}]
animal_info = [{'HYC2'},{2},{[2,3]};...
    {'HYC3'},{8},{1};...
    {'HYC3'},{9},{1};...
    {'PPP4'},{8},{3};...
    {'PPP4'},{10},{5};...
    {'PPP4'},{11},{[1,3,5]};...
    {'PPP7'},{8},{[3,5]};...
    {'PPP7'},{12},{[1,4]};...
    {'PPP7'},{14},{1};...
    {'PPP8'},{7},{[1,3,5]};...
    {'PPP8'},{8},{[1,3]}];
%% define parameters
nshuffles = 500; % number of shuffles
savedata = 1; % save data

si_pymod = py.importlib.import_module('structure_index'); % import the python toolbox for structure index calculation
%% animal-day loop
for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    eps = animal_info{session_list,3};
    daystring = num2str(day);
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];
    
    % load LFP features
    lfp_ripples_all = [];
    pupil_label_all = [];
    for ep = eps
        disp([animalname,' Day-',num2str(day),' Ep-',num2str(ep)])
        load(fullfile(animaldir,[prefix,'.RippleLFPfeatures-NoRS-','EP',num2str(ep),'.LFP.mat']));
        lfp_ripples_all  = [lfp_ripples_all;lfp_ripples];
        pupil_label_all = [pupil_label_all;ripple_triggerZ_pupil_mean];
    end
    %% calculate structure index
    LFP_features_trans = lfp_ripples_all';
    features_py = py.numpy.array(LFP_features_trans(:).');
    features_py  = py.numpy.reshape(features_py ,int64(size(LFP_features_trans')));

    label_py = py.numpy.array(pupil_label_all(:).');

    result = si_pymod.compute_structure_index(pyargs(...
        'data',features_py,'label',label_py,...
        'num_shuffles',int16(nshuffles),'n_neighbors',int16(100)));
    Structure_index = double(result(1));
    %% save results
    shuf_SI_vals = result{4};
    shuf_SI_vals = double(py.array.array('d',py.numpy.nditer(shuf_SI_vals)));

    if savedata
        save(fullfile(animaldir,[prefix,'.rawStructureIndex-RipplePupil-NoRS','.mat']),'lfp_ripples_all','pupil_label_all','Structure_index','shuf_SI_vals');
    end
    close all
end
