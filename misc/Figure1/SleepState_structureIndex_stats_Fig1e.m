clc
clear all
close all
%% add codes to path
addpath(genpath('/Users/wenbotang/Documents/Remote_HMM_files/AYALab_Code/'))
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

NREMorALL = 1; % 1 for NREM only, 0 for ALL states
%% gather result
structure_index_all = [];
for session_list = 1:length(animal_info)
    animalname = animal_info{session_list,1};
    day = animal_info{session_list,2};
    daystring = num2str(day);
    animaldir = [dir,animalname,'/day',daystring,'/'];
    prefix = ['day',daystring];
    
    % load structure index
    if NREMorALL
        % NREM only
        % load([animaldir,prefix,'.StructureIndex-','.mat']) % for UMAP
        load([animaldir,prefix,'.rawStructureIndex-','.mat']) % for raw trace

    else
        % all states
        % load([animaldir,prefix,'.StructureIndex-ALLStatesN20','.mat']) % for UMAP
        load([animaldir,prefix,'.rawStructureIndex-ALLStatesN20','.mat']) % for raw trace

    end
    structure_index_all = [structure_index_all;Structure_index, prctile(shuf_SI_vals,95)];
end
    
