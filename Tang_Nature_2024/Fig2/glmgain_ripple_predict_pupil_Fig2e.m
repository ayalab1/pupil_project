% this is the main script to generate Fig. 2e
clc
clear all
close all
%% load data
load('GLMgain_SWRproperties_pupil_all.mat') % preprocess results
auxiliary_variables = ripple_info(:,[1,3,5]); % [amplitude,RS,replay Pval]
auxiliary_variables = zscore(auxiliary_variables);
pupil_infoZ_sum = sum(pupil_info,2);
num_features  = length(auxiliary_variables(1,:));
%% define parameters
data_prc = 0.8;% randomly take 80% for  training, and remaining 20% for testing
run_num = 500; % run 500 time
data_num = round(data_prc * length(pupil_infoZ_sum));
%% loop
for run = 1:run_num
    run
    randid_label = randperm(length(pupil_info(:,1))); % randomly permute data 
    data_sample_ids = randid_label(1:data_num); 
    [gain,shGain,~,~,~,~] = GLMgain(auxiliary_variables(data_sample_ids,:),pupil_infoZ_sum(data_sample_ids),'dist','normal','link','identity');
    gain_z(run) = (gain' - nanmean(shGain))./nanstd(shGain);
    shGain_95th_z(run) = (prctile(shGain,95)- nanmean(shGain))./nanstd(shGain);
    %% blate one variable and test again
    for removed_dim = 1:num_features
        [gain_blated,shGain_blated,predictions,er,sh,w] = GLMgain_control(auxiliary_variables(data_sample_ids,:),pupil_infoZ_sum(data_sample_ids),removed_dim,'dist','normal','link','identity');
        gain_blated_z = (gain_blated' - nanmean(shGain))./nanstd(shGain);
        gain_blated_gather(run, removed_dim)= gain_blated_z;
    end
end
%% save data
save('GLMgain_SWRproperties_pupil_all_results.mat')