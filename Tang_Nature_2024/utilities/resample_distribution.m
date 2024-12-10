function [indices, y_resample] = resample_distribution(x,y)
    % Resample y to fit distribution of x
    % when sample number of y is larger than that of x
    % Wenbo Tang, 2024
    bin_min = min(min(x),min(y));
    bin_max = max(max(x),max(y));
    binnum = 100;
    bin_limits = linspace(bin_min, bin_max, binnum);
    %%
    [count_x,~]= histc(x, bin_limits);
    [count_y,binnumber_y]= histc(y, bin_limits);
    %%
    prob = count_x(binnumber_y);
    prob2 = count_y(binnumber_y);
    probabilities = prob ./ prob2;
    %%
    indices = randsample(1:length(y),length(x),true,probabilities);
    y_resample = y(indices);
end
