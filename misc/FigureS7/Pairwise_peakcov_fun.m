function [timebase,Zcrosscov_sm,peakcov,peakcov_time] = Pairwise_peakcov_fun(spikes1,spikes2,timelist)
% profile on
bin = 0.01; % 10 ms
tmax = 0.5; % +/- 500ms for corrln
sw1 = bin*3; % for smoothing corrln. Necessary?

% bincc = 0.1; % 100ms bin for cc
tmaxcc = 0.2; % -200 to 200 ms extent

if ~isempty(timelist)
    timedur = sum(timelist(:,2)-timelist(:,1));

    spikes1 = spikes1(find(isExcluded(spikes1,timelist)));
    spikes2 = spikes2(find(isExcluded(spikes2,timelist)));

    if (~isempty(spikes1) && ~isempty(spikes2))
%     if (length(spikes1)> 10 && length(spikes2) > 10)
           spkstime1 = double(spikes1);
           spkstime2 = double(spikes2);

           xc = spikexcorr(spkstime1, spkstime2, bin, tmax);
           p1 = xc.nspikes1/timedur; p2 = xc.nspikes2/timedur; % Fir rate in Hz
           exp_p = p1*p2; % per sec
           crosscov = (xc.c1vsc2 ./ (bin*timedur))-exp_p;
           % Convert to Z-score
           factor = sqrt((bin*timedur) / exp_p);
           Zcrosscov = crosscov .* (factor);
           rawcorr = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2);

           nspikes1 = xc.nspikes1;
           nspikes2 = xc.nspikes2;
           nstd=round(sw1/(xc.time(2) - xc.time(1))); % will be 3 std
           g1 = gaussian(nstd, nstd);
           timebase = xc.time;
           bins_run = find(abs(timebase) <= tmaxcc); % +/- Corrln window

           % Zcrosscov_sm = smoothvect(Zcrosscov, g1);% target
           Zcrosscov_sm = smoothvect(rawcorr, g1);% target


           currthetacorr_sm = nanmean(Zcrosscov_sm,1);
           [peakcov,peakcov_id] = nanmax(currthetacorr_sm(bins_run)); % Already smoothened, or can take +/-3 bins around peak
           peakcov_time = timebase(bins_run(peakcov_id));
%            thetacov = 0.5*log((1+thetacov)/(1-thetacov));
    %        tmpcorr = spikexcorr(spkstime1, spkstime2, bincc, tmaxcc);
    %        normcorr_cal  = tmpcorr.c1vsc2 ./ sqrt(tmpcorr.nspikes1 * tmpcorr.nspikes2);
    %        thetacov = max(normcorr_cal); % no bias-correction (see Jones and Wilson 2005)
    else
        timebase = nan;
        Zcrosscov_sm = nan;
        peakcov = nan;
        peakcov_time = nan;
    end
else
    timebase = nan;
    Zcrosscov_sm = nan;
    peakcov = nan;
    peakcov_time = nan;
end




