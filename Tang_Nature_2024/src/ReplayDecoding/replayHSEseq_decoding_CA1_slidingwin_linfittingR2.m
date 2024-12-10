function replaytrajactory = replayHSEseq_decoding_CA1_slidingwin_linfittingR2(dir,animalname,animalprefix,day,ep,spikes,linfields,linfields_pos,cellcountthresh,figopt,savedata)
%---------------------------------------------------------------%
%  This is the function for replay decoding                     %
%  -- Wenbo Tang (Oct 16, 2023)                                 %
%  See also Shin*, Tang*, Jadhav, Neuron, 2019                  %
%---------------------------------------------------------------%
%
% INPUTS:
%
%    dir = data directory.
%    animalname = animal name.
%    animalprefix = animal prefix, usually the same as experimental day in string.
%    day = experimental day string.
%    ep = epoch.
%    spikes = spike times for each cell.
%    linfields = linearized place fields.
%    linfields_pos = position labels for linearized place fields.
%    cellcountthresh = mininum number of cells active for considering as a
%                      cadidate event, usually = 5
%    figopt = 1 or 0, 1 = plot the decoding result
%    savedata = 1 or 0, 1 = save the decoding results
%
% OUTPUTS:
%
%    replaytrajactory = decoding information for each SWR event
%% define save directory
savedir = '/Users/wt248/Replay_decoding_20231026/';
%% define parameters
tBinSz = 10; %default temporal bin in ms [typically 5-10ms]
tBinSz_sm = 20; %default temporal bin in ms [typically 20-40ms]
sBinSz = 2;%Spatial bin isze in cm
wellcutoff = 0; %cm, remove the edges?

%-----line fitting parameters----%
spdBnd = [100,10000]; % cm/s
nstep = 100;
spdBnd = log10(spdBnd);
spd2Test = logspace(spdBnd(1),spdBnd(2),nstep); % use logspace to reduce dim
spd2Test = [-1.*spd2Test(end:-1:1), spd2Test];
spd2Test        =spd2Test(~(spd2Test > -spdBnd(1) & spd2Test < spdBnd(1))); %Remove slow speeds
spd2Test        =(spd2Test ./sBinSz) .*tBinSz/1000;

Yinterc = -400:2:400;
halfwidth = 2;%8cm in total, 2cm/bin

nperm = 500;%shuffling 500 times
%% %-----create the ratemaps [nPosBin x nHPCells]-----%
hpnum = length(spikes.UID);
track_num = length(linfields{1}); %trajectory number
trackpos = linfields_pos';

rm = []; % ratemap matrix
pm = []; % position matrix
tm = []; % track matrix
trackid = [];
cellidxm = [];

for i = 1:hpnum
    linfield1 = linfields{i};
    linfield_hp = [];
    lintrack_hp = [];
    pos_hp = [];
    validtrack = zeros(1,track_num);
    for track = 1:track_num
        occnormrate1 = linfield1{track}';
        occnormrate1(find(occnormrate1 <0)) = nan;

        lintrack1 = ones(size(trackpos))*track;
        if max(occnormrate1) > 2 %Hz
            validtrack(track) = 1;
        else
            validtrack(track) = 0;
        end
        linfield_hp = [linfield_hp;occnormrate1];
        pos_hp = [pos_hp;trackpos];
        lintrack_hp = [lintrack_hp;lintrack1];
    end

    if (max(linfield_hp) >= 3) % peak firing rate max larger than 3 Hz
        a = find(isnan(linfield_hp));
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
        rm = [rm;linfield_hp'];
        pm = [pm;pos_hp'];
        tm = [tm;lintrack_hp'];
        trackid = [trackid;validtrack];
        cellidxm = [cellidxm; i];
    end
end

rm = rm'; %[nPosBin x nHPCells]
pm = pm';
tm = tm';
for i = 1:track_num
    pm_traj = pm(find(tm == i));
    maxpos = max(max(pm_traj));
    rm(find(tm == i & pm <= wellcutoff)) = 0;
    rm(find(tm == i & pm >= maxpos-wellcutoff)) = 0;
end
rm = rm+ (eps.^8); %Add a small number so there are no zeros
expecSpk =rm.*tBinSz_sm./1000; %[nPos x nCells] Expected number of spikes per bin
hpnum = length(rm(1,:)); % update
%% %-----create the event matrix during theta-----%
% get rippleHSE time
load(sprintf('%s/%srippleHSEs.events.mat', dir, animalprefix));% rippleHSE info
riptimes = [rippleHSEs{day}{ep}.starttime',rippleHSEs{day}{ep}.endtime'];

dur = 1000*(riptimes(:,2) - riptimes(:,1));
keepidx = find(dur >= 5*tBinSz);%at least 5 bins, 50 ms for 10ms bins
riptimes = riptimes(keepidx,:);

if ~isempty(riptimes)
    %% calculate candidate events
    celldata = [];
    spikecounts = [];
    for cellcount = 1:hpnum
        index = cellidxm(cellcount) ;
        validtrack = trackid(cellcount,:);
        if ~isempty(spikes.times{index})
            spiketimes = spikes.times{index};
        else
            spiketimes = [];
        end

        spikebins = periodAssign(spiketimes, riptimes(:,[1 2]));
        if ~isempty(spiketimes)
            validspikes = find(spikebins);
            spiketimes = spiketimes(validspikes);
            spikebins = spikebins(validspikes);
            tmpcelldata = [spiketimes spikebins];
        end
        if ~isempty(spiketimes)
            tmpcelldata(:,3) = cellcount;
            for t = 1:track_num
                tmpcelldata(:,3+t) = validtrack(t);
            end
        else
            tmpcelldata = [0 0 cellcount,validtrack];
        end
        celldata = [celldata; tmpcelldata];
        spikecount = zeros(1,size(riptimes,1));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        spikecounts = [spikecounts; spikecount];
    end

    [cellid,cyclenum] = find(spikecounts > 0);
    eventindex = [];
    for c = 1:length(spikecounts(1,:))
        cycid = find(cyclenum == c);
        trajinfo_tmp = sum(trackid(cellid(cycid),:),1);
        validtraj = find(trajinfo_tmp >= cellcountthresh);
        if ~isempty(validtraj)
            eventindex = [eventindex,c];
            validtrajinfo{length(eventindex)} = validtraj;
            for tj = validtraj
                ii = find(trackid(cellid(cycid),tj) > 0);
                activecellinfo{length(eventindex)}{tj} = cellid(cycid(ii));
            end
        end
    end

    %% looping events
    pvalue = nan(length(eventindex),4);
    pvalue_r = nan(length(eventindex),4);

    for event = 1:length(eventindex)
        event
        validtraj = validtrajinfo{event};
        for traj = validtraj
            cellsi = activecellinfo{event}{traj};
            event_cellSeq = cellsi;
            tmpind = find(celldata(:,2) == eventindex(event));
            spiketimes = celldata(tmpind,1);
            cellindex = celldata(tmpind,3);
            %-----create the event matrix during a theta cycle (spkT{cells}.spiketimes) -----%
            for ce = event_cellSeq'
                validspikeidx = find(cellindex == ce);
                spkT{ce} = spiketimes(validspikeidx).*1000;
            end
            %--------calculate the posterior probability on an event-by-event basis--------%
            startevent = riptimes(eventindex(event),1).*1000;
            endevent = riptimes(eventindex(event),2).*1000;
            eventtimes(event,1:2) = [startevent/1000,endevent/1000];

            timebins = startevent:tBinSz:endevent; % timebins are the binedges
            nTBin = length(timebins)-1;
            nCell = hpnum;
            spkPerBin = zeros(1,nTBin, nCell); % keep the inactive cells as 0s.
            for nn  = 1:hpnum
                cellInd = nn; %current cell
                if length(spkT) >= cellInd
                    if ~isempty(spkT{cellInd})
                        temp = histc(spkT{cellInd}, timebins); %[1 x nTBin x nCell]
                        temp2 = temp(1:end-1)+temp(2:end);
                        spkPerBin(1,:,cellInd) = temp2;
                    end
                end
            end
            nSpkPerTBin = squeeze(sum(spkPerBin,3)); %[nTBin x 1] number of spikes in tBin
            nPBin = size(rm,1); %N positional bin
            expecSpk  = reshape(expecSpk, [nPBin,1, nCell]); %[nPos x 1 x nCell]
            expon = exp(-expecSpk); %Exponent of equation.
            factSpkPerBin = factorial(spkPerBin); %Factorial to divide by

            wrking = bsxfun(@power, expecSpk, spkPerBin); %[nPos x nTbin x nCell]
            wrking = bsxfun(@rdivide, wrking, factSpkPerBin); %[nPos x nTbin x nCell]
            wrking = bsxfun(@times,wrking, expon); %[nPos x nTbin x nCell]
            post = prod(wrking,3); %Non normalised prob [nPos x Tbin]
            post(:,nSpkPerTBin==0)  =0; % FO changed this 27/5/15 so the posterior matrix can be smoothed.
            post(isnan(post)) = 0;
            trajinfo = mean(tm,2);% trajactory number
            figure(4)

            trajidx = find(trajinfo == traj);
            %pMat: [nPos, nTbin]
            pMat = post(trajidx,:);% create a posterior matrix for each traj
            distvector = pm(trajidx,1)';

            szPM1 = size(pMat,1);
            szPM2 = size(pMat,2);
            positions = 1:szPM1;

            for i = 1:szPM2; if (sum(pMat(:,i))>0)
                    pMat(:,i) = pMat(:,i)./sum(pMat(:,i));end;end % normalized across positions to 1 for each time bin
            nonzerobins = find(nSpkPerTBin > 0);
            timevec = 1:length(pMat(1,:));

            % calculate weighted correlation
            mloc = sum(sum(pMat.*positions'))/sum(sum(pMat));
            mt = sum(sum(pMat.*timevec))/sum(sum(pMat));
            dloc = positions'-mloc;
            dt = timevec - mt;
            cov_loc_t = sum(sum(pMat.*dloc.*dt))/sum(sum(pMat));
            cov_loc = sum(sum(pMat.*(dloc.^2)))/sum(sum(pMat));
            cov_t = sum(sum(pMat.*(dt.^2)))/sum(sum(pMat));
            rvalues(event,traj) = cov_loc_t/sqrt(cov_loc*cov_t);

            % line fitting
            tmp = spd2Test'* timevec;
            distvector = repmat((tmp)',1,length(Yinterc)) + reshape(repmat(Yinterc,length(tmp(1,:))*length(tmp(:,1)),1),length(tmp(1,:)),[]);
            aa = repmat((spd2Test'* ones(size(timevec)))',1,length(Yinterc));
            bb = reshape(repmat(Yinterc,length(tmp(1,:))*length(tmp(:,1)),1),length(tmp(1,:)),[]);            parameter_space = [aa(1,:);bb(1,:)];
            clear aa, clear bb;
            distvector_LB = distvector-halfwidth;
            distvector_HB = distvector+halfwidth;
            [pmax,~] = max(pMat);
            tid = find(pmax > 0);
            prob_sum = zeros(1,length(distvector(1,:)));
            icc_total = cell(1,length(nonzerobins));
            ibic_total = cell(1,length(nonzerobins));
            nn = 0;
            for t = nonzerobins
                nn = nn + 1;
                pMat_tmp = pMat(:,t);
                pos_ROI = round(distvector(t,:)'+(-halfwidth:halfwidth))';
                id = find(pos_ROI >= positions(1) & pos_ROI <= positions(end));
                ib = zeros(size(pos_ROI));
                ib(id) = pos_ROI(id);

                pmat_total = zeros(size(ib));
                ib = int32(ib);
                icc = find(ib > 0);
                ibic_total{nn} = ib(icc);
                icc_total{nn} = icc;
                pmat_total(icc) = pMat_tmp(ib(icc));

                prob_sum = prob_sum + sum(pmat_total);
            end
            probmax = max(prob_sum);
            maxid = find(prob_sum == probmax);

            % if there are multiple optimal parameters, pick the first one
            bestdistvector = parameter_space(1,maxid(1)).* (timevec) + parameter_space(2,maxid(1));
            Goodness_of_fit(event,traj) =  probmax./length(nonzerobins);
            bestparameters = parameter_space(:,maxid(1));%[speed or slope, interc]

            clear prob_sum

            % plot pMat and line fit
            if figopt
                subplot(1,track_num,traj)
                imagesc(1:szPM2,positions,pMat);
                colormap(colormapformulae([34,35,36])); % it is same as 'black-red-yellow-white'
                caxis([0 0.1])

                hold on
                plot(timevec,bestdistvector, 'cyan','linewidth',3,'linestyle',':');

                hold off
            end

            pMat_cell{event}{traj}.pMat = pMat;
            pMat_cell{event}{traj}.timevec = timevec;
            pMat_cell{event}{traj}.posvec = positions;
            pMat_cell{event}{traj}.timebinsz = tBinSz;
            pMat_cell{event}{traj}.spatialbinsz = sBinSz;
            pMat_cell{event}{traj}.bstline = bestdistvector;
            pMat_cell{event}{traj}.bstSpd = bestparameters(1);
            pMat_cell{event}{traj}.bstInterc = bestparameters(2);

            %% %-------Shuffling positions to get the pvalue for each traj------%
            permbins = nonzerobins;
            shiftVal = randi([1,szPM1-1],nperm,length(permbins));

            % initializing
            GoF_shuffle = [];
            srvalues = [];
            pmat0 =  zeros(size(ib));
            prob_sum0 = zeros(1,length(distvector(1,:)));

            % looping
            for iteration = 1:nperm

                % circularly shift position columns
                pMatShuf = zeros(size(pMat));=
                for t = 1:length(permbins)
                    pMatShuf(:,permbins(t)) = circshift(pMat(:,permbins(t)),shiftVal(iteration,t));
                end

                % calculate weighted correlation
                mloc = sum(sum(pMatShuf.*positions'))/sum(sum(pMatShuf));
                mt = sum(sum(pMatShuf.*timevec))/sum(sum(pMatShuf));
                dloc = positions'-mloc;
                dt = timevec - mt;
                cov_loc_t = sum(sum(pMatShuf.*dloc.*dt))/sum(sum(pMatShuf));
                cov_loc = sum(sum(pMatShuf.*(dloc.^2)))/sum(sum(pMatShuf));
                cov_t = sum(sum(pMatShuf.*(dt.^2)))/sum(sum(pMatShuf));
                srvalues = [srvalues;cov_loc_t/sqrt(cov_loc*cov_t)];

                prob_sum = prob_sum0;
                for t = 1:length(nonzerobins)
                    pMat_tmp = pMatShuf(:,nonzerobins(t));

                    ibic = ibic_total{t};
                    icc = icc_total{t};

                    pmat_total = pmat0;
                    pmat_total(icc) = pMat_tmp(ibic);

                    prob_sum = prob_sum + sum(pmat_total);

                end

                probmax = max(prob_sum);
                GoF_shuffle =  [GoF_shuffle;probmax./length(nonzerobins)];
            end
            if rvalues(event,traj) > 0
                pvalue_r(event,traj) =sum(rvalues(event,traj) < srvalues)/length(srvalues);
            else
                pvalue_r(event,traj) =sum(rvalues(event,traj) > srvalues)/length(srvalues);
            end
            pvalue(event,traj) =sum(Goodness_of_fit(event,traj) < GoF_shuffle)/length(GoF_shuffle);

            shuffle_GoF{event}{traj} = GoF_shuffle;
            shuffle_rvalues{event}{traj} = srvalues;
            clear ib_total
        end

        dtraj_id = find(pvalue(event,:) < 0.025 & pvalue_r(event,:) < 0.025);
        if ~isempty(dtraj_id)
            [~,tidx] = max(abs(rvalues(event,dtraj_id)));
            decode_traj(event) = dtraj_id(tidx);
        else
            decode_traj(event) = 0;% no significant traj
        end
        activecell{event} = cellsi;
        activecellidx{event} = cellidxm(cellsi,:);

        clear wrking factSpkPerBin expon

        title(sprintf('Decoded Traj =  %d;', decode_traj(event)));
        if figopt && (decode_traj(event) > 0)% only save the significant sequences
            figfile = [savedir,animalname,'_CA1_replayHSEseq','-Ep',num2str(ep),'_event',num2str(event)];
            print('-djpeg', figfile);
        end
        pause(0.05)
        clf
    end
end
%% save data into a structure
replaytraj.eventtimes = eventtimes;
replaytraj.pMat = pMat_cell;
replaytraj.Goodness_of_fit = Goodness_of_fit;
replaytraj.rvalues = rvalues;
replaytraj.eventidx = eventindex;
replaytraj.besttraj = decode_traj;
replaytraj.shuffle_GoF = shuffle_GoF;
replaytraj.shuffle_rvalues = shuffle_rvalues;
replaytraj.pvalue = pvalue;
replaytraj.pvalue_r = pvalue_r;
replaytraj.activecell = activecell;
replaytraj.activecellidx = activecellidx;
replaytraj.sigeventprc = length(find(decode_traj~=0))./length(decode_traj);
replaytraj.sigeventnum = length(find(decode_traj~=0));
replaytraj.candeventnum = length(decode_traj);
replaytraj.wellcutoff = wellcutoff;
replaytraj.cellcountthresh = cellcountthresh;
replaytraj.halfwidth = halfwidth;
replaytraj.wellcutoff = wellcutoff;

replaytrajactory{ep} = replaytraj;
%% save data
if savedata
    save(sprintf('%s%sreplayHSEseqdecode_CA1_%02d.mat', savedir,animalname,ep), 'replaytrajactory');
end
