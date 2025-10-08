
addpath(genpath('/Users/sld33/Dropbox/fb-hhmm/Matlab-preprocessing'))
addpath(genpath('/Users/sld33/Dropbox/animaltags/tagtools_matlab'), '-begin')
% don't add the directory "compatibility" as it is for users of older matlab versions
rmpath(genpath('/Users/sld33/Dropbox/animaltags/tagtools_matlab/compatibility'))
% don't add R crap or git crap
rmpath(genpath('/Users/sld33/Dropbox/animaltags/tagtools_matlab/.Rproj.user'))
rmpath(genpath('/Users/sld33/Dropbox/animaltags/tagtools_matlab/.git'))
savepath()

smrt_path = '/Users/sld33/Dropbox/FBdata/';

whales = {
    % 'Zica-20180113-173188';% Lander2
    % 'Zica-20180330-173188'; % Lander2
    % 'Zica-20180331-173187'; % Lander2
    % 'Zica-20190111-173186'; % Lander2
    % 'Zica-20190113-151361'; % SMRT (timing issue in processing)
    % 'Zica-20191012-144029'; % SMRT (timing issue in processing)
    % 'Zica-20191012-145101'; % SMRT (timing issue in processing)
    % 'Zica-20191111-94810';  % SMRT (timing issue in processing)
    % 'Zica-20191117-195993'; % SMRT (timing issue in processing)
    % 'Zica-20200106-195994'; % SMRT (timing issue in processing)
    % 'Zica-20220112-195994'; % newest SMRT
    % 'Zica-20211113-195993'; % newest SMRT
    % 'Zica-20211112-94819' ;
    'Zica-20211001-220819'; % Guadalupe 
    'Zica-20220930-195992'; % no mfas. Guadalupe
    'Zica-20221002-237889'; % Guadalupe
    'Zica-20230518-233391'; % has individual file
    'Zica-20230519-232950'; % has individual file
    'Zica-20230723-233394';
    'Zica-20230723-233395';
    'Zica-20231206-232951'; % Guadalupe
    'Zica-20240227-233396' ; 
    'Zica-20240227-240128' ;
    'Zica-20240526-233390' % Guadalupe
    } ;% newest SMRT

tag_types = [
    %1 1 1 1 2 2 2 2 2 2 3 3 3
    3 3 3 3 3 3 3 3 3 3 3]; 
% 1 = Lander2, 2 = SMRT (timing), 3 = SMRT (newest)

    
for w = 1:length(whales)
    disp(['Processing whale ', whales{w}])
    
    % Read in NC file
    %##########################################################################
    tag = whales{w};
    load_nc([smrt_path, tag, '-cal.nc'])
    if exist('Aw', 'var')
        Aa = Aw;
    end
    if exist('Mw', 'var')
        Ma = Mw;
    end
    if exist('depth_corrected', 'var')
        depth = depth_corrected;
    end
    if exist('CorrectedDepth', 'var')
        depth = CorrectedDepth;
    end
% code below will skip and Lander2 tags
%     if ~exist('Ma', 'var')
%         clearvars -except whales w
%         continue
%     end

Aa5 = decdc(Aa, floor(Aa.sampling_rate / 5)) ;
[pitch, roll] = a2pr(Aa5) ;

    % Compute fluke rate/amplitude using LMML method
    %##########################################################################
    
    % -----------------------------------------------------------------------
    % find dominant stroke freq
    % -----------------------------------------------------------------------
    ff = 0.7;
    fc = 0.4*ff;    % filter cut-off frequency in Hz - about 0.4 x fluking rate
    
    % -----------------------------------------------------------------------
    % Use magnetometer data to estimate body rotations 
    % -----------------------------------------------------------------------
    if exist('Ma', 'var')
        Ma5 = decdc(Ma, Ma.sampling_rate / 5) ;
        ph = get_br(Ma5, fc) ;
    
    % -----------------------------------------------------------------------
    % stroke detection
    % -----------------------------------------------------------------------
    
    % find zero-crossings
        thr = 2; % threshold in degrees
        [K,s] = zero_crossings(ph, thr*pi/180, Ma5.sampling_rate/0.2) ;
    
    % classify:  Up-strokes were defined as the period between a positive and a
    % negative zero-crossing and down-strokes from negative to a positive
    % zerocrossing. Intervals >2/fr between zero crossings were defined as
    % glides
    
    % stroke types: 1 = up, -1 = down, 0 = glide
        stroke_type = NaN*ones(length(ph),1);
        stroke_int = diff(K/Ma.sampling_rate); % time between detected strokes in seconds
    
        for i = 1:(length(K)-1)
            % indices of this period (from one detected zero-crossing to the next)
            ix = floor(K(i)):ceil(K(i+1));
            if stroke_int > 2/ff % glides are long intervals between zero-crossings
                stroke_type(ix) = 0;
            elseif stroke_int(i) <= 2/ff && s(i)*s(i+1) == -1 % strokes are times from -1 to 1 or 1 to -1 zero-crossings
                stroke_type(ix) = s(i);
            end
        end
    
        stroke_logical = zeros(length(ph), 1);
        stroke_logical(round(K),1) = 1;
    else
        stroke_logical = NaN*zeros(size(Aa5.data, 1), 1) ;
        stroke_type = stroke_logical ;
        stroke_int = stroke_logical ;
    end % end of "if Ma exists"
    
    % MSA
    MSA = msa(Aa);
    MSA5 = msa(Aa5);
    
    % Vertical (downward) velocity )
    % note + = increasing depth
    % and = - decreasing depth
    % Velocity (horizontal)
    % issue, interp2length not returning a struct?
    Aa_lo = decdc(Aa, Aa.sampling_rate / depth.sampling_rate);
    Aa_lo = interp2length(Aa_lo, depth);
    [s, vv] = speed_from_depth(depth.data, Aa_lo, depth.sampling_rate);
    s = abs(s); % speed is signed and mostly negative...just no
    
    
    % dive wiggliness
    
    % Windowed calculations
    %##########################################################################
    % duration of time-step (averaging window) in seconds
    win_sec = 5*60;
    % number of windows there will be for this tagout
    n_win = ceil((length(depth.data) / depth.sampling_rate) / win_sec);
    % times of window centers in second since start
    win_ctr_sec = [win_sec / 2 : win_sec : n_win * (win_sec)];
    
    % variance of pitch, roll, head
    pitch_buff = buffer(pitch, win_sec*Aa5.sampling_rate);
    roll_buff = buffer(roll, win_sec*Aa5.sampling_rate);
    if exist('Ma', 'var')
        [head, junk, junk2] = m2h(Ma5, Aa5) ;
        head_buff = buffer(head, win_sec*Aa5.sampling_rate);
    else
        head_buff = NaN * zeros(size(pitch_buff)) ;
    end
    
    [junk, std_pitch] = circ_std(pitch_buff);
    [junk, std_roll] = circ_std(roll_buff);
    [junk, std_head] = circ_std(head_buff);
    
    % mean MSA
    msa_buff = buffer(MSA, win_sec*Aa.sampling_rate);
    mean_msa = mean(msa_buff);

    msa_buff5 = buffer(MSA5, win_sec*Aa5.sampling_rate); % 5 Hz
    mean_msa5 = mean(msa_buff5); % 5 Hz
    
    % 90% MSA
    pct90_msa = quantile(msa_buff, 0.9);
    pct90_msa5 = quantile(msa_buff5, 0.9); % 5 Hz
    
    % fit EVD to -max(MSA each second), compute and report -mean of fitted distro
    % OR fit gEVD fitdist(data, 'GeneralizedExtremeValue') to max(MSA each second) using fitdist, then mean(result) to
    % report
    
    %allocate space
    evd_med_msa = NaN*ones(1, size(msa_buff, 2));
    evd_med_msa5 = NaN*ones(1, size(msa_buff5, 2));
    
    ops = statset('MaxIter', 2000, 'MaxFunEvals', 4000);
    
    for c = 1:size(msa_buff, 2)
        this_max_buff = buffer(msa_buff(:,c), 5*Aa.sampling_rate);
        maxes = max(this_max_buff);
        % NaN stuff needed for tags with gap in sensor data that is NaN
        % filled
        maxes(isnan(maxes)) = [];
        if (~isempty(maxes))
            warning('off','all')
            gevd_fit = fitdist(maxes(:), 'GeneralizedExtremeValue', 'options', ops);
            evd_med_msa(c) = median(gevd_fit);
            warning('on','all')
        else
            evd_med_msa(c) = NaN;
        end
        clear gevd_fit this_max_buff maxes
    end

    % at 5 Hz
    for c = 1:size(msa_buff5, 2)
        this_max_buff5 = buffer(msa_buff5(:,c), floor(5*Aa5.sampling_rate));
        maxes5 = max(this_max_buff5);
        % NaN stuff needed for tags with gap in sensor data that is NaN
        % filled
        maxes5(isnan(maxes5)) = [];
        if (~isempty(maxes5) && sum(size(maxes5)) > 30)
            warning('off','all')
            gevd_fit5 = fitdist(maxes5(:), 'GeneralizedExtremeValue', 'options', ops);
            evd_med_msa5(c) = median(gevd_fit5);
            warning('on','all')
        else
            evd_med_msa5(c) = NaN;
        end
        clear gevd_fit5 this_max_buff5 maxes5
    end
    
    
    % wiggliness (inflections / length of vector)
    pbuff = buffer([0; depth.data], win_sec*depth.sampling_rate);
    if (size(pbuff,2) > n_win)
        pbuff = pbuff(:, 1:n_win) ;
    end
        
    wiggles = sum(diff(sign(diff(pbuff))) ~= 0) ./ size(pbuff,1);
    
    % min and max depth
    min_depth = min(pbuff, [], 1, 'omitnan');
    max_depth = max(pbuff, [], 1, 'omitnan');
    med_depth = median(pbuff, 1, 'omitnan'); 
    
    % median speed
    speed_buff = buffer(s, win_sec*depth.sampling_rate);
    med_speed = median(speed_buff, 'omitnan');
    
    % median vertical velocity
    vv_buff = buffer(vv, win_sec*depth.sampling_rate);
    med_vv = median(vv_buff);
  
    if exist('Ma', 'var')
        % mean stroke rate
        sr_buff = buffer(stroke_logical, win_sec*Ma5.sampling_rate);
        strokes_per_sec = sum(sr_buff) / win_sec;
        % rms stroke amplitude
        sa_buff = buffer(ph(:,1), win_sec*Ma5.sampling_rate);
        rms_stroke = sqrt(mean(sa_buff.^2, 'omitnan'));
        peak_stroke = max(abs(sa_buff));
    else
        sa_buff = NaN * win_ctr_sec ;
        sr_buff = sa_buff;
        strokes_per_sec = sa_buff;
        rms_stroke = sa_buff;
        peak_stroke = sa_buff;
    end
    
    % Tagon time information (for converting with UTC times)
    % tag start time UTC as matlab datenum
    stimenum = datenum(info.dephist_device_datetime_start, ...
        info.dephist_device_regset);
    
    
    % ___________________________________________________________________
    % buzzes 
    % in comparison to before, for 2022, there are not "corr_times" files
    % and the files are in .csv rather than .xlsx format
    % ___________________________________________________________________
    % echolocation buzzes
    % Need to use BUZZ files and not EVENTS files for buzzes, per SC email
    buzzevent_infoA = dir([smrt_path 'ae-files/' tag '*PG_BUZZ_Events*.xlsx']);
    buzzevent_infoB = dir([smrt_path 'ae-files/' tag '*PG_BUZZ_Events*.csv']);
    buzzevent_info = [buzzevent_infoA; buzzevent_infoB];
    
    if (~isempty(buzzevent_info))
        for ckf = 1:length(buzzevent_info)
            if ckf == 1
                all_buzz_events = readtable([buzzevent_info(ckf).folder, '/', buzzevent_info(ckf).name]);
                % need to rename variables in "old corr" smrt files to
                % match 2022 ones
                if tag_types(w) == 2;
                    all_buzz_events = renamevars(all_buzz_events, ...
                        ["buzz_start_UTC_corrected","buzz_duration_s", ...
                        "note", "event_type_label"], ...
                        ["UTC", "Duration", "Note", "Label"]);
                end
                all_buzz_events = all_buzz_events(:, {'Note', 'Label', 'UTC', 'Duration'});
                % add sec_since_tagon variable
                if tag == 'Zica-20230723-233395'
                    all_buzz_events.sec_since_tagon = ...
                        (datenum(all_buzz_events.UTC, 'yyyy-mm-ddTHH:MM:SSZ') - stimenum)*24*60*60;
                else
                    all_buzz_events.sec_since_tagon = ...
                        (datenum(all_buzz_events.UTC) - stimenum)*24*60*60;
                end

            else
                more_buzz_events = readtable([buzzevent_info(ckf).folder, '/', buzzevent_info(ckf).name]);
                if tag_types(w) == 2;
                    more_buzz_events = renamevars(more_buzz_events, ...
                        ["buzz_start_UTC_corrected","buzz_duration_s", ...
                        "note", "event_type_label"], ...
                        ["UTC", "Duration", "Note", "Label"]);
                end
                more_buzz_events = more_buzz_events(:, {'Note', 'Label', 'UTC', 'Duration'});
                % add sec_since_tagon variable
                more_buzz_events.sec_since_tagon = ...
                    (datenum(more_buzz_events.UTC) - stimenum)*24*60*60;
                all_buzz_events = [all_buzz_events; more_buzz_events];
            end
        end
        
        % per SC keep ones that are labelled "BW Buzz" but not "Poss" or
        % "Possible"
        focal_buzzs = all_buzz_events(startsWith(all_buzz_events.Label, 'BW Buzz') | ...
            startsWith(all_buzz_events.Note, 'BW Buzz'), :);
        focal_buzzs = focal_buzzs(~contains(focal_buzzs.Label, 'Poss', 'IgnoreCase', true) & ...
            ~contains(focal_buzzs.Note, 'Poss', 'IgnoreCase', true), :);
        
        % convert all var names to lower case
        focal_buzzs.Properties.VariableNames = lower(focal_buzzs.Properties.VariableNames);
        % so names are: note, label, utc, duration, sec_since_tagon
    end
    % 
    
    % ___________________________________________________________________
    % clicks by tagged whale
    % need to use "post processed events" files for 2022 since audits were
    % not done at the level of logging every click so there are no
    % "allclicks" files
    % also in comparison to before, files are .csv and not .xlsx
    % ___________________________________________________________________
    allclick_infoA = dir([smrt_path 'ae-files/' tag '*PG_Post_Processed_Events*.xlsx']);
    allclick_infoB = dir([smrt_path 'ae-files/' tag '*PG_Post_Processed_Events*.csv']);
    allclick_info = [allclick_infoA; allclick_infoB];
    
    if (~isempty(allclick_info))
        for ckf = 1:length(allclick_info)
            if ckf == 1
                if tag == 'Zica-20211001-220819'
                    opts = detectImportOptions([allclick_info(ckf).folder, '/', allclick_info(ckf).name]);
                    opts = setvaropts(opts, 'UTC','InputFormat','MM/dd/uuuu HH:mm:ss.SSS');
                    all_clicks = readtable([allclick_info(ckf).folder, '/', allclick_info(ckf).name], opts);
                else
                    all_clicks = readtable([allclick_info(ckf).folder, '/', allclick_info(ckf).name]);
                end
                if tag_types(w) == 2
                   all_clicks = renamevars(all_clicks, ...
                        ["event_UID", "event_start_UTC_corrected", "event_end_UTC_corrected"], ...
                        ["UID", "UTC", "event_end"]);
                else
                   all_clicks = renamevars(all_clicks, ...
                       ["eventType", "EventEnd"], ["event_type", "event_end"]);
                end 
                all_clicks = all_clicks(:, {'UID', 'UTC', 'event_type', 'event_end'}) ;
                   % these vars no longer available in later smrt
                   % deployments:
                   % 'sec_since_file_start', 'samp_since_file_start', 'wav_file', ...
               % add sec_since_tagon variable
                all_clicks.sec_since_tagon = ...
                    (datenum(all_clicks.UTC) - stimenum)*24*60*60;
               % add duration in seconds
               all_clicks.duration = etime(datevec(datenum(all_clicks.event_end)), ...
                   datevec(datenum(all_clicks.UTC)));
            else
                more_clicks = readtable([allclick_info(ckf).folder, '/', allclick_info(ckf).name]);
                if tag_types(w) == 2
                   more_clicks = renamevars(more_clicks, ...
                        ["event_UID","event_start_UTC_corrected", "event_end_UTC_corrected"], ...
                        ["UID", "UTC", "event_end"]);
                else
                   more_clicks = renamevars(more_clicks, ...
                       ['eventType', 'EventEnd'], ['event_type', 'event_end']);
                end 
                more_clicks = more_clicks(:, {'UID', 'UTC', 'event_type', 'event_end'});
               % add sec_since_tagon variable
                more_clicks.sec_since_tagon = ...
                    (datenum(more_clicks.UTC) - stimenum) * 24*60*60;
                % add duration in seconds
               more_clicks.duration = etime(datevec(datenum(more_clicks.event_end)), ...
                   datevec(datenum(more_clicks.UTC)));
                all_clicks = [all_clicks; more_clicks];
            end
            
        end
        focal_clicks = all_clicks(strcmp(all_clicks.event_type, 'FD'),:);
        % convert all var names to lower case
        focal_clicks.Properties.VariableNames = lower(focal_clicks.Properties.VariableNames);
        % so names are: uid, utc, event_type, sec_since_tagon, event_end,
        % duration
        
        % ___________________________________________________________________
        % presence/absence of tagged, conspecific clicks?
        % ___________________________________________________________________
        nf_clicks = all_clicks(strcmp(all_clicks.event_type, 'Other BW'), :);
    end
    
    % ___________________________________________________________________
    % count click related events in bins
    % we CANNOT count clicks (since not individually marked)
    % we CAN still count buzzes
    % we CAN still mark presence/absence of clicks (ignoring pauses)
    % in previous tags we had times of each click. now we have start and
    % dur for each clicking period.
    % ___________________________________________________________________
    
    if (exist('focal_clicks', 'var'))
        % allocate space
        n_focal_buzzes = zeros(length(win_ctr_sec),1) ;
        clicks_present = n_focal_buzzes ;
        clicks_n = NaN * n_focal_buzzes ;
        nonfocal_clicks_present = n_focal_buzzes ;
        nonfocal_clicks_n = NaN * n_focal_buzzes ;
        
        for b = 1:length(win_ctr_sec)
            clear fci nfci bi
            % find clicking periods that start before the end of the window
            % AND end after the start of the the window
            fci = focal_clicks.sec_since_tagon < (win_ctr_sec(b) + win_sec/2) & ...
                focal_clicks.sec_since_tagon + focal_clicks.duration >= (win_ctr_sec(b) - win_sec/2) ;
            nfci = nf_clicks.sec_since_tagon < (win_ctr_sec(b) + win_sec/2) & ...
                nf_clicks.sec_since_tagon + nf_clicks.duration >= (win_ctr_sec(b) - win_sec/2) ;
            % logical index vector to buzzes in this time window
            bi = focal_buzzs.sec_since_tagon < (win_ctr_sec(b) + win_sec/2) & ...
                focal_buzzs.sec_since_tagon >= (win_ctr_sec(b) - win_sec/2) ;
            n_focal_buzzes(b) = sum(bi) ;
            % need to check whether this window b is in a clicking bout
            if sum(fci, 'omitnan') > 0
                clicks_present(b) = 1 ;
            end
            % same for nonfocal clicks
            if sum(nfci, 'omitnan') > 0
                nonfocal_clicks_present(b) = 1 ;
            end
        end
        
    else
        n_focal_buzzes = NaN * zeros(length(win_ctr_sec),1) ;
        clicks_present = n_focal_buzzes ;
        clicks_n = n_focal_buzzes ;
        nonfocal_clicks_present = n_focal_buzzes ;
        nonfocal_clicks_n = n_focal_buzzes ;
    end
    
    % Save data for this whale
    %##########################################################################
    whale_ID = repmat(tag, length(peak_stroke), 1);
    win_ctr_sec = win_ctr_sec(:);
    std_pitch = std_pitch(:);
    std_roll = std_roll(:);
    std_head = std_head(:);
    mean_msa = mean_msa(:);
    mean_msa5 = mean_msa5(:);
    pct90_msa = pct90_msa(:);
    pct90_msa5 = pct90_msa5(:);
    evd_med_msa = evd_med_msa(:);
    evd_med_msa5 = evd_med_msa5(:);
    wiggles = wiggles(:);
    strokes_per_sec = strokes_per_sec(:);
    rms_stroke = rms_stroke(:);
    peak_stroke = peak_stroke(:);
    med_speed = med_speed(:);
    med_vv = med_vv(:);
    n_focal_buzzes = n_focal_buzzes(:) ;
    clicks_present = clicks_present(:) ;
    % clicks_n = clicks_n(:) ;
    nonfocal_clicks_present = nonfocal_clicks_present(:) ;
    % nonfocal_clicks_n = nonfocal_clicks_n(:) ;
    min_depth = min_depth(:);
    max_depth = max_depth(:);
    med_depth = med_depth(:);
    
    
    
    this_whale_dat = table(whale_ID, win_ctr_sec, ...
        std_pitch, std_roll, std_head, ...
        mean_msa, pct90_msa, evd_med_msa, ...
        mean_msa5, pct90_msa5, evd_med_msa5, ...
        wiggles, med_speed, med_vv, ...
        strokes_per_sec, rms_stroke, peak_stroke, ...
        n_focal_buzzes, clicks_present, ... % clicks_n, ...
        nonfocal_clicks_present, ... % nonfocal_clicks_n,
        min_depth, max_depth, med_depth);
    
    writetable(this_whale_dat, [smrt_path, 'preproc-fine-coarse/', tag, '-finescale.csv']);
    clearvars -except whales w smrt_path tag_types
end