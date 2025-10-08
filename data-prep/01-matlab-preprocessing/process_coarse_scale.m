% matlab paths and file paths
%*************************************************************************
addpath(genpath('/Users/sld33/Dropbox/fb-hhmm/Matlab-preprocessing'))
addpath(genpath('/Users/sld33/Dropbox/animaltags/tagtools_matlab'), '-begin')
% don't add the directory "compatibility" as it is for users of older matlab versions
rmpath(genpath('/Users/sld33/Dropbox/animaltags/tagtools_matlab/compatibility'))
% don't add R crap or git crap
rmpath(genpath('/Users/sld33/Dropbox/animaltags/tagtools_matlab/.Rproj.user'))
rmpath(genpath('/Users/sld33/Dropbox/animaltags/tagtools_matlab/.git'))
savepath()

smrt_path = '/Users/sld33/Dropbox/FBdata/';
%*************************************************************************

% output from DAS model for foraging predictions
%*************************************************************************
foraging_model_preds = readtable(strcat(smrt_path, "MLdiveclassifications_20210810.csv"));
foraging_model_preds = foraging_model_preds(:, ["TagID", "Start", "End", "Foraging", "preds_class_bottom_tree"]);
foraging_model_preds = renamevars(foraging_model_preds, "preds_class_bottom_tree", "ForagingModPrediction");

foraging_model_preds22 = readtable(strcat(smrt_path, "AllSMRTDives_ForagingPredicted_2022.csv"));
foraging_model_preds22 = foraging_model_preds22(:, ["TagID", "Start", "End", "Foraging", "ForagingModPrediction"]);

% combine "lander2 and old smrt" output with "all smrt" output, remove dups
FMpred = [foraging_model_preds; foraging_model_preds22];
FMpred = unique(FMpred,'rows');

% add variable to table that is "Foraging" if known from acoustics, 
% and "ForagingModPrediction" if not
% this will be the dive class for determining foraging dive cycles
HHMMForageClass = FMpred.Foraging;
unk_dives = find(HHMMForageClass == "Unknown");
HHMMForageClass(unk_dives) = FMpred.ForagingModPrediction(unk_dives);
FMpred = addvars(FMpred, HHMMForageClass);

% there is a dive for whale Zica-20190111-173186 with NaN for most data
% and no acous so the foraging status is NaN. But need to treat it as a
% foraging cycle or else its time will be added on to the previous one!
FMpred_mini = cell2table({'Zica-20190111-173186',  474441, 478928, ...
    'Unknown', 'Unknown', 'Yes'});
FMpred_mini.Properties.VariableNames = FMpred.Properties.VariableNames;
FMpred = [FMpred; FMpred_mini];

% arrange by TagID then Start time
FMpred = sortrows(FMpred, {'TagID', 'Start'});

%*************************************************************************

% info about tag ID strings and tag types
%*************************************************************************
whales = {
    % 'Zica-20180113-173188';% Lander2
    % 'Zica-20180330-173188'; % Lander2
    % 'Zica-20180331-173187'; % Lander2
    % 'Zica-20190111-173186'; % Lander2
    'Zica-20190113-151361'; % SMRT (timing issue in processing)
    'Zica-20191012-144029'; % SMRT (timing issue in processing)
    'Zica-20191012-145101'; % SMRT (timing issue in processing)
    'Zica-20191111-94810';  % SMRT (timing issue in processing)
    'Zica-20191117-195993'; % SMRT (timing issue in processing)
    % skip this whale: has only one deep dive, no complete dive cycles
    % 'Zica-20200106-195994'; % SMRT (timing issue in processing)
    'Zica-20220112-195994'; % newest SMRT
    'Zica-20211113-195993'; % newest SMRT
    'Zica-20211112-94819';
    'Zica-20220112-195994'; % has individual file
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
    } ;

tag_types = [% 1 1 1 1 
    2 2 2 2 2 ... % 2 (Zica-20200106-195994)
    3 3 3 ...
    3 3 3 3 3 3 3 3 3 3 3 3]; 
% 1 = Lander2, 2 = SMRT (timing issues), 3 = SMRT (newest)
%*************************************************************************

% loop over whales and process coarse scale data
%*************************************************************************
for w = 1:length(whales)
    clear tinfo
    disp(['Processing whale ', whales{w}])
    
    % Read in NC file
    %*************************************************************************
    tag = whales{w};
    load_nc([smrt_path, tag, '-cal']);
    
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

    if (exist('wet', 'var'))
        if max(wet.data) > 80 % lander2 wet data goes up to 255 - it's not threshold-detected, it's raw
            wet.data = wet.data >= 80; % 80 is the threshold used in SMRT tags
            wet.conductivity_threshold = 80;
        end
    end
    
    % if ~exist('Ma', 'var')
    %     clearvars -except whales w
    %     continue
    % end
    
    % Tagon time information (for converting with UTC times)
    % tag start time UTC as matlab datenum
    stimenum = datenum(info.dephist_device_datetime_start, info.dephist_device_regset);
    
    
    % Click audit/detector data
    %*************************************************************************
    % read in events:
    % note some are excel and some are csv
    clickevent_infoA = dir([smrt_path 'ae-files/' tag '*Post_Processed_Events*.xlsx']);
    clickevent_infoB = dir([smrt_path 'ae-files/' tag '*Post_Processed_Events*.csv']);
    clickevent_info = [clickevent_infoA; clickevent_infoB];
    
    if (~isempty(clickevent_info))
        for ckf = 1:length(clickevent_info)
            if ckf == 1
                if strcmp(tag,'Zica-20211001-220819')
                    opts = detectImportOptions([clickevent_info(ckf).folder, '/', clickevent_info(ckf).name]);
                    opts = setvaropts(opts, 'UTC','InputFormat','MM/dd/uuuu HH:mm:ss.SSS');
                    all_ck_events = readtable([clickevent_info(ckf).folder, '/', clickevent_info(ckf).name], opts);
                else
                    all_ck_events = readtable([clickevent_info(ckf).folder, '/', clickevent_info(ckf).name]);
                end
                if tag_types(w) == 2
                   all_ck_events = renamevars(all_ck_events, ...
                        ["event_UID", "event_start_UTC_corrected", "event_end_UTC_corrected"], ...
                        ["UID", "UTC", "event_end"]);
                else
                   all_ck_events = renamevars(all_ck_events, ...
                       ["eventType", "EventEnd"], ["event_type", "event_end"]);
                end 
                if tag_types(w) ~= 2
                    all_ck_events = all_ck_events(:, {'UID', 'UTC', 'event_type', 'event_end'}) ;
                else
                    all_ck_events = all_ck_events(:, {'UID', 'UTC', 'event_type', 'event_end', 'breath_count'});
                end
            else
                more_clicks = readtable([clickevent_info(ckf).folder, '/', clickevent_info(ckf).name]);
                if tag_types(w) == 2
                   more_clicks = renamevars(more_clicks, ...
                        ["event_UID", "event_start_UTC_corrected", "event_end_UTC_corrected"], ...
                        ["UID", "UTC", "event_end"]);
                   more_clicks = more_clicks(:, {'UID', 'UTC', 'event_type', 'event_end', 'breath_count'}) ;
                else
                   more_clicks = renamevars(more_clicks, ...
                       ["eventType", "EventEnd"], ["event_type", "event_end"]);
                   more_clicks = more_clicks(:, {'UID', 'UTC', 'event_type', 'event_end'}) ;
                end 
                
                all_ck_events = [all_ck_events; more_clicks];
            end % end of "if file 2+"  
        end % end of loop over files
    
    % add sec_since_tagon variable
    all_ck_events.start_sec_since_tagon = ...
        (datenum(all_ck_events.UTC) - stimenum) * 24*60*60;
    % add duration in seconds
    all_ck_events.duration = etime(datevec(datenum(all_ck_events.event_end)), ...
        datevec(datenum(all_ck_events.UTC)));
    % and end time in seconds
    all_ck_events.end_sec_since_tagon = all_ck_events.start_sec_since_tagon + ...
        all_ck_events.duration ;
    
    focal_clicks = all_ck_events(strcmp(all_ck_events.event_type, 'FD'),:);
    focal_clicks.Properties.VariableNames = lower(focal_clicks.Properties.VariableNames);
        % so names are: uid, utc, event_type, sec_since_tagon, event_end,
        % duration
    end % end of "if there are click files
        
    
    % echolocation buzzes
    %*************************************************************************
    % Need to use BUZZ files and not EVENTS files for buzzes, see SC email
    % also: argh, some csv and some xlsx
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
            else
                more_buzz_events = readtable([buzzevent_info(ckf).folder, '/', buzzevent_info(ckf).name]);
                if tag_types(w) == 2
                    more_buzz_events = renamevars(more_buzz_events, ...
                        ["buzz_start_UTC_corrected","buzz_duration_s", ...
                        "note", "event_type_label"], ...
                        ["UTC", "Duration", "Note", "Label"]);
                end
                more_buzz_events = more_buzz_events(:, {'Note', 'Label', 'UTC', 'Duration'});
                all_buzz_events = [all_buzz_events; more_buzz_events];
            end
            % add sec_since_tagon variable
            if strcmp(tag, 'Zica-20230723-233395')
                all_buzz_events.sec_since_tagon = ...
                    (datenum(all_buzz_events.UTC, 'yyyy-mm-ddTHH:MM:SSZ') - stimenum)*24*60*60;
            else
                all_buzz_events.sec_since_tagon = ...
                    (datenum(all_buzz_events.UTC) - stimenum)*24*60*60;
            end
        end
        
        % per SC keep ones that are labelled "BW Buzz" but not "Poss" or
        % "Possible"
        focal_buzzs = all_buzz_events(startsWith(all_buzz_events.Label, 'BW Buzz') | ...
            startsWith(all_buzz_events.Note, 'BW Buzz'), :);
        focal_buzzs = focal_buzzs(~contains(focal_buzzs.Label, 'Poss', 'IgnoreCase', true) & ...
            ~contains(focal_buzzs.Note, 'Poss', 'IgnoreCase', true), :);
        % sec_since_tagon, Duration
        % convert all var names to lower case
        focal_buzzs.Properties.VariableNames = lower(focal_buzzs.Properties.VariableNames);
        % so names are: note, label, utc, duration, sec_since_tagon
    
        nf_clicks = all_ck_events(strcmp(all_ck_events.event_type, 'Other BW'), :);
        if any("breath_count" == string(all_ck_events.Properties.VariableNames))
            man_breaths = all_ck_events(strcmp(all_ck_events.event_type, 'Surf'), :);
        end

    end % end "if there are buzz files"
    
    % Calculate dive statistics
    %*************************************************************************
    % Detect dives
    D = find_dives(depth, 50);
    % Dive durations in seconds
    D.dur = D.end - D.start;
    D.post_surf_dur = [D.start(2:end) - D.end(1:(end-1)); NaN];
    
    
    % Dive depths
    % are in D.max
    
    % Classify foraging/not
    %*************************************************************************
    % using data from DAS model
    % not applied to later tags so...
    % this_FMpred = FMpred(contains(FMpred.TagID,tag), :);
    % Dtab = struct2table(D);
    % Dtab = renamevars(Dtab, ["start", "end"], ["Start", "End"]);
    % this_FMpredD = join(Dtab, this_FMpred);

    % old way using K-means here
    % [idx, C] = kmeans([D.max, D.dur], 2);
    % make sure 1 is always deep and 2 is shallow
    % [maxs, imax] = max(C);

    % values will be "deep" and "shal"
    % numeric version, 2 is "deep" and 1 is "shal"
    % use this_FMpredD.HHMMForageClass
    D.dtype = repmat("not modeled", length(D.max), 1 );
    D.dtype_num = ones(size(D.dtype)); % "shal"
    % D.dtype(contains(this_FMpredD.HHMMForageClass, 'Yes')) = "deep";
    % D.dtype_num(contains(this_FMpredD.HHMMForageClass, 'Yes')) = 2;
    
    %*************************************************************************
    
    % more dive stats
    %*************************************************************************
    % allocate space for audit/click related vars
    D.click_dur = NaN * zeros(length(D.dur), 1); 
    D.click_end_sec = NaN * zeros(length(D.dur), 1);
    D.click_start_sec = NaN * zeros(length(D.dur), 1); 
%     D.num_clicks = NaN * zeros(length(D.dur), 1); % didn't detect w/later tags
    D.click_start_depth = NaN * zeros(length(D.dur), 1);
    D.click_end_depth = NaN * zeros(length(D.dur), 1);
    D.click_max_depth = NaN * zeros(length(D.dur), 1);
    D.click_min_depth = NaN * zeros(length(D.dur), 1);
    % non focal clicks
    D.nf_click_dur = NaN * zeros(length(D.dur), 1); 
    D.nf_click_start_depth = NaN * zeros(length(D.dur), 1);
    D.nf_click_end_depth = NaN * zeros(length(D.dur), 1);
    % buzzes
    D.num_buzzes = NaN * zeros(length(D.dur), 1); 
    D.buzz_dur_total = NaN * zeros(length(D.dur), 1); 
    D.buzz_dur_mean = NaN * zeros(length(D.dur), 1); 
    D.buzz_dur_median = NaN * zeros(length(D.dur), 1); 
    D.buzz_dur_q25 = NaN * zeros(length(D.dur), 1); 
    D.buzz_dur_q75 = NaN * zeros(length(D.dur), 1); 
    % breath_count with event_type == 'Surf' in this_events
    D.breath_count_manual = NaN * zeros(length(D.dur), 1); 
    D.buzz_depth_median = NaN * zeros(length(D.dur), 1); 
    D.buzz_depth_mean = NaN * zeros(length(D.dur), 1); 
    D.buzz_depth_max = NaN * zeros(length(D.dur), 1); 
    D.buzz_depth_min = NaN * zeros(length(D.dur), 1); 
    D.buzz_depth_q25 = NaN * zeros(length(D.dur), 1); 
    D.buzz_depth_q75 = NaN * zeros(length(D.dur), 1); 
    
    % determine if foraging present for each dive
    if (exist("focal_clicks", "var"))
    for d = 1:length(D.dur)
        % find any foraging click periods that are start during this dive
        % record total duration and total number of clicks
        ck_id = focal_clicks.start_sec_since_tagon > D.start(d) & focal_clicks.start_sec_since_tagon < D.end(d);
        nf_ck_id = nf_clicks.start_sec_since_tagon > D.start(d) & nf_clicks.start_sec_since_tagon < D.end(d);
        
        if sum(ck_id) > 0
            D.click_dur(d) = sum(focal_clicks.end_sec_since_tagon(ck_id) - ...
                focal_clicks.start_sec_since_tagon(ck_id));
            D.click_end_sec(d) = max(focal_clicks.end_sec_since_tagon(ck_id));
            D.click_start_sec(d) = min(focal_clicks.start_sec_since_tagon(ck_id));
            D.click_start_depth(d) = depth.data(round(min(focal_clicks.start_sec_since_tagon(ck_id)) * ...
                depth.sampling_rate)) ;
            D.click_end_depth(d) = depth.data(round(max(focal_clicks.end_sec_since_tagon(ck_id)) * ...
                depth.sampling_rate)) ;
            this_focal_ix = round(min(focal_clicks.start_sec_since_tagon(ck_id)) * depth.sampling_rate) : ...
                round(max(focal_clicks.end_sec_since_tagon(ck_id)) * depth.sampling_rate);
            D.click_max_depth(d) = max(depth.data(this_focal_ix));
            D.click_min_depth(d) = min(depth.data(this_focal_ix));
            % non focal clicks
            D.nf_click_dur(d) = sum(nf_clicks.end_sec_since_tagon(nf_ck_id) - ...
                nf_clicks.start_sec_since_tagon(nf_ck_id));
            if D.nf_click_dur(d) > 0
                D.nf_click_start_depth(d) = depth.data(round(min(nf_clicks.start_sec_since_tagon(nf_ck_id)) * ...
                depth.sampling_rate)) ;
                D.nf_click_end_depth(d) = depth.data(round(max(nf_clicks.end_sec_since_tagon(nf_ck_id)) * ...
                depth.sampling_rate)) ;
            end
            % D.num_clicks(d) = sum(focal_clicks.nClicks(ck_id));
            % find any buzzes that are in this dive
            bz_id = focal_buzzs.sec_since_tagon > D.start(d) & focal_buzzs.sec_since_tagon < D.end(d);
            bz_z = depth.data(round(depth.sampling_rate * focal_buzzs.sec_since_tagon(bz_id))); % depth @ start of buzz
            D.num_buzzes(d) = sum(bz_id, 'omitnan');
            D.buzz_dur_total(d) = sum(focal_buzzs.duration(bz_id), 'omitnan'); % in seconds
            D.buzz_dur_mean(d) = mean(focal_buzzs.duration(bz_id), 'omitnan');
            D.buzz_dur_median(d) = median(focal_buzzs.duration(bz_id), 'omitnan');
            D.buzz_dur_q25(d) = quantile(focal_buzzs.duration(bz_id), 0.25);
            D.buzz_dur_q75(d) = quantile(focal_buzzs.duration(bz_id), 0.75);
            D.buzz_depth_median(d) = median(bz_z, 'omitnan');
            D.buzz_depth_mean(d) = mean(bz_z, 'omitnan');
            D.buzz_depth_max(d) = max([NaN; bz_z], [], 'all', 'omitnan');
            D.buzz_depth_min(d) = min([NaN; bz_z], [], 'all', 'omitnan');
            D.buzz_depth_q25(d) = quantile(bz_z, 0.25);
            D.buzz_depth_q75(d) = quantile(bz_z, 0.75);
            if exist('man_breaths', 'var')
                % manually counted breaths
                if d < length(D.dur)
                    br_id = man_breaths.start_sec_since_tagon > D.start(d) & man_breaths.start_sec_since_tagon < D.start(d+1);
                    D.breath_count_manual(d) = sum(br_id, 'omitnan');
                end
            end
        end % end "if there are clicks in the dive     
    end % end loop over dives
% if there are no click audit info available they will stay NaN
    end
    
    % Detect IFDIs (time steps) -- one row per IFDI
    % instead of using D.click_dur > 0, need to use D.dtype_num == 2
    % to make sure it works for Lander2 also
    D.dtype_num(D.click_dur > 0) = 2;
    nDC = sum(D.dtype_num == 2);
    var_names = {'start_sec', 'end_sec', 'dur_sec', ...
                         'fd_dur_sec', 'nshal', 'nbreath', 'surf_sec', ...
                         'click_dur', ... % 'num_clicks', ...
                         'click_start_sec', 'click_end_sec', ...
                         'click_start_depth', 'click_end_depth', ...
                         'click_max_depth', 'click_min_depth', ...
                         'num_buzzes', 'buzz_dur_total', 'buzz_dur_mean', ...
                         'buzz_dur_median', 'buzz_dur_q25', 'buzz_dur_q75', ...
                         'buzz_depth_median', 'buzz_depth_max', 'buzz_depth_min', ...
                         'buzz_depth_q25', 'buzz_depth_q75', ...
                         'breath_count_manual', ...
                         'nf_click_bouts', 'nf_click_dur', ...
                         'nf_click_start_depth', 'nf_click_end_depth', ...
                         'path', 'step', 'tort', ...
                         'longest_dur_sec', 'longest_nf_dur_sec'};
    nVars = length(var_names);
    DC = table('Size',[nDC, nVars], ...
        'VariableTypes',cellstr(repmat('double', nVars, 1)), ...
        'VariableName', var_names);
    
    % note these will all be NA for dives w/no clicks AND NO ACOUS
    DC.start_sec = D.start(D.dtype_num == 2);
    DC.max_depth = D.max(D.dtype_num == 2);
    DC.end_sec = [DC.start_sec(2:end); NaN];
    DC.dur_sec = DC.end_sec - DC.start_sec;
    DC.fd_dur_sec = D.dur(D.dtype_num == 2);
    % echolocation (tag whale)
    DC.click_dur = D.click_dur(D.dtype_num == 2);
    DC.click_start_sec = D.click_start_sec(D.dtype_num == 2);
    DC.click_end_sec = D.click_end_sec(D.dtype_num == 2);
    DC.click_start_depth = D.click_start_depth(D.dtype_num == 2);
    DC.click_end_depth = D.click_end_depth(D.dtype_num == 2);
    DC.click_max_depth = D.click_max_depth(D.dtype_num == 2);
    DC.click_min_depth = D.click_max_depth(D.dtype_num == 2);
    % non focal clicks
    DC.nf_click_start_depth = D.nf_click_start_depth(D.dtype_num == 2);
    DC.nf_click_end_depth = D.nf_click_end_depth(D.dtype_num == 2);
    % buzzes
    DC.num_buzzes = D.num_buzzes(D.dtype_num == 2);
    DC.buzz_dur_total = D.buzz_dur_total(D.dtype_num == 2);
    DC.buzz_dur_mean = D.buzz_dur_mean(D.dtype_num == 2);
    DC.buzz_dur_median = D.buzz_dur_median(D.dtype_num == 2);
    DC.buzz_dur_q25 = D.buzz_dur_q25(D.dtype_num == 2);
    DC.buzz_dur_q75 = D.buzz_dur_q75(D.dtype_num == 2);
    DC.buzz_depth_median = D.buzz_depth_median(D.dtype_num == 2);
    DC.buzz_depth_mean = D.buzz_depth_mean(D.dtype_num == 2); 
    DC.buzz_depth_max = D.buzz_depth_max(D.dtype_num == 2);
    DC.buzz_depth_min = D.buzz_depth_min(D.dtype_num == 2);
    DC.buzz_depth_q25 = D.buzz_depth_q25(D.dtype_num == 2);
    DC.buzz_depth_q75 = D.buzz_depth_q75(D.dtype_num == 2);
    % breath_count with event_type == 'Surf' in this_events
    if exist('man_breaths', 'var')
    for dc = 1:size(DC,1)
        br_id = man_breaths.start_sec_since_tagon > DC.start_sec(dc) & man_breaths.start_sec_since_tagon < DC.end_sec(dc);
        DC.breath_count_manual(dc) = sum(man_breaths.breath_count(br_id), 'omitnan');
    end
    else
        DC.breath_count_manual = NaN * DC.start_sec;
    end
    
    
    Aa_lo = decdc(Aa, Aa.sampling_rate / depth.sampling_rate);
    Aa_lo = interp2length(Aa_lo, depth);
    Aa_lo(isnan(Aa_lo)) = 0.05*rand([sum(sum(isnan(Aa_lo))),1]);
    [s, vv] = speed_from_depth(depth.data, Aa_lo, depth.sampling_rate);
    x = 1:length(s) ;
    s(isnan(s)) = interp1(x(~isnan(s)),s(~isnan(s)),x(isnan(s)), 'linear', 'extrap') ;

    tort_dur = 900;
    if exist("Ma", "var")
        Ma_lo = decdc(Ma, Ma.sampling_rate / depth.sampling_rate);
        Ma_lo = interp2length(Ma_lo, depth);
        Ma_lo(isnan(Ma_lo)) = 0.05*rand([sum(sum(isnan(Ma_lo))),1]);
        [this_ptrack, pe] = ptrack(Aa_lo, Ma_lo, s, depth.sampling_rate);

        tort = tortuosity(this_ptrack, depth.sampling_rate, tort_dur);
        tort = tort(:,1);
    else
        this_ptrack = NaN * Aa_lo(:, 1:2);
        tort = NaN * Aa_lo(:,1);
    end
    
    % midpoints in seconds for tortuosity calc
    tort_sec = tort_dur/2 + [0:(-1 + length(tort))]*tort_dur;
    
    %add criterion for timing error? (emailed MarEcoTel to ask what rule to use on 7/14/2023)
    % per initial DAS reply, abs(error) <= 3 seconds. confirmed by GSS
    hi_gps = GPS_position.data(GPS_satellites.data(:,2) >= 4 & ...
                                ~isnan(GPS_satellites.data(:,2)) & ...
                                ~isnan(GPS_residual.data(:,2)) & ...
                                GPS_residual.data(:,2) < 35 & ...
                                abs(GPS_time_error.data(:,2)) <= 3 , :);
    loc_gps = lalo2llf(hi_gps(:, 2:3));

    [combo_track, current] = fit_tracks(loc_gps, hi_gps(:,1), ...
        [this_ptrack(:,2), this_ptrack(:,1)], depth.sampling_rate);
    
    for dc = 1:height(DC)
        c_start = DC.start_sec(dc);
        c_end = DC.end_sec(dc);
        if isnan(c_end)
            % for last dive cycle end time is unknown, use all the data
            % there is
            c_end = length(depth.data) / depth.sampling_rate;
        end
        % Number of surface dives
        % using D.dtype_num == 1 is shallow
        DC.nshal(dc) = sum(D.start >= c_start & D.start < c_end & D.dtype_num == 1 );
        DC.longest_dur_sec(dc) = max(D.dur(D.start >= c_start & D.start < c_end));
        DC.longest_nf_dur_sec(dc) = max([0; D.dur(D.start >= c_start & D.start < c_end & D.dtype_num == 1)]);
        % Number of breaths/surfacings
        this_depth = depth.data(round(depth.sampling_rate * c_start) : round(depth.sampling_rate * c_end));
        this_wet = wet.data(round(wet.sampling_rate * c_start) : round(wet.sampling_rate * c_end));
        % find where wet/dry sensor value crosses 0.5
        [K, s, KK] = zero_crossings(this_wet - 0.5, 0.1);
        % s is sign of crossing; find times in seconds of all the "toward
        % dry" crossings
         DC.nbreath(dc) = sum(s == 1);
         
         % surfacing - total seconds at less than 5m depth
        DC.surf_sec(dc) = sum(this_depth < 5, "omitnan") / depth.sampling_rate;
         
        
        % conspecific clicking (number of events; and duration)
        % need to make sure is NaN and not 0 for times after end of acous
        % so if click_dur is NaN this should be also
        if exist("nf_clicks", "var")
        nfci = nf_clicks.start_sec_since_tagon > c_start & ...
            nf_clicks.start_sec_since_tagon < c_end ;
        DC.nf_click_bouts(dc) = sum(nfci);
        DC.nf_click_dur(dc) = sum(nf_clicks.end_sec_since_tagon(nfci) - ...
            nf_clicks.start_sec_since_tagon(nfci)) ; 
        else
            DC.nf_click_bouts(dc) = NaN;
            DC.nf_click_dur(dc) = NaN;
        end

        DC.nf_click_bouts(isnan(DC.click_dur)) = NaN;
        DC.nf_click_dur(isnan(DC.click_dur)) = NaN;

        % step length
        % total distance / duration
        if exist("Ma", "var") % this uses the combo_track thus the ptrack
            this_ctrack = combo_track(round(c_start * depth.sampling_rate): ...
                round(c_end * depth.sampling_rate), :);
            diffz   = diff(this_ctrack, 1);
            distz   = sqrt(sum(diffz .* diffz, 2));
            se_track = this_ctrack(distz < 20, :);
            % se_track = this_ctrack([1,end],:);
            se_diff = se_track(end, :) - se_track(1, :);
            se_dist   = sqrt(sum(se_diff .* se_diff, 2));
            DC.path(dc) = sum(distz(distz < 20)) / (c_end - c_start) / 1000 * 3600; % km/h
            DC.step(dc) = se_dist(1) / (c_end - c_start) / 1000 * 3600; % km/h
            % turning angle
            % SKIP FOR NOW - ISSUES WITH IRREG SAMP
        else
            DC.path(dc) = NaN;
            DC.step(dc) = NaN;
        end
        
        % tortuosity (from ptrack)
        DC.tort(dc) = median(tort(tort_sec >= c_start & tort_sec <= c_end) , 'omitnan');
    
    end % end loop over Dive cycles
       
         % Write output to file
    %######################################################################
    writetable(DC, [smrt_path, 'preproc-fine-coarse/', tag, '-coarsescale.csv']);
    clearvars -except whales w smrt_path FMpred HHMMForageClass tag_types
end % loop over whales
