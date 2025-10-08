% this is an example of how RLs were computed. 
% tags were processed as data were collected and audits completed
% set up paths and files
addpath(genpath('/Users/sld33/Dropbox/dtagtools/dtag2'))
addpath(genpath('/Users/sld33/Dropbox/dtag/dtag4'), '-begin')

savepath()
% acoustic audit info
tags = { % 'Zica-20180113-173188'; % lander
    % 'Zica-20180330-173188'; % lander
    % 'Zica-20180331-173187'; % lander
    % 'Zica-20190111-173186'; % lander
    % 'Zica-20190113-151361'; % SMRT but not in qping log, no sounds?
    % 'Zica-20191012-144029'; % in qPing_log
    % 'Zica-20191012-145101'; % in qPing_log
    % 'Zica-20191111-94810'; % in qPing_log
    % 'Zica-20191117-195993'; % in qPing_log
    % % 'Zica-20200106-195994'; % excluded b/c has only one deep dive, no complete dive cycles
    % 'Zica-20211113-195993'; % has individual file
    % 'Zica-20211112-94819'; % has individual file
    % 'Zica-20220112-195994'; % has individual file
    % 'Zica-20230519-232950'; % has individual file
    % 'Zica-20230518-233391'; % has individual file 
    % (note, commented ones have been previously run)
   % 'Bp-20211107-180745'; % 
   % NOT DONE YET! 'Bp-20211107-180747'; % NO NC AVAIL YET -- NOT DONE YET -- must do sep w/o z
   % 'Zica-20211001-220819';% 
    %'Zica-20220930-195992'; % no mfas
   % 'Zica-20221002-237889'; 
   %  'Zica-20230519-232950'; 
    %'Zica-20240227-233396' ; 
    % 'Zica-20240227-240128' ;
    'Zica-20230723-233394';
    'Zica-20230723-233395'
    }; 

% audit_path = '/Volumes/TOSHIBA EXT/FreakinBeakin/';
audit_path = '/Volumes/Boomin/';
% final tag nc-file won't fit on FAT32 formatted external drive :( so on
% mac
nc_path = '/Users/sld33/Dropbox/FBdata/';

pingfiles = {
    %'/Users/sld33/Dropbox/FBdata/RLs/qPing_log_corr_times_master.csv';
    %'Bp-20211107-180745_Individual_anthro_pings.csv';
    %'Bp-20211107-180747_Individual_anthro_pings.csv';
    %'Zica-20211001-220819_Individual_anthro_pings.csv';
    %'Zica-20211112-94819_Individual_MFA_Pings.csv';
    %'Zica-20211113-195993_Individual_MFA_Pings.csv';
    %'Zica-20220112-195994_Individual_MFA_Pings.csv';
    %'Zica-20220930-195992_Individual_anthro_pings.csv';
    %'Zica-20221002-237889_Individual_anthro_pings.csv';
    %'Zica-20230518-233391_Individual_MFA_Pings.csv';
    %'Zica-20230518-233391_Individual_anthro_pings.csv';
    %%'Zica-20230519-232950_Individual_MFA_Pings.csv';
    %'Zica-20230519-232950_Individual_anthro_pings.csv';
    %'Zica-20240227-233396_Individual_anthro_pings.csv';
    %Zica-20240227-240128_Individual_anthro_pings.csv';
    'Zica-20230723-233394_Individual_anthro_pings.csv';
    'Zica-20230723-233395_Individual_anthro_pings.csv'};

for t = 2 %1:length(tags)
    clear this_depth tinfo audit_data ping_times ping_durs
    % get tag record start time
    % try
    %     tinfo = load_nc([audit_path, tags{t}, '/', tags{t}, '-cal'], ...
    %     {'info', 'depth'});
    % catch ME
    %     tinfo = load_nc([nc_path, tags{t}, '-cal'], {'info', 'depth'});
    % end
    % if isempty(tinfo)
        % seems like if the file isn't found it's not an "error"
        tinfo = load_nc([nc_path, tags{t}, '-cal'], {'info', 'depth'});
    % end
    if isempty(tinfo)
        tinfo = load_nc([nc_path, tags{t}, '-working2cal'], {'info', 'P'});
        tinfo.depth = tinfo.P;
    end

    if ~isempty(tinfo)
        this_depth = tinfo.depth;
        tinfo = tinfo.info;
        % tag start time UTC as matlab datenum
        stimenum = datenum(tinfo.dephist_device_datetime_start, tinfo.dephist_device_regset);
    else
        disp(['No nc file found for: ', tags{t}]);
        continue % if no nc file
    end
    
    % read in audit data with ping start times in UTC strings
    
    if any(contains(pingfiles, tags{t}))
        single_fname = pingfiles{contains(pingfiles, tags{t})};
        audit_data = readtable(single_fname);
        mfas_echo_ix = contains(lower(audit_data.Label), 'mfa') | contains(lower(audit_data.Label), 'echo');
        audit_data = audit_data(mfas_echo_ix,:);
        % SORT audit data chronologically
        audit_data = sortrows(audit_data, "UTC");
        % convert vector of ping times to seconds since tag start
        ping_times = (datenum(audit_data.UTC, 'yyyy-mm-ddTHH:MM:SSZ') - stimenum)*24*60*60;
        ping_durs = audit_data.Duration;
        signal_types = audit_data.Label;
        notes = audit_data.Note;
    else
        % read in many-whale combo ping log data
        audit_data = readtable(pingfiles{1});
        audit_data = sortrows(audit_data, {'TagID','sec_since_tagon'},{'ascend','ascend'});
        % allow for filtering of ping log data
        arf = rowfilter(audit_data);
        % filter to just the current whale
        audit_data = audit_data(arf.TagID == tags{t}, :);
        if size(audit_data,1) == 0
            % if this whale has no pings then go on to next one
            continue
        end
        ping_times = audit_data.sec_since_tagon;
        ping_durs = audit_data.Seconds;
        audit_data.UTC = audit_data.Start_Ping;
        signal_types = audit_data.Type;
        notes = audit_data.Note;
    end
    % make min ping dur 0.4 seconds
    ping_durs(ping_durs < 0.4) = 0.4;
    
    % compute RLs and save results as mat and csv files in tag main dir 
    % note "-" is needed in "tag name" FOR LATER # TAGS since it was inserted in dtg file
    % names.
    if matches(tags{t}, {'Zica-20211113-195993'; % has individual file
                'Zica-20211112-94819'; % has individual file
                'Zica-20220112-195994'})
        this_dtg_name = [tags{t}, '-'];
    elseif matches(tags{t}, 'Zica-20221002-237889')
        this_dtg_name = 'Zica-20221003-237889'; % typo in dtg file names
    else
        this_dtg_name = tags{t};
    end

    % not using old multi-tag smrt_rls function b/c of matlab
    % weirdness with char vs string arrays for vectors of paths and depids
    smrt_rls_onetag(this_dtg_name, ... % depid for dtg/wav files
        [audit_path, tags{t}, '/SMFiles/'], ... % recdir
        ping_times, audit_data.UTC, ... % start times in sec since tagon and UTC
        ping_durs, ... % durations
        this_depth, ... % depth data
        6, [],[], 1, '3obank', 'both', ... % minSNR, nst, ndur, transElim, filter, output format
        [nc_path, 'RLs/', tags{t}], ... % output filename
        signal_types,[], [], 0, [], [],...
        notes) ; % signal type, wav clips, f0, plot, sens, progress, notes
end