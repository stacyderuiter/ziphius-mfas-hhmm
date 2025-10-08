function smrt_rls_onetag(depid, recdir, st, st_UTC, dur, depth, minSNR, ...
    nst, ndur, trans_elim, filter_bank, output, out_fname, ...
    signal_type, wav_clips, f0, plot, sens, progress, notes)
%Received level determinations for controlled exposure experiments 
%for freaking beakin data from SMRT tags 
%
%    INPUTS are:
%1.     depid           tag ID string; should match the prefix of the wav files.
%                       Either one string (if all pings are
%                       from the same tag) or a vector the same length as st.
%2.     recdir          path to the directory where the wav files are stored. 
%                       Can be one string (if all pings are from same tag) 
%                       or a vector the same lengths as st
%3.     st              a vector of start times for the
%                       exposure sounds, in seconds since start of tag record.
%                       ASSUMED TO BE IN CHRONOLOGICAL ORDER WITHIN TAG RECORD.
%                       (please sort before calling this function if
%                       needed)
%                       For example, if you exposed the tagged
%                       whale to MFA sonar sounds every 25 seconds starting
%                       3500 seconds after tagon, st would be
%                       something like st = [3500; 3525; 3550; 3575; ...].
%3.1    st_UTC          passed along for inclusion in output table
%4.    dur             optional vector of (maximum) durations of all the CEE sounds in st (in seconds). 
%                       if not provided, the max duration for any sound will be 20 seconds. 
%5.    minSNR          Default: 6dB. SNR threshold used to define signal
%                      duration in a band. The duration of the signal (for 
%                      SEL calculations) is the period from the
%                      first to the last time that the in-band RMS SNR exceeds minSNR, within dur.
%6.    nst            start time of noise clips. default: 1 sec before st
%7.    ndur            duration of noise clips. default: 600 msec (0.6 sec)
%8.    trans_elim    Logical. Default is 1.
%                       0 to measure RL WITHOUT any transient elimination (will
%                       give inaccurate results if any whale clicks overlap the CEE sound).
%                       1 to apply transient elimination algorithm - allows
%                       measurement of RL even if whale clicks overlap the CEE
%                       sound, but NOT useful if the CEE sound includes
%                       transients/clicks.  This algorithm will NOT compensate
%                       for CEE sounds that are covered by longer noises 
%                       such as surface splashing.
%9.     filter_bank     DEFAULT: '3obank', use a filter bank of ANSI standard 
%                       1/3 oct filters
%                       spanning the range 1000 Hz to 40 kHz.
%                       '3o' to use a 1/3 oct filter (if f0 is not specified it will be centered at 3.4 kHz,
%                       spanning 3029 - 3815 Hz.  It is an FIR filter, see
%                       section 2 for details.)
%                       'bb' to make a broadband measurement - bandpass filter
%                       from 1000 Hz to 40 kHz (or afs/2) will be applied.  
%                       See Section 2 for details of filter design.
%10.     output          Default: 'csv' to save a flat text file with results
%                       'mat' to save results as a mat-file
%                       'both' for both .mat and .csv
%11.     out_fname       file name (including path to directory, but not file suffix) in which to save results.
%                       defaults to current working directory and file name "RLs".
%11.5   signal_type      string. signal type for which RL is being
%                       measured. defaults to 'mfas'
%12.     wav_clips       Default: 'none', will NOT save or use wav clip (pull from original data files)
%                       'save' to create a short wav file clip of each ping
%                       and noise window. User will be prompted to choose a
%                       directory in which to save the files.
%13.     f0              optional; center freq **IN HERTZ** for single 1/3 oct filter.
%                       ignored if filter_bank ~= '3o'.
%                       if filter_bank is '3o' and f0 is not provided it defaults to 3.4 kHz.
%                         f0 is the center frequency around which a two element vector will be created that specifies
%                         the lower and upper bounds of a 1/3 octave filter encompassing the signal of interest.  
%                         **IN HERTZ**
%                         Previous values: 
%                         SOCAL RHIBS-only B Jul 9 and 10 2013 -- 3.4 kHz center frequency, [3029 3815]
%                         SOCAL RHIBS-only B Jul 12 2013 -- 2.8 kHz center frequency
%                         SOCAL RHIBS-only B Jul 12 2013 -- 2.7 KHz center frequency
%                         SOCAL scaled source -- 3.7 kHz center frequency
%14.     plot           optional (defaults to 1). If it is 1 a spectrogram
%                       of each ping and noise clip will be saved for later inspection. User will
%                       be prompted to choose the directory in which to save these files.
%15.     sens           tag sensitivity (clip level dB re uPa 0-pk). 
%                       Defaults to 172 (from tag spec for SMRT tag with gain on - 
%                       - if you have calibrated your device(s) use those
%                       values. If sens is a single value, it will be used
%                       for all pings, but a vector the same length as st
%                       can also be provided.
%16.     progress       Logical - default 1. Print info about progress to console?
%        notes          vector of "Notes" from audit file (passed thru to
%                       output)
%
%
% NOTES:
% the rms level reported is the max rms level observed in any 200 ms 
% window during the duration of the signal.
% also reported are peak level, SEL and SNR.  For SNR, rms noise level is
% calculated as for the rms signal level, but using a clip beginning two
% seconds before the start of the CEE sound and ending 1 second before the
% start of the CEE sound.  SEL is calculated over a window where SNR > minSNR 
%
% if notrans = 1, this code implements a version of Walter Zimmer's 
%(NATO NURC La Spezia) "RMS_Filter"
% script, to eliminate transients like clicks.
%
% Stacy DeRuiter, WHOI, May 2011 
% modified @ U of St Andrews, June 2012,
% modified for user-defined input (for naval transmissions), July 2013 A Stimpert
% modified for SMRT tags and csv output, Aug 2020 SDR

%##############################################################################################
%  1.  Preliminaries: enter values for constants, check tag sampling rate,
%       etc.
%##############################################################################################
if nargin < 3 || isempty(st)
    error('st (start times of signals to measure), recdir, and depid are required inputs \n');
end

if nargin < 19 || isempty(progress)
    progress = 1; % change to 0 to disable printing messages to command window during calculations
end

if nargin < 18 || isempty(sens)
    sens  = 172; %from tag spec. for SMRT tags
end

if nargin < 17 || isempty(plot)
    plot = 1; % default is to save plots of each ping
end

if nargin < 16 || isempty(f0)
    f0 = 3400;
end

if nargin < 15 || isempty(wav_clips)
    wav_clips = 'none';
end

if nargin < 14 || isempty(signal_type)
    signal_type = 'mfas';
end

if nargin < 12 || isempty(output)
    output = 'csv';
end

if nargin < 11 || isempty(filter_bank)
    filter_bank = '3obank';
end

if nargin < 13 || isempty(out_fname)
    out_fname = [pwd '\RLs_' filter_bank ];
end


if nargin < 10 || isempty(trans_elim)
    trans_elim = 1;
end

if nargin < 9 || isempty(ndur)
    ndur = 0.6;
end

if nargin < 8 || isempty(nst)
    nst = st - 2;
end

if nargin < 7 || isempty(minSNR)
    minSNR = 6;
end

if nargin < 5 || isempty(dur)
    dur = 20;
end


% for inputs that have a single value, replicate to length(st) if needed
st = check_length(st, st); % will convert from row to col vec if needed
% st_UTC = check_length(st_UTC, st);
dur = check_length(dur, st);
ndur = check_length(ndur, st);
nst = check_length(nst, st);
sens = check_length(sens, st);
signal_type = check_length(signal_type, st);
notes = check_length(notes, st);

% if some marked pings are after tag is off whale, we won't want to measure
% those
post_tagoff = st*depth.sampling_rate > size(depth.data, 1); 
st(post_tagoff) = [];
st_UTC(post_tagoff) = [];
dur(post_tagoff) = [];
ndur(post_tagoff) = [];
nst(post_tagoff) = [];
sens(post_tagoff) = [];
signal_type(post_tagoff) = [];
notes(post_tagoff) = [];

% FOR MFAS make duration max of: 20 and reported duration
% but also min of: 20 and time until next ping
% except if depth is less than 10m, then keep reported dur
dur0 = dur;
max_dur = min([repmat(20, size(dur)), [st(2:end) - st(1:(end-1)); 20]], [], 2);
surf_ix = find(depth.data(round(st * depth.sampling_rate)) < 10);
surf_ipis = [st(2:end) - st(1:(end-1)); 20];
surf_ipis = surf_ipis(surf_ix);
max_dur(surf_ix) = min([dur(surf_ix) , surf_ipis], [], 2);
dur(surf_ix) = max_dur(surf_ix);
dur = max([dur, max_dur], [], 2);
% for echosounder use 
dur(signal_type == "ECHO") = 3 * dur0(signal_type == "ECHO");
% strcat('Durations used:')
% dur

% if using a 1/3 oct filter with user input center freq
% check with user that correct band is being selected
if strcmp(filter_bank, '3o')
    fb = [f0/(2^(1/6)) f0*(2^(1/6))];
    fprintf('Analysis filter band from %1.3f to %1.3f\n **HERTZ**', fb);
    % fprintf('  Accept filter band? Type y or n... ') ;
    [s] = input('  Accept filter band? Type y or n... ','s') ;
    s = char(s) ;
    % fprintf('\n') ;
    if lower(s(1))~='y'
        fprintf('  Rejecting filter band (no RL calculations completed)\n') ;
        return
    end
end

% if saving plots, ask user where to store them
if plot
    fig_dir = uigetdir(pwd, 'Choose a directory to store spectrograms');
end

% check FIRST audio sampling rate
% if iscell(recdir)
%     rd = recdir{1};
% else
%     rd = recdir(1);
% end
% 
% if iscell(depid)
%     did = depid{1};
% else 
%     did = depid(1);
% end

[x,afs] = d3wavread([st(1), st(1) + 1], recdir, depid);
clear x;

% calculations/definitions for constants
peakclip_pa = 10.^(sens/20); %tag peak clip level, converted to Pascals
nf0 = 41; %number of 10 msec periods to average over during application of smoothing filter (to use as baseline for detection of transients)-leads to 200msec window)
thr = 6; %threshold (in dB) above baseline for detection of transients

%##############################################################################################
%  2.  Set up filters 
%##############################################################################################

[B, A, cf, FL] = get_filter(filter_bank, afs, f0);
% note: will need to re-fetch filter coefficients if we move to a new tag
% with different afs.


%##############################################################################################
%  3.  Preallocate space for results
%##############################################################################################
SEL = zeros(length(st),size(B,1)); %SEL in dB re 1 uPa^2*sec
SPL_pk = zeros(length(st),size(B,1)); %peak SPL in dB re 1uPa
SPL_rms = zeros(length(st),size(B,1)); %rms SPL in dB re 1uPa, max obs in any 200 ms window
SNR = zeros(length(st),size(B,1)); % max rms signal to noise ratio in any 10 ms window
hiSNRdur = zeros(length(st),size(B,1)); %duration in sec when  rms SNR in 10 ms windows > minSNR and within 20dB of max level
Ts = zeros(length(st),size(B,1)); %start time of hiSNRdur
Te = zeros(length(st),size(B,1)); %end time of hiSNRdur
noise_rms = zeros(length(st),size(B,1)); %rms noise level...
noise_pk = zeros(length(st),size(B,1)); %peak noise level...
%**********************************************************************************************


%##########################################################################
%  4.  A Read in signal and noise clips
%      B Plot and save spectrograms (if required)
%      C Apply filters
%      D Do transient elimination (if required)
%      E Measure levels
%##########################################################################

if strcmp(wav_clips, 'save')
    workdir = pwd;
    clipdir = uigetdir(workdir, 'Choose the directory where wav clips should be stored.');
    disp('Sorry - saving of wav clips not implemented yet.');
end

%--------------------------------------------------------------------------
% A. Read in signal and noise clips
%--------------------------------------------------------------------------

for pingnum = 1:length(st) %loop over signals to measure
    clear x x_filt x_filt2 xn_filt xn_filt2 tend t1 ww vv3 w v3 v2 vv2
    if progress
        disp(['Measuring RLs (clip ' num2str(pingnum) ' of ' num2str(length(st)) ')']);
    end
    
    afs_init = afs;
    
    %read in signal clip
%     if iscell(recdir)
%         rd = recdir{pingnum};
%     else
%         rd = recdir(pingnum);
%     end
%     
%     if iscell(depid)
%         did = depid{pingnum};
%     else 
%         did = depid(pingnum);
%     end
    
    [x, afs]= d3wavread([st(pingnum), st(pingnum) + dur(pingnum)], recdir, depid);

    
    % read in noise clip
    if pingnum > 1 &&  st(pingnum-1) + dur(pingnum-1) >= nst(pingnum)
        %if previous (ping + buffer) bleeds into noise clip
        % then use previous one
        nst(pingnum) = nstart;
    else 
        % read in this noise clip if not overlapping
        clear xn
        [xn, afsn]= d3wavread([nst(pingnum), nst(pingnum) + ndur(pingnum)], recdir, depid);
        nstart = nst(pingnum);
    end
    %--------------------------------------------------------------------------
    % B.Plot and save spectrograms (if required)
    %--------------------------------------------------------------------------
    if plot
        if dur > 0.005 % if to short to plot, skip
            if dur > 0.05
                BL = 4*512 ;
            else
                BL = 256; % specgram (fft) block size
            end
            window = hanning(BL);
            noverlap = length(window)/2;
            
            figure(1), clf
            set(gcf, 'color', [1 1 1]);
            
            AXm = axes('position',[0.25,0.11,0.60,0.8]) ;
            AXn = axes('position',[0.11,0.11,0.10,0.8]) ;
            
            [s,f,t] = spectrogram(x, window, noverlap, BL, afs);
            [sn,fn,tn] = spectrogram(xn, window, noverlap, BL, afsn);
            
            axes(AXm);
            imagesc(t , f/1000, 20*log10(abs(s)));
            set(AXm,'ydir','normal');
            xlabel('Time (sec)')
            ylabel('')
            title(['Tag: ', depid, '; Time: ', num2str(st(pingnum)), ' sec.']);
            s = []; f = []; t = [];
            
            axes(AXn);
            imagesc(tn - 2, fn/1000, 20*log10(abs(sn)));
            set(AXn,'ydir','normal');
            xlabel('')
            ylabel('Frequency (kHz)')
            title([num2str(nstart) ' sec.']);
            sn = []; fn = []; tn = [];
            
            
            
            figfile = [fig_dir, '/', depid, '-cst', num2str(st(pingnum))];
            pngfile = [figfile '.png'];
            
            % savefig(gcf, figfile);
            % supposedly faster?
            % screencapture(gcf, [], pngfile);
            print(gcf, pngfile, '-r300', '-dpng');
            
        end
    end
    
    %--------------------------------------------------------------------------
    % C. (Re-compute if needed and) Apply filters
    %--------------------------------------------------------------------------
    
    % check if afs is the same as before, and if not, re-do filter design
    if (afs_init ~= afs)
        [B, A, cf, FL] = get_filter(filter_bank, afs, f0);
    end
    
    for k2 = 1:size(B,1) %and loop over filters to be applied (for 3obank filter_bank option)
        clear x_filt xn_filt x_filt2 xn_filt2
        if ~isempty(A) %butterworth bandpass filter
            x_filt = filter(B(k2,:), A(k2,:), x(:,1)); %apply the bandpass butterworth filter to the signal
            xn_filt = filter(B(k2,:), A(k2,:), xn(:,1));
        elseif ~isempty(FL) % FIR filter
            x_filt = fftfilt(B,[x(:,1); zeros(FL,1)]); %apply fft filter
            xn_filt = fftfilt(B,[xn(:,1); zeros(FL,1)]); %apply fft filter
        end
        
        x_filt2 = x_filt.^2; %convert from amplitude to intensity
        xn_filt2 = xn_filt.^2; %convert from amplitude to intensity
        
        % apply smoothing filter...
        % to signal
        nf = nf0;
        nf1 = 1 + (nf-1) / 2;
        vv = zeros(length(ceil(afs*(0:0.01:length(x_filt)/afs - 0.2))));
        for p = 1:length(vv)
            vv(p) = (mean(x_filt2((p-1)*(afs/200)+(1:(2*afs/100)))));
        end
        vv2 = filter(ones(nf, 1) / nf, 1, [vv; zeros(nf, size(vv,2))]);
        vv2(1 : nf1, :) = [];
        vv2 = vv2(1 : size(vv, 1), :);
        
        % apply same smoothing filter...
        % to noise clip
        v = zeros(length(ceil(afs*(0:0.01:length(xn_filt)/afs - 0.2))));
        for p = 1:length(v);
            v(p) = (mean(xn_filt2((p-1)*(afs/200)+(1:(2*afs/100))))); %average over 10 msec durations, with start of interval stepping forward by 5 msec
        end
        % smooth v with a filter spanning 41 sampled periods = 200 msec
        v2 = filter(ones(nf, 1) / nf, 1, [v; zeros(nf, size(v, 2))] );
        v2(1 : nf1, :) = [];
        v2 = v2(1 : size(v,1), :);
        
        %--------------------------------------------------------------------------
        % D. Transient elimination
        %--------------------------------------------------------------------------
        if trans_elim
            % SIGNAL CLIP
            itr = vv > vv2 * 10^(thr/10); % index of where entries of vv > vv2*10^(thr/10), thr = 4
            vv3 = vv .* (1 - itr); % vv with zeros where there are transients
            %widen holes to where sample values are same as noise level
            for ii = 1:size(vv3, 2)
                for jj = 2:size(vv3, 1) - 1
                    if (vv3(jj - 1, ii) == 0) && (vv3(jj, ii) > vv3(jj + 1, ii))
                        vv3(jj, ii) = vv3(jj - 1, ii);
                    end
                end
                for jj = size(vv3, 1) - 1:-1:2
                    if (vv3(jj + 1, ii) == 0) && (vv3(jj, ii) > vv3(jj - 1, ii))
                        vv3(jj, ii) = vv3(jj + 1, ii);
                    end
                end
            end
            %fill up holes with sample @ background/previous level
            for ii = 1:size(vv3, 2)
                for jj = 2:size(vv3, 1)
                    if vv3(jj, ii) == 0,
                        vv3(jj, ii) = vv3(jj - 1, ii);
                    end
                end
            end
            
            % NOISE CLIP
            itrn = v > v2 * 10^(thr/10); %index of where transients are = where entries of v > v2*10^(thr/10)
            v3 = v .* (1 - itrn); % fill vv with zeros where there are transients
            %widen holes to where sample values are same as noise level
            for ii = 1 : size(v3, 2)
                for jj = 2 : size(v3, 1) - 1
                    if (v3(jj - 1, ii) == 0) && (v3(jj, ii) > v3(jj + 1, ii))
                        v3(jj, ii) = v3(jj - 1, ii);
                    end
                end
                for jj = size(v3, 1) -1:-1:2
                    if (v3(jj + 1, ii) == 0) && (v3(jj, ii) > v3(jj - 1, ii))
                        v3(jj, ii) = v3(jj + 1, ii);
                    end
                end
            end
            %fill up holes with sample @ background/previous level
            for ii = 1:size(v3, 2)
                for jj = 2:size(v3, 1)
                    if v3(jj, ii) == 0
                        v3(jj, ii) = v3(jj - 1, ii);
                    end
                end
            end
        else
            v3 = v;
            vv3 = vv;
        end % end of transient elimination
        
        % signal clip
        ww = filter(ones(nf, 1) / nf, 1, [vv3; zeros(nf, size(vv3, 2))]); % vv3 filtered w/200msec smoothing filter
        ww(1:nf1, :) = [];
        ww = ww(1:size(vv3, 1), :);
        % noise clip
        wn = filter(ones(nf, 1) / nf, 1, [v3; zeros(nf, size(v3, 2))]); %v3 filtered w/100msec smoothing filter
        wn(1:nf1, :) = [];
        wn = wn(1:size(v3, 1), :);
        
        %--------------------------------------------------------------------------
        % E. Compute levels
        %--------------------------------------------------------------------------
        noise_rms(pingnum,k2) = max( 20*log10(sqrt(mean(wn))) + sens(pingnum));
        noise_pk(pingnum,k2) = max(10*log10(max(abs(v3))) + sens(pingnum));
        
        % cee clip level determination -- using SMOOTHED, transient-removed data (ww)
        %Calculate the SEL for the entire time for which SNR is > minSNR
        %first determine SNR in 10 milli second windows ...
        %find the time window in which SNR > minSNR
        %Calculate the SEL for the entire time for which SNR is > minSNR
        %first determine SNR in 10 milli second windows
        q = 1:length(vv3);
        %find the time window in which SNR > minSNR
        snrs = zeros(length(q),1);
        for tt = 1:length(q)
            snrs(tt) =  max(10*log10(max(max(ww(tt:min( max(q) , tt+3 ),:)))) + sens(pingnum) - noise_rms(pingnum,k2));
        end
        t1 = q(snrs > minSNR); %indices of the 10-msec (overlapping) windows with SNR > minSNR
        if ~isempty(t1)
            hsw = snrs(t1); %vector of high SNRs from time of first SNR>minSNR to last SNR>minSNR
            THdfm = 20; %threshold (decline from max observed) for picking the end of the hi SNR window.
            % shorten hiSNR window so it spans time from first time minSNR is exceeded
            % until the LAST time the level is within THdfm dB below the highest measured level
            endwind = find(max(hsw) - hsw < THdfm,  1, 'last');
            if ~isempty(endwind) %if SNR never falls THdfm dB below max within the hiSNR window, then keep the original hiSNR window.
                t1 = t1(1:endwind);
            end
        end
        
        % % get rid of any hisnrdur parts that follow a gap of more than 40 samples = 0.005*40 sec = 0.2 sec
        % if ~isempty(t1)
        %     gap = find(diff([(t1(1)-1);t1(:)]) > 40, 1, 'first');
        %     if ~isempty(gap)
        %         t1(gap:end) = [];
        %     end
        % end
        
        % not sure when this would help??
        % if ~isempty(t1) && t1(1)*0.005 > 0.3 % if "start time " of signal is more than 0.3 sec in, discard.
        %     t1 = [];
        % end
        
        % note start/end times of hi SNR window
        if ~isempty(t1)
            tst = t1(1); %start index of hiSNR window
            tend = t1(end); %end index of hi SNR window
            Ts(pingnum,k2) = st(pingnum) + tst*0.005; %time steps in vv3=q=t1 are 5 ms apart
            Te(pingnum,k2) = st(pingnum) + tend*0.005;
            
        else
            tst = 0; tend = 0;
            Ts(pingnum,k2) = 0;
            Te(pingnum,k2) = 0;
        end
        
        %record SNR (max found in any 10 msec window)
        try
            % this will throw an error if the sound is super short like
            % 0.05 sec (in that case report NaN)
            SNR(pingnum,k2) = max(snrs);
        catch
        end
        
        %next calculate the SEL and peak & rms SPL for the time when rms SNR is greater than minSNR
        if tst == tend %if there is no signal with SNR > minSNR
            SEL(pingnum,k2) = NaN;
            SPL_pk(pingnum,k2) = NaN;
            SPL_rms(pingnum,k2) = NaN;
        else
            %SEL in the hi snr window
            SEL(pingnum, k2) = max( 10*log10(sum(((ww(tst:tend, :) .* peakclip_pa(pingnum)^2)) ./ 200)) );
            %rms spl in the hi snr window
            SPL_rms(pingnum, k2) = (10 * log10(max(max(abs(ww(tst:tend, :))))) + sens(pingnum));
            %peak level
            SPL_pk(pingnum, k2) = max(10 * log10(max(abs(vv3))) + sens(pingnum));
        end
        
        hiSNRdur(pingnum, k2) = (tend - tst)*0.005; %duration of the signal with SNR > minSNR, in seconds
    end
end

%##########################################################################
%  5.  Save results.
%##########################################################################

if size(SPL_rms, 2) == 1
    % if there is only one filter there is one column per measurement
    % label it with the fc used anyway
    levels = table(SEL, SPL_rms, SPL_pk, noise_rms, noise_pk, SNR, hiSNRdur, Ts, Te);
    levels.Properties.VariableNames = strcat(levels.Properties.VariableNames, ...
        '_fc_', num2str(cf));
else
    % if there is a bank of 1/3 oct filters, there is one column per filter
    % per metric
    cf_cell = compose('%d', cf);
    levels = [array2table(SEL, 'VariableNames', strcat('SEL_fc_', cf_cell)), ...
                    array2table(SPL_rms, 'VariableNames', strcat('SPL_rms_fc_', cf_cell)), ...
                    array2table(SPL_pk, 'VariableNames', strcat('SPL_pk_fc_', cf_cell)), ...
                    array2table(noise_rms, 'VariableNames', strcat('noise_rms_fc_', cf_cell)), ...
                    array2table(noise_pk, 'VariableNames', strcat('noise_pk_fc_', cf_cell)), ...
                    array2table(SNR, 'VariableNames', strcat('SNR_fc_', cf_cell)), ...
                    array2table(hiSNRdur, 'VariableNames', strcat('hiSNRdur_fc_', cf_cell)), ...
                    array2table(Ts, 'VariableNames', strcat('Ts_fc_', cf_cell)), ...
                    array2table(Te, 'VariableNames', strcat('Te_fc_', cf_cell)) ];
end

if size(depid, 1) == 1
    % remove the extra - from tag ID if it was added to facilitate wavread
    if depid(end) == '-'
        depid = depid(1:(end-1));
    end
    depid = repmat(depid, size(st,1), 1);
end

meta_table = table(depid, st, st_UTC, nst, dur, signal_type, notes);
RL_table = [meta_table, levels];
% add more metadata about how calcs were done
RL_table.transient_elim = repmat(trans_elim, size(RL_table, 1), 1);
RL_table.minSNR = repmat(minSNR, size(RL_table, 1), 1);
RL_table.filter_type = repmat(filter_bank, size(RL_table, 1), 1);


% save file(s)
if strcmp(output,'mat') || strcmp(output, 'both')
    save([out_fname '.mat'], 'RL_table');
end

if strcmp(output, 'csv') || strcmp(output, 'both')
    writetable(RL_table, [out_fname '.csv']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Helper functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function x = check_length(x, st)
        % check length of input to make sure it's 1, or matches st.
        % if it's length one, replicate to same size as st.
        
        % if input is one row and ncol ~= n entries in st
        if size(x, 1) == 1
            if size(x, 2) ~= length(st)
                x = repmat(x, length(st), 1);
            else
                % change row vec to col vec
                x = x';
            end
        end
        
        if size(x, 1) ~= length(st)
            error('tag, recdir, dur, and sens should either be single values or vectors the same length as st.');
        end
     end

    function [B,A, cf, FL] = get_filter(filter_bank, afs, f0)
        clear B A cf FL% get rid of old values each time this function is called
        if strcmp(filter_bank, '3o')
            %generate FIR 1/3 oct filter
            FL = 512; %fir filter length
            fbb = [f0/(2^(1/6)) f0*(2^(1/6))];
            B = fir1(FL, fbb./(afs/2)); %FIR filter - %centered at f0
            A = [];
            cf = f0;
        elseif strcmp(filter_bank, '3obank')
            %   A. determine center frequencies for 1/3 octave filters 
            % which satisfy ANSI S1.6-R2006 and ANSI S1.11-R2014.
            top_f = min([40000, 0.75*afs/2]);
            [fexact_lcu, fnominal_lcu, fstrings] = fract_oct_freq_band(3, 1000, top_f);
            good_freqs = find(fexact_lcu(:,3) <= afs/2);
            fexact_lcu = fexact_lcu(good_freqs,:);
            fnominal_lcu = fnominal_lcu(good_freqs,:);
            fstrings = fstrings(good_freqs,:);
            %   B. design 1/3 oct filters with the above center freqs (given in Hz), to
            %   satisfy ANSI
            cf = fnominal_lcu(:,2);
            %preallocate space for filter coeffs
            B = zeros(length(cf),7);
            A = zeros(length(cf),7);
            FL = [];
            %generate 1/3 oct filter coefficients
            for g = 1:length(cf)
                [B(g,:),A(g,:)] = oct3dsgn(fexact_lcu(g,2), afs, 3); %6th order butterworth bandpass filter
            end
        elseif strcmp(filter_bank, 'bb')
            top_f = min(40000, afs/2);
            %generate bandpass "broadband" filter coefficients
            [B,A] = butter(3,[1000/(afs/2) top_f/(afs/2)]); %6th order butterworth bandpass filter from 250 Hz to 15 kHz
            cf = [];
            FL = [];
        else error('Unknown value for input filter_bank - please consult help for this script');
        end
    end

    function [fexact_l_c_u, fnominal_l_c_u, fnominal_str_l_c_u] = fract_oct_freq_band(N, min_f, max_f, SI_prefixes, base10, fr)
% Copyright notice for "Nth octave freq bands" functions:
% Copyright (c) 2019, Edward Zechmann
% Copyright (c) 1997, Christophe COUVREUR
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of NIOSH nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

        % % fract_oct_freq_band: Calculates the 1/nth octave band center, lower, and upper bandedge frequency limits
        % %
        % % Syntax;
        % %
        % % [fexact_l_c_u, fnominal_l_c_u, fnominal_str_l_c_u] = fract_oct_freq_band(N, min_f, max_f, SI_prefixes, base10, fr);
        % %
        % % **********************************************************************
        % %
        % % Description
        % %
        % % This program calculates the fractional, 1/N octave band center
        % % frequencies and the lower and upper bandedge limits of every
        % % frequency band from min_f to max_f (Hz).
        % %
        % % The exact and nominal values of the frequencies are output.
        % % The exact values are useful for calculations and making filters.
        % % The nominal values are useful for tables, plots, figures, and graphs.
        % %
        % % The frequency band nominal values are rounded to three or more
        % % significant digits based on the number of bands per octave.
        % % Full octave and third octave bands are rounded to three significant digits.
        % %
        % % The exact and nominal frequency values and character strings are output.
        % %
        % % fr default 1000 Hz
        % % fr Hz is a reference frequency when N is odd, odd-fractional-octave-bands.
        % % fr Hz is an edge frequency when N is even, even-fractional-octave-bands.
        % %
        % % fract_oct_freq_band produces exact and nominal frequency bands which
        % % satisfy ANSI S1.6-R2006 and ANSI S1.11-R2014.
        % %
        % %
        % % fract_oct_freq_band is a modification of centr_freq
        % % centr_freq can be found on Matlab Central File Exchange
        % % The Matlab ID is 17590.
        % %
        % % **********************************************************************
        % %
        % % Input Variables
        % %
        % % N=3;            % one-third-octave-bands which have three bands per octave.
        % %                 % 1 one band per octave.
        % %                 % N is the number of frequency bands per octave.
        % %                 % N can be any integrer > 0.
        % %                 % Default is 3 for one-third-octave-bands.
        % %
        % % min_f=20;       % min_f is the minimum frequency band to calculate (Hz).
        % %                 % min_f > 0. Must be graeater than 0.
        % %                 % default is 20;
        % %
        % % max_f=20000;    % max_f is the maximum frequency band to calculate (Hz).
        % %                 % max_f > 0. Must be graeater than 0.
        % %                 % default is 20000;
        % %
        % % SI_prefixes=0;  % This parameter controls the output for the character
        % %                 % strings of the nominal frequencies.
        % %                 % 1 for using SI prefixes (i.e. k for 1000)
        % %                 % 0 for not using prefixes.
        % %                 % default is SI_prefixes=0;
        % %
        % % base10=1;       % 1 for using base 10 filter frequencies
        % %                 % 0 for using base 2 filter frequencies
        % %                 % otherwise  base 10 filter frequencies
        % %                 % default base 10 filter frequencies
        % %
        % % fr=1000;        % fr is the preferred central frequency for
        % %                 % calculating all of the other center frequencies.
        % %                 % default is 1000;
        % %
        % % **********************************************************************
        % %
        % % Output Variables
        % %
        % % fexact_l_c_u is a 2D array of exact frequency bands in Hz.
        % % fexact_l_c_u( frequecy bands, [lower=1, center=2, upper =3;] );
        % %
        % % fnominal_l_c_u is a 2D array of nominal frequency bands in Hz.
        % % fnominal_l_c_u{ frequecy bands, [lower=1, center=2, upper =3;] };
        % %
        % % fnominal_str_l_c_u is a 2D cell array of nominal frequency band strings in Hz.
        % % fnominal_str_l_c_ufnominal_l_c_u{ frequecy bands, [lower=1, center=2, upper =3;] };
        % %
        % % **********************************************************************
        %
        %
        % Example='1';
        %
        % % Full-octave-band center frequencies from 20 Hz to 20000 Hz
        % close('all');
        % [fexact_l_c_u, fnominal_l_c_u, fnominal_str_l_c_u] = fract_oct_freq_band(1, 20, 20000);
        % figure(1);
        % semilogy(fexact_l_c_u(:, 1), ':ro', 'linewidth', 4, 'markersize', 10);
        % hold on;
        % semilogy(fexact_l_c_u(:, 2), '-gp', 'linewidth', 4, 'markersize', 10);
        % semilogy(fexact_l_c_u(:, 3), '-.cs', 'linewidth', 4, 'markersize', 10);
        % semilogy(fnominal_l_c_u(:, 1), '-.kx', 'linewidth', 3, 'markersize', 14);
        % semilogy(fexact_l_c_u(:, 2), ':m*', 'linewidth', 3, 'markersize', 14);
        % semilogy(fexact_l_c_u(:, 3), '--b+', 'linewidth', 3, 'markersize', 10);
        % legend({'Exact lower band', 'Exact center', 'Exact upper band', ...
        %         'Nominal lower band','Nominalcenter', 'Nominal upper band'}, ...
        %         'location', 'northwest');
        % set(gca, 'fontsize', 24);
        % maximize(gcf);
        %
        %
        % % one-third-octave-band center frequencies from 20 Hz to 20000 Hz
        % [fexact_l_c_u, fnominal_l_c_u, fnominal_str_l_c_u] = fract_oct_freq_band(3, 20, 20000);
        %
        % % one-twelveth-octave-band center frequencies from 100 Hz to 10000 Hz
        % [fexact_l_c_u, fnominal_l_c_u, fnominal_str_l_c_u] = fract_oct_freq_band(12, 100, 10000);
        %
        % % twenty-fourth-octave-band center frequencies from 0.001 Hz to 10000000 Hz
        % [fexact_l_c_u, fnominal_l_c_u, fnominal_str_l_c_u] = fract_oct_freq_band(24, 0.001, 10000000);
        %
        %
        % % **********************************************************************
        % %
        % % References
        % %
        % % 1)  ANSI S1.6-R2006 Preferred Frequencies, Frequency Levels, and
        % %     Band Numbers for Acoustical Measurements, 1984.
        % %
        % % 2)  ANSI S1.11-2014, (2014). American National Standard Electroacoustics
        % %     - Octave-band and Fractional-octave-band Filters - Part 1:
        % %     Specifications (a nationally adopted international standard).
        % %     American National Standards Institute, Acoustical Society of
        % %     America, New York
        % %
        % %
        % % **********************************************************************
        % %
        % % fract_oct_freq_band is a modification of centr_freq
        % % centr_freq can be found on Matlab Central File Exchange
        % % The Matlab ID is 17590.
        % %
        % % Written by   Author  Eleftheria
        % %              E-mail  elegeor@gmail.com
        % %              Company/University: University of Patras
        % %
        % %
        % % **********************************************************************
        % %
        % % List of Dependent Subprograms for
        % % fract_oct_freq_band
        % %
        % % FEX ID# is the File ID on the Matlab Central File Exchange
        % %
        % %
        % % Program Name   Author   FEX ID#
        % % 1) sd_round		Edward L. Zechmann
        % %
        % %
        % % **********************************************************************
        % %
        % %
        % % Program Modified by Edward L. Zechmann
        % %
        % % modified  3 March       2008    Original modification of program
        % %                                 updated comments
        % %
        % % modified 13 August      2008    Updated comments.
        % %
        % % modified 18 August      2008    Added rounding last digit to nearest
        % %                                 multiple of 5.  Added Examples.
        % %                                 Updated comments.
        % %
        % % modified 21 August      2008    Fixed a bug in usign min_f and max_f
        % %                                 which does not include 1000 Hz.
        % %                                 Zipped the depended program sd_round.
        % %                                 Updated comments.
        % %
        % % modified 18 November    2008   	Added additional rounding
        % %
        % % modified  8 December    2008    Updated comments.
        % %
        % % modified 18 December    2008    Updated comments.
        % %
        % % modified  4 January     2008    Changed rounding of the lower and upper
        % %                                 bandedge frequency limits.
        % %
        % % modified  6 October     2009    Updated comments
        % %
        % % modified 22 January     2010    Modified the number of significant
        % %                                 digits for rounding. The number of
        % %                                 significnat digits increases as the
        % %                                 number of bands per octave increases.
        % %                                 This supports high resolution octave
        % %                                 band analysis.
        % %
        % % modified  4 October     2014    Modified the even-octave-bands to have
        % %                                 1000 Hz as an edge frequency.
        % %                                 Changed the number of significant
        % %                                 digits calculation for rounding.
        % %                                 Limited the number of bands per octave
        % %                                 from 1 to 43.
        % %
        % % modified 21 February    2019    Revised the number of significant
        % %                                 digits calculation for rounding.
        % %                                 Limited the number of bands per octave
        % %                                 from 1 to 10E12.
        % %
        % %                                 A warning is issued if the number of
        % %                                 bands is greater than 10E6.
        % %
        % % modified 22 February    2019    Added the reference frequency, fr as an
        % %                                 input variable and optimized the
        % %                                 removal of fractional-octave-bands
        % %                                 outside the range of min_f to max_f.
        % %
        % %
        % % **********************************************************************
        % %
        % % Please Feel Free to Modify This Program
        % %
        % % See Also: centr_freq, sd_round, nth_freq_band
        % %
        if (nargin < 1 || isempty(N)) || ~isnumeric(N)
            N=3;
        end
        % N must be an integrer;
        N=round(N);
        % N is limited to between 1 and 10^12.
        % Computer seems to lock up for N > 10^12
        N(N < 1)=1;
        N(N > 10^12)=10^12;
        if (nargin < 2 || isempty(min_f)) || (logical(min_f < 0) || ~isnumeric(min_f))
            min_f=20;
        end
        if (nargin < 3 || isempty(max_f)) || (logical(max_f < 0) || ~isnumeric(max_f))
            max_f=20000;
        end
        if (nargin < 4 || isempty(SI_prefixes)) ||  ~isnumeric(SI_prefixes)
            SI_prefixes=0;
        end
        if (nargin < 5 || isempty(base10)) || ~isnumeric(base10)
            base10=1;
        end
        if base10 ~= 0
            base10=1;
        end
        if (nargin < 6 || isempty(fr)) || ~isnumeric(fr)
            fr=1000;
        end
        % This program uses different symbols than in the standard.
        %
        % N is used as the band designator.
        % In ANSI S1.11-2014 the letter b is the band designator.
        % Determine if N is odd or even.
        % Noe stands for N-odd-even
        Noe=mod(N-1, 2);
        % Noe equals 0 when N is odd
        % Noe equals 1 when N is even
        % Now the nominal frequencies can be calculated
        % Calculate the base number G
        if base10 == 1
            G=10^(3/10);
        else
            G=2;
        end
        % Estimate the number of bands being output and return a warning if the
        % array is too large.
        Nbands_total=ceil(N.*(log(max_f)-log(min_f))./log(G));
        if Nbands_total > 10^6
            warning('output array estimated to be larger than 10E6 elements')
        end
        if fr > min_f
            % going down to get just below min_f
            % a minus sign indicates going down
            Num_bands_root_to_fmin = -ceil(N.*(log(fr)-log(min_f))./log(G));
        else
            % going up to get just below min_f
            % a plus sign indicates going up
            Num_bands_root_to_fmin = +floor(N.*(log(min_f)-log(fr))./log(G));
        end
        % check if band is just below min_f
        log_fmin_est=log(fr)-(Noe/((1+Noe)*N)).*log(G)+Num_bands_root_to_fmin./N.*log(G);
        while log_fmin_est < log(min_f+eps)
            Num_bands_root_to_fmin=Num_bands_root_to_fmin+1;
            log_fmin_est=log(fr)-(Noe/((1+Noe)*N)).*log(G)+Num_bands_root_to_fmin./N.*log(G);
        end
        while log_fmin_est > log(min_f+eps)
            Num_bands_root_to_fmin=Num_bands_root_to_fmin-1;
            log_fmin_est=log(fr)-(Noe/((1+Noe)*N)).*log(G)+Num_bands_root_to_fmin./N.*log(G);
        end
        if fr > max_f
            % going down to get just above max_f
            % a minus sign indicates going down
            Num_bands_root_to_fmax = -floor(N.*(log(fr)-log(max_f))./log(G));
        else
            % going up to get just above max_f
            % a plus sign indicates going up
            Num_bands_root_to_fmax = +ceil(N.*(log(max_f)-log(fr))./log(G));
        end
        % check if band is just above max_f
        log_fmax_est=log(fr)+(Noe/((1+Noe)*N)).*log(G)+Num_bands_root_to_fmax/N.*log(G);
        while log_fmax_est > log(max_f+eps)
            Num_bands_root_to_fmax=Num_bands_root_to_fmax-1;
            log_fmax_est=log(fr)+(Noe/((1+Noe)*N)).*log(G)+Num_bands_root_to_fmax/N.*log(G);
        end
        while log_fmax_est < log(max_f-eps)
            Num_bands_root_to_fmax=Num_bands_root_to_fmax+1;
            log_fmax_est=log(fr)+(Noe/((1+Noe)*N)).*log(G)+Num_bands_root_to_fmax/N.*log(G);
        end
        log_fc = log_fmin_est:(log(G)/N):log_fmax_est;
        % Remove bands above the extended max_f limit
        log_fc = log_fc( log_fc < log(max_f)+0.5*log(G)/N );
        % Remove bands below the extended min_f limit
        log_fc = log_fc( log_fc > log(min_f)-0.5/N*log(G) );
        % Make center frequency bands unique
        log_fc = unique(log_fc);
        log_fc = log_fc(:);
        % Calculate the number of center frequency bands
        num_bands=length(log_fc);
        fc=exp(log_fc);
        % fc is the array of exact center frequencies;
        fc_exact=fc;
        % ***********************************************************************
        %
        % Calculate the lower and upper bounds for the lower and upper
        % band edge frequencies
        flu = exp([ log_fc(1)-1./(2.*N)*log(G); log_fc+(1./(2.*N)).*log(G) ]);
        % Form the numeric arrays for the exact lower frequency band edge limits.
        fl_exact=flu(1:num_bands);
        % Form the numeric arrays for the exact upper frequency band edge limits.
        fu_exact=flu(1+(1:num_bands));
        % ************************************************************************
        %
        % Calculate the number of significant digits for rounding the exact
        % frequencies into the nominal frequencies.
        %
        % Full-octave-band and one-third-octave-band must have 3 significant digits.
        % the number of significant digits increases as log10(N).
        %
        num_sd_digits=floor(3+log10(N));
        m=10.^(num_sd_digits-1);
        % ************************************************************************
        %
        % Calculations for nominal frequencies
        %
        % Apply appropriate rounding to the center frequencies
        [fc,   fc_str]   = sd_round(fc, num_sd_digits, 1, 5, SI_prefixes);
        [fc2,  fc_str2]  = sd_round(fc, num_sd_digits, 1, 100, SI_prefixes);
        % If the center frequency rounded to 3 significant digits and the last
        % digit rounded to the nearest multiple of 5 is within 1% of the center
        % frequency rounded to 1 significant digit, then round to 1 significant
        % digit.
        %
        ix=find(abs(m*(1-exp( log(fc)-log(fc2) ))) < 1);
        fc(ix)=fc2(ix);
        fc_str(ix)=fc_str2(ix);
        % Nominal frequencies
        fc_nom=fc;
        fc_nom_str=fc_str;
        % Apply the same rounding technique to the lower and upper frequency
        % bandedge limits as was applied to the center frequencies.
        [flu,   flu_str]   = sd_round(flu, num_sd_digits, 1, 5, SI_prefixes);
        [flu2,  flu_str2]  = sd_round(flu, num_sd_digits, 1, 100, SI_prefixes);
        % If the center frequency rounded to 3 significant digits and the last
        % digit rounded to the nearest multiple of 5 is within 1% of the center
        % frequency rounded to 1 significant digit, then round to 1 significant
        % digit.
        %
        ix=find(abs(m*(1-exp( log(flu)-log(flu2) ))) < 1);
        flu(ix)=flu2(ix);
        flu_str(ix)=flu_str2(ix);
        % Form the numeric and string arrays for the nominal lower frequency
        % band edge limits.
        fl_nom=flu(1:num_bands);
        fl_nom_str=flu_str(1:num_bands);
        % Form the numeric and string arrays for the nominal upper frequency
        % band edge limits.
        fu_nom=flu(1+(1:num_bands));
        fu_nom_str=flu_str(1+(1:num_bands));
        % Concatenate the outputs
        fexact_l_c_u=[fl_exact, fc_exact, fu_exact];
        fnominal_l_c_u=[fl_nom, fc_nom, fu_nom];
        fnominal_str_l_c_u=cell(num_bands, 3);
        for e1=1:num_bands
            fnominal_str_l_c_u{e1, 1}=fl_nom_str{e1};
            fnominal_str_l_c_u{e1, 2}=fc_nom_str{e1};
            fnominal_str_l_c_u{e1, 3}=fu_nom_str{e1};
        end
    end

    function [A2, A_str, real_digitsL, real_digitsR, imag_digitsL, imag_digitsR]=sd_round(A, N, flag2, mult, SI_prefixes)
        % % sd_round: Rounds an array to a specified number of Significant Digits, significant figures, digits of precision
        % %
        % % *********************************************************************
        % %
        % % Syntax;
        % %
        % % [A2, A_str, real_digitsL, real_digitsR, imag_digitsL,...
        % % imag_digitsR]=sd_round(A, N, flag2);
        % %
        % % *********************************************************************
        % %
        % % Description
        % %
        % % sd_round stands for "Significant Digits Round".
        % %
        % % This program rounds a 2-d matrix of numbers to a specified number
        % % of significant digits.
        % %
        % % This program support five different styles of rounding the last digit:
        % % to the nearest integer, up, down, toward zero, and away from zero.
        % %
        % % This program supports real and complex numbers.
        % %
        % % The program outputs the rounded array, a cell string of the
        % % rounded matrix the number of digits, to the left and right of the
        % % decimal place.
        % %
        % % This program is useful for presenting scientific data that
        % % requires rounding to a specific number of significant digits
        % % for publication.
        % %
        % % Significant digits are counted starting from the first non-zero
        % % digit from the left.
        % %
        % %
        % % *********************************************************************
        % %
        % % Input variables
        % %
        % % A is the input matrix of number to be rounded.
        % %      default is empty array A=[];.
        % %
        % % N is the number of significant digits.
        % %      default is N=3;
        % %
        % % flag2 specifiies the style of rounding.
        % %      This program supports four different styles of rounding.
        % %      flag2 == 1 rounds to the nearest integer
        % %      flag2 == 2 rounds up
        % %      flag2 == 3 rounds down
        % %      flag2 == 4 rounds toward zero
        % %      flag2 == 5 rounds away from zero
        % %      otherwise round to the nearest integer
        % %      default is round to the nearest integer
        % %
        % % mult is a whole number.  The program rounds the last digit to mult.
        % %      It is preferred that mult be between 1 and 9; however, all whole
        % %      numberS >= 1 are valid input.  The program rounds mult to the
        % %      nearest integer and makes sure the value is at least 1.
        % %      default is mult=1;
        % %
        % % SI_prefixes=0;  % 1 for using SI prefixes (i.e. K for 1000)
        % %                 % 0 for not using prefixes.
        % %                 % default is SI_prefixes=0;
        % %
        % %
        % % *********************************************************************
        % %
        % % Output variables
        % %
        % % A2 is the rounded array.
        % %
        % % A_str           % The rounded array is converted to a cell
        % %                 % string format with the specified rounding and showing
        % %                 %  the trainling zeros.
        % %                 % This is convenient for publishing tables in a tab
        % %                 % delimited string format
        % %
        % % real_digitsL    % The number of real digits to the left of the decimal
        % %                 % point
        % %
        % % real_digitsR    % The number of real digits to the right of the decimal
        % %                 % point
        % %
        % % imag_digitsL    % The number of imaginary digits to the left of the
        % %                 % decimal point
        % %
        % % imag_digitsR 	% The number of imaginary digits to the right of the
        % %                 % decimal point
        % %
        % % *********************************************************************
        % %
        %
        % Example1='';
        % D1=pi;        % Double or Complex two dimensional array of numbers
        % N=3;          % Number of significant digits.  3 is the default
        %
        % [P1, P1_str]=sd_round(pi, N, 4);
        %
        % % P1 should be 3.14 which has 3 significant digits.
        %
        % Example2='';
        % D1=pi/1000000;    % Double much smaller than 1
        % N=3;              % Number of significant digits.  3 is the default
        % flag2=1;           % round to the nearest digit
        % mult=5;           % round to a multiple 5
        %
        % [P1, P1_str]=sd_round(D1, N, 1, 5);
        %
        % % P1_str should be 0.00000315 which has 3 significant digits.
        % % and the last digit is rounded to the nearest multiple of 5.
        %
        % Example3='';
        % N=4;                          % N is the number of significant digits
        % D2=10.^randn(10, 100);        % D2 is the unrounded array
        % [P2, P2_str]=sd_round(D2, N); % P2 is the rounded array
        %                               % of real numbers
        %                               % P2_str is the cell array of strings of
        %                               % rounded real numbers
        % Example4='';
        % D3=log10(randn(10, 100));     % D3 is an unrounded array of complex
        %                               % numbers
        % [P3, P3_str]=sd_round(D3, 4); % P3 is the rounded array of
        %                               % complex numbers
        %                               % P3_str is the cell array of strings of
        %                               % rounded complex numbers
        %
        % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %
        % % Program Written by Edward L. Zechmann
        % %
        % %      date 23 November 2007
        % %
        % %  modified 26 November 2007   updated comments
        % %
        % %  modified 17 December 2007   added outputs
        % %                              real_digitsL
        % %                              real_digitsR
        % %                              imag_digitsL
        % %                              imag_digitsR
        % %
        % %  modified 28 December 2007   added string output
        % %                              fixed bug in real_digitsR
        % %                              fixed bug in imag_digitsR
        % %
        % %  modified  7 January  2008   fixed bug in real_digitsR
        % %                              fixed bug in imag_digitsR
        % %                              sped up program by only
        % %                              converting the array to strings
        % %                              if it is an output argument
        % %
        % %  modified  1 March    2008   Added support for rounding
        % %                              to nearest integer, up, down,
        % %                              and toward zero
        % %
        % %  modified  3 March    2008   updated comments
        % %
        % %  modified 16 March    2008   Changed Program name from
        % %                              p_round to sd_round.
        % %                              Added another rounding style
        % %                              flag2 =5; (away from 0).
        % %                              Updated comments.
        % %
        % %  modified 18 August   2008   Added option to round last digit to a
        % %                              multiple of a given number.
        % %                              Fixed a bug in rounding powers of 10.
        % %                              Improved examples.
        % %
        % %  modified 21 August   2008   Fixed a bug in rounding numbers less
        % %                              than 1.  Added an example.
        % %
        % %  modified 25 August   2008   Modified program to recalculate the
        % %                              number of digits after rounding,
        % %                              because rounding can change the number
        % %                              of digits to the left and right of the
        % %                              decimal place. Updated Comments
        % %
        % % modified 16 August   2010   Changed K to k for SI prefixes.
        % %
        % %
        % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %
        % % Please Feel Free to Modify This Program
        % %
        % % See Also: pow10_round, round, ceil, floor, fix, fix2, round2, roundd
        % %
        flag1=0;
        if (nargin < 1 || isempty(A)) || ~isnumeric(A)
            flag1=1;
            A=[];
            A2=[];
            warning('Error in p_round did not Input Matrix A');
        end
        if isempty(A)
            flag1=1;
            A=[];
            A2=[];
            warning('Error in p_round Matrix A is Empty');
        end
        if (nargin < 2 || isempty(N)) || ~isnumeric(N)
            N=3;
        end
        if (nargin < 3 || isempty(flag2)) || ~isnumeric(flag2)
            flag2=1;
        end
        if (nargin < 4 || isempty(mult)) || ~isnumeric(mult)
            mult=1;
        end
        mult=round(mult);
        if mult < 1
            mult=1;
        end
        if (nargin < 5 || isempty(SI_prefixes)) || ~isnumeric(SI_prefixes)
            SI_prefixes=0;
        end
        SI_prefixes=SI_prefixes(1);
        % Number of digits to keep
        N=round(N);
        if N < 1
            N=1;
        end
        if ~isempty(A)
            real_digitsL=zeros(size(A));
            real_digitsR=zeros(size(A));
            imag_digitsL=zeros(size(A));
            imag_digitsR=zeros(size(A));
        else
            real_digitsL=[];
            real_digitsR=[];
            imag_digitsL=[];
            imag_digitsR=[];
        end
        letters={'y', 'z', 'a', 'f', 'p', 'n', '\mu', 'm', '', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'};
        if isequal(flag1, 0)
            if isnumeric(A)
                if isreal(A)
                    B=isinf(A);
                    A(B)=sign(A(B)).*10^15;
                    
                    
                    % Digit furthest to the left of the decimal point
                    D1=ceil(log10(abs(A)));
                    buf1=D1( abs(A)-10.^D1 == 0)+1;
                    D1( abs(A)-10.^D1 == 0)=buf1;
                    % rounding factor
                    dec=10.^(N-D1);
                    % Rounding Computation
                    % This program supports five different styles of rounding.
                    % flag2 == 1 rounds to the nearest integer
                    % flag2 == 2 rounds up
                    % flag2 == 3 rounds down
                    % flag2 == 4 rounds toward zero
                    % flag2 == 5 rounds away from zero
                    % otherwise round to the nearest integer
                    buf=dec./mult;
                    switch flag2
                        case 1
                            A2=1./buf.*round(buf.*A);
                        case 2
                            A2=1./buf.*ceil(buf.*A);
                        case 3
                            A2=1./buf.*floor(buf.*A);
                        case 4
                            A2=1./buf.*fix(buf.*A);
                        case 5
                            A2=sign(A)./buf.*ceil(buf.*abs(A));
                        otherwise
                            A2=1./buf.*round(buf.*A);
                    end
                    A2(A==0)=0;
                    if isequal( SI_prefixes, 1)
                        regime=floor(log10(abs(A2))./3);
                        regime(regime < -9)=0;
                        regime(regime > 9)=0;
                    else
                        regime=zeros(size(A));
                    end
                    A2=A2.*10.^(-3.*regime);
                    % After rounding recalculate the number of significant digits.
                    % The number of digits to the left and right of the decimal
                    % place can be effected by rounding.
                    % Digit furthest to the left of the decimal point
                    D1=ceil(log10(abs(A2)));
                    buf1=D1( abs(A2)-10.^D1 == 0)+1;
                    D1( abs(A2)-10.^D1 == 0)=buf1;
                    % Number of digits to the left of the decimal place
                    real_digitsL=max(D1, 0);
                    real_digitsL(A2==0)=0;
                    % rounding factor
                    dec=10.^(N-D1);
                    % Number of digits to the right of the decimal place
                    real_digitsR=max(N-D1, 0);
                    real_digitsR(A2==0)=N;
                else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Round the real part
                    Ar=real(A);
                    
                    B=isinf(Ar);
                    Ar(B)=sign(Ar(B)).*10^15;
                    
                    % Digit furthest to the left of the decimal point
                    D1=ceil(log10(abs(Ar)));
                    buf1=D1( abs(Ar)-10.^D1 == 0)+1;
                    D1( abs(Ar)-10.^D1 == 0)=buf1;
                    % rounding factor
                    dec=10.^(N-D1);
                    % Rounding Computation
                    % This program supports five different styles of rounding.
                    % flag2 == 1 rounds to the nearest integer
                    % flag2 == 2 rounds up
                    % flag2 == 3 rounds down
                    % flag2 == 4 rounds toward zero
                    % flag2 == 5 rounds away from zero
                    % otherwise round to the nearest integer
                    buf=dec./mult;
                    switch flag2
                        case 1
                            A2r=1./buf.*round(buf.*Ar);
                        case 2
                            A2r=1./buf.*ceil(buf.*Ar);
                        case 3
                            A2r=1./buf.*floor(buf.*Ar);
                        case 4
                            A2r=1./buf.*fix(buf.*Ar);
                        case 5
                            A2r=sign(Ar)./buf.*ceil(buf.*abs(Ar));
                        otherwise
                            A2r=1./buf.*round(buf.*Ar);
                    end
                    A2r(Ar==0)=0;
                    if isequal( SI_prefixes, 1)
                        regimeR=floor(log10(abs(A2r))./3);
                        regimeR(regimeR < -9)=0;
                        regimeR(regimeR > 9)=0;
                    else
                        regimeR=zeros(size(A));
                    end
                    A2r=A2r.*10.^(-3.*regimeR);
                    % After rounding recalculate the number of significant digits.
                    % The number of digits to the left and right of the decimal
                    % place can be effected by rounding.
                    % Digit furthest to the left of the decimal point
                    D1=ceil(log10(abs(A2r)));
                    buf1=D1( abs(A2r)-10.^D1 == 0)+1;
                    D1( abs(A2r)-10.^D1 == 0)=buf1;
                    % Number of digits to the left of the decimal place
                    real_digitsL=max(D1, 0);
                    real_digitsL(A2r==0)=0;
                    % rounding factor
                    dec=10.^(N-D1);
                    % Number of digits to the right of the decimal place
                    real_digitsR=max(N-D1, 0);
                    real_digitsR(A2r==0)=N;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Round the imaginary part
                    Ai=imag(A);
                    B=isinf(Ai);
                    Ai(B)=sign(Ai(B)).*10^15;
                    
                    
                    % Digit furthest to the left of the decimal point
                    D1=ceil(log10(abs(Ai)));
                    buf1=D1( abs(Ai)-10.^D1 == 0)+1;
                    D1( abs(Ai)-10.^D1 == 0)=buf1;
                    % rounding factor
                    dec=10.^(N-D1);
                    % Rounding Computation
                    % This program supports five different styles of rounding.
                    % flag2 == 1 rounds to the nearest integer
                    % flag2 == 2 rounds up
                    % flag2 == 3 rounds down
                    % flag2 == 4 rounds toward zero
                    % flag2 == 5 rounds away from zero
                    % otherwise round to the nearest integer
                    buf=dec./mult;
                    switch flag2
                        case 1
                            A2i=1./buf.*round(buf.*Ai);
                        case 2
                            A2i=1./buf.*ceil(buf.*Ai);
                        case 3
                            A2i=1./buf.*floor(buf.*Ai);
                        case 4
                            A2i=1./buf.*fix(buf.*Ai);
                        case 5
                            A2i=sign(Ai)./buf.*ceil(buf.*abs(Ai));
                        otherwise
                            A2i=1./buf.*round(buf.*Ai);
                    end
                    A2i(Ai==0)=0;
                    if isequal( SI_prefixes, 1)
                        regimeI=floor(log10(abs(A2i))./3);
                        regimeI(regimeI < -9)=0;
                        regimeI(regimeI > 9)=0;
                    else
                        regimeI=zeros(size(A));
                    end
                    A2i=A2i.*10.^(-3.*regimeI);
                    % After rounding recalculate the number of significant digits.
                    % The number of digits to the left and right of the decimal
                    % place can be effected by rounding.
                    % Digit furthest to the left of the decimal point
                    D1=ceil(log10(abs(A2i)));
                    buf1=D1( abs(A2i)-10.^D1 == 0)+1;
                    D1( abs(A2i)-10.^D1 == 0)=buf1;
                    % Number of digits to the left of the decimal place
                    imag_digitsL=max(D1, 0);
                    imag_digitsL(A2i==0)=0;
                    % rounding factor
                    dec=10.^(N-D1);
                    % Number of digits to the right of the decimal place
                    imag_digitsR=max(N-D1, 0);
                    imag_digitsR(A2i==0)=N;
                end
            else
                warningdlg('Error in sd_round Input Matrix A is not numeric');
                A2=A;
            end
        end
        % % Convert the rounded array to string format with specified
        % % number of significant digits.
        [m1, n1]=size(A);
        A_str=cell(size(A));
        real_digitsL(isinf(real_digitsL))=15;
        real_digitsR(isinf(real_digitsR))=15;
        imag_digitsL(isinf(imag_digitsL))=15;
        imag_digitsR(isinf(imag_digitsR))=15;
        rtd=round(real_digitsL+real_digitsR);
        itd=round(imag_digitsL+imag_digitsR);
        % Output the cell array of strings if it is in the output argument list
        if nargout > 1 && isequal(flag1, 0)
            % This code formats the rounded numbers into a cell array
            % of strings.
            if isreal(A)
                for e1=1:m1
                    for e2=1:n1;
                        aa2=num2str(A2(e1, e2), ['%', int2str(rtd(e1, e2)), '.', int2str(real_digitsR(e1, e2)), 'f' ]);
                        if length(aa2) < N && ~isequal( SI_prefixes, 1)
                            aa2=num2str(A2(e1, e2),['%', int2str(N), '.', int2str(N), 'f']);
                        end
                        if isequal( SI_prefixes, 1)
                            A_str{e1, e2}=[aa2 letters{regime(e1, e2)+9}];
                        else
                            A_str{e1, e2}=aa2;
                        end
                    end
                end
                A2=A2.*10.^(3.*regime);
            else
                for e1=1:m1
                    for e2=1:n1;
                        aa1=num2str(A2r(e1, e2), ['%', int2str(rtd(e1, e2)),'.', int2str(real_digitsR(e1, e2)), 'f' ]);
                        if length(aa1) < N && ~isequal( SI_prefixes, 1)
                            aa1=num2str(A2r(e1, e2),['%', int2str(N), '.', int2str(N), 'f']);
                        end
                        aa2=num2str(abs(A2i(e1, e2)), ['%', int2str(itd(e1, e2)),'.', int2str(imag_digitsR(e1, e2)), 'f' ]);
                        if length(aa2) < N && ~isequal( SI_prefixes, 1)
                            aa2=num2str(abs(A2i(e1, e2)),['%', int2str(N), '.', int2str(N), 'f']);
                        end
                        if imag(A2i(e1, e2)) >= 0
                            aa=' + ';
                        else
                            aa=' - ';
                        end
                        if isequal( SI_prefixes, 1)
                            A_str{e1, e2}= [aa1, letters{regimeR(e1, e2)+9} aa, aa2, letters{regimeI(e1, e2)+9} 'i'];
                        else
                            A_str{e1, e2}= [aa1, aa, aa2, 'i'];
                        end
                    end
                end
                % Add the real and imaginary parts together
                A2=A2r.*10.^(3.*regimeR)+i*A2i.*10.^(3.*regimeI);
            end
        else
            A_str={};
        end

    end

%##############################################################################################
%  Auxiliary function -- oct3dsgn from octave toolbox.
%  octave, by Christophe COUVREUR, 29 Dec 1997
%  from http://www.mathworks.com/matlabcentral/fileexchange/69-octave
%##############################################################################################

function [B,A] = oct3dsgn(Fc,Fs,N); 
% OCT3DSGN  Design of a one-third-octave filter.
%    [B,A] = OCT3DSGN(Fc,Fs,N) designs a digital 1/3-octave filter with 
%    center frequency Fc for sampling frequency Fs. 
%    The filter is designed according to the Order-N specification 
%    of the ANSI S1.1-1986 standard. Default value for N is 3. 
%    Warning: for meaningful design results, center frequency used
%    should preferably be in range Fs/200 < Fc < Fs/5.
%    Usage of the filter: Y = FILTER(B,A,X). 
%
%    Requires the Signal Processing Toolbox. 
%
%    See also OCT3SPEC, OCTDSGN, OCTSPEC.

% Author: Christophe Couvreur, Faculte Polytechnique de Mons (Belgium)
%         couvreur@thor.fpms.ac.be
% Last modification: Aug. 25, 1997, 2:00pm.

% References: 
%    [1] ANSI S1.1-1986 (ASA 65-1986): Specifications for
%        Octave-Band and Fractional-Octave-Band Analog and
%        Digital Filters, 1993.

if (nargin > 3) | (nargin < 2)
  error('Invalid number of arguments.');
end
if (nargin == 2)
  N = 3; 
end
if (Fc > 0.88*(Fs/2))
  error('Design not possible. Check frequencies.');
end
  
% Design Butterworth 2Nth-order one-third-octave filter 
% Note: BUTTER is based on a bilinear transformation, as suggested in [1]. 
pi = 3.14159265358979;
f1 = Fc/(2^(1/6)); 
f2 = Fc*(2^(1/6)); 
Qr = Fc/(f2-f1); 
Qd = (pi/2/N)/(sin(pi/2/N))*Qr;
alpha = (1 + sqrt(1+4*Qd^2))/2/Qd; 
W1 = Fc/(Fs/2)/alpha; 
W2 = Fc/(Fs/2)*alpha;
[B,A] = butter(N,[W1,W2]); 
end

function imageData = screencapture(varargin)
% screencapture - get a screen-capture of a figure frame, component handle, or screen area rectangle
%
% ScreenCapture gets a screen-capture of any Matlab GUI handle (including desktop, 
% figure, axes, image or uicontrol), or a specified area rectangle located relative to
% the specified handle. Screen area capture is possible by specifying the root (desktop)
% handle (=0). The output can be either to an image file or to a Matlab matrix (useful
% for displaying via imshow() or for further processing) or to the system clipboard.
% This utility also enables adding a toolbar button for easy interactive screen-capture.
%
% Syntax:
%    imageData = screencapture(handle, position, target, 'PropName',PropValue, ...)
%
% Input Parameters:
%    handle   - optional handle to be used for screen-capture origin.
%                 If empty/unsupplied then current figure (gcf) will be used.
%    position - optional position array in pixels: [x,y,width,height].
%                 If empty/unsupplied then the handle's position vector will be used.
%                 If both handle and position are empty/unsupplied then the position
%                   will be retrieved via interactive mouse-selection.
%                 If handle is an image, then position is in data (not pixel) units, so the
%                   captured region remains the same after figure/axes resize (like imcrop)
%    target   - optional filename for storing the screen-capture, or the
%               'clipboard'/'printer' strings.
%                 If empty/unsupplied then no output to file will be done.
%                 The file format will be determined from the extension (JPG/PNG/...).
%                 Supported formats are those supported by the imwrite function.
%    'PropName',PropValue - 
%               optional list of property pairs (e.g., 'target','myImage.png','pos',[10,20,30,40],'handle',gca)
%               PropNames may be abbreviated and are case-insensitive.
%               PropNames may also be given in whichever order.
%               Supported PropNames are:
%                 - 'handle'    (default: gcf handle)
%                 - 'position'  (default: gcf position array)
%                 - 'target'    (default: '')
%                 - 'toolbar'   (figure handle; default: gcf)
%                      this adds a screen-capture button to the figure's toolbar
%                      If this parameter is specified, then no screen-capture
%                        will take place and the returned imageData will be [].
%
% Output parameters:
%    imageData - image data in a format acceptable by the imshow function
%                  If neither target nor imageData were specified, the user will be
%                    asked to interactively specify the output file.
%
% Examples:
%    imageData = screencapture;  % interactively select screen-capture rectangle
%    imageData = screencapture(hListbox);  % capture image of a uicontrol
%    imageData = screencapture(0,  [20,30,40,50]);  % capture a small desktop region
%    imageData = screencapture(gcf,[20,30,40,50]);  % capture a small figure region
%    imageData = screencapture(gca,[10,20,30,40]);  % capture a small axes region
%      imshow(imageData);  % display the captured image in a matlab figure
%      imwrite(imageData,'myImage.png');  % save the captured image to file
%    img = imread('cameraman.tif');
%      hImg = imshow(img);
%      screencapture(hImg,[60,35,140,80]);  % capture a region of an image
%    screencapture(gcf,[],'myFigure.jpg');  % capture the entire figure into file
%    screencapture(gcf,[],'clipboard');     % capture the entire figure into clipboard
%    screencapture(gcf,[],'printer');       % print the entire figure
%    screencapture('handle',gcf,'target','myFigure.jpg'); % same as previous, save to file
%    screencapture('handle',gcf,'target','clipboard');    % same as previous, copy to clipboard
%    screencapture('handle',gcf,'target','printer');      % same as previous, send to printer
%    screencapture('toolbar',gcf);  % adds a screen-capture button to gcf's toolbar
%    screencapture('toolbar',[],'target','sc.bmp'); % same with default output filename
%
% Technical description:
%    http://UndocumentedMatlab.com/blog/screencapture-utility/
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% See also:
%    imshow, imwrite, print
%
% Release history:
%    1.17 2016-05-16: Fix annoying warning about JavaFrame property becoming obsolete someday (yes, we know...)
%    1.16 2016-04-19: Fix for deployed application suggested by Dwight Bartholomew
%    1.10 2014-11-25: Added the 'print' target
%    1.9  2014-11-25: Fix for saving GIF files
%    1.8  2014-11-16: Fixes for R2014b
%    1.7  2014-04-28: Fixed bug when capturing interactive selection
%    1.6  2014-04-22: Only enable image formats when saving to an unspecified file via uiputfile
%    1.5  2013-04-18: Fixed bug in capture of non-square image; fixes for Win64
%    1.4  2013-01-27: Fixed capture of Desktop (root); enabled rbbox anywhere on desktop (not necesarily in a Matlab figure); enabled output to clipboard (based on Jiro Doke's imclipboard utility); edge-case fixes; added Java compatibility check
%    1.3  2012-07-23: Capture current object (uicontrol/axes/figure) if w=h=0 (e.g., by clicking a single point); extra input args sanity checks; fix for docked windows and image axes; include axes labels & ticks by default when capturing axes; use data-units position vector when capturing images; many edge-case fixes
%    1.2  2011-01-16: another performance boost (thanks to Jan Simon); some compatibility fixes for Matlab 6.5 (untested)
%    1.1  2009-06-03: Handle missing output format; performance boost (thanks to Urs); fix minor root-handle bug; added toolbar button option
%    1.0  2009-06-02: First version posted on <a href="http://www.mathworks.com/matlabcentral/fileexchange/authors/27420">MathWorks File Exchange</a>
% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.
% Programmed and Copyright by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.17 $  $Date: 2016/05/16 17:59:36 $
    % Ensure that java awt is enabled...
    if ~usejava('awt')
        error('YMA:screencapture:NeedAwt','ScreenCapture requires Java to run.');
    end
    % Ensure that our Java version supports the Robot class (requires JVM 1.3+)
    try
        robot = java.awt.Robot; %#ok<NASGU>
    catch
        uiwait(msgbox({['Your Matlab installation is so old that its Java engine (' version('-java') ...
                        ') does not have a java.awt.Robot class. '], ' ', ...
                        'Without this class, taking a screen-capture is impossible.', ' ', ...
                        'So, either install JVM 1.3 or higher, or use a newer Matlab release.'}, ...
                        'ScreenCapture', 'warn'));
        if nargout, imageData = [];  end
        return;
    end
    % Process optional arguments
    paramsStruct = processArgs(varargin{:});
    
    % If toolbar button requested, add it and exit
    if ~isempty(paramsStruct.toolbar)
        
        % Add the toolbar button
        addToolbarButton(paramsStruct);
        
        % Return the figure to its pre-undocked state (when relevant)
        redockFigureIfRelevant(paramsStruct);
        
        % Exit immediately (do NOT take a screen-capture)
        if nargout,  imageData = [];  end
        return;
    end
    
    % Convert position from handle-relative to desktop Java-based pixels
    [paramsStruct, msgStr] = convertPos(paramsStruct);
    
    % Capture the requested screen rectangle using java.awt.Robot
    imgData = getScreenCaptureImageData(paramsStruct.position);
    
    % Return the figure to its pre-undocked state (when relevant)
    redockFigureIfRelevant(paramsStruct);
    
    % Save image data in file or clipboard, if specified
    if ~isempty(paramsStruct.target)
        if strcmpi(paramsStruct.target,'clipboard')
            if ~isempty(imgData)
                imclipboard(imgData);
            else
                msgbox('No image area selected - not copying image to clipboard','ScreenCapture','warn');
            end
        elseif strncmpi(paramsStruct.target,'print',5)  % 'print' or 'printer'
            if ~isempty(imgData)
                hNewFig = figure('visible','off');
                imshow(imgData);
                print(hNewFig);
                delete(hNewFig);
            else
                msgbox('No image area selected - not printing screenshot','ScreenCapture','warn');
            end
        else  % real filename
            if ~isempty(imgData)
                imwrite(imgData,paramsStruct.target);
            else
                msgbox(['No image area selected - not saving image file ' paramsStruct.target],'ScreenCapture','warn');
            end
        end
    end
    % Return image raster data to user, if requested
    if nargout
        imageData = imgData;
        
    % If neither output formats was specified (neither target nor output data)
    elseif isempty(paramsStruct.target) & ~isempty(imgData)  %#ok ML6
        % Ask the user to specify a file
        %error('YMA:screencapture:noOutput','No output specified for ScreenCapture: specify the output filename and/or output data');
        %format = '*.*';
        formats = imformats;
        for idx = 1 : numel(formats)
            ext = sprintf('*.%s;',formats(idx).ext{:});
            format(idx,1:2) = {ext(1:end-1), formats(idx).description}; %#ok<AGROW>
        end
        [filename,pathname] = uiputfile(format,'Save screen capture as');
        if ~isequal(filename,0) & ~isequal(pathname,0)  %#ok Matlab6 compatibility
            try
                filename = fullfile(pathname,filename);
                imwrite(imgData,filename);
            catch  % possibly a GIF file that requires indexed colors
                [imgData,map] = rgb2ind(imgData,256);
                imwrite(imgData,map,filename);
            end
        else
            % TODO - copy to clipboard
        end
    end
    % Display msgStr, if relevant
    if ~isempty(msgStr)
        uiwait(msgbox(msgStr,'ScreenCapture'));
        drawnow; pause(0.05);  % time for the msgbox to disappear
    end
    return;  % debug breakpoint
%% Process optional arguments
function paramsStruct = processArgs(varargin)
    % Get the properties in either direct or P-V format
    [regParams, pvPairs] = parseparams(varargin);
    % Now process the optional P-V params
    try
        % Initialize
        paramName = [];
        paramsStruct = [];
        paramsStruct.handle = [];
        paramsStruct.position = [];
        paramsStruct.target = '';
        paramsStruct.toolbar = [];
        paramsStruct.wasDocked = 0;       % no false available in ML6
        paramsStruct.wasInteractive = 0;  % no false available in ML6
        % Parse the regular (non-named) params in recption order
        if ~isempty(regParams) & (isempty(regParams{1}) | ishandle(regParams{1}(1)))  %#ok ML6
            paramsStruct.handle = regParams{1};
            regParams(1) = [];
        end
        if ~isempty(regParams) & isnumeric(regParams{1}) & (length(regParams{1}) == 4)  %#ok ML6
            paramsStruct.position = regParams{1};
            regParams(1) = [];
        end
        if ~isempty(regParams) & ischar(regParams{1})  %#ok ML6
            paramsStruct.target = regParams{1};
        end
        % Parse the optional param PV pairs
        supportedArgs = {'handle','position','target','toolbar'};
        while ~isempty(pvPairs)
            % Disregard empty propNames (may be due to users mis-interpretting the syntax help)
            while ~isempty(pvPairs) & isempty(pvPairs{1})  %#ok ML6
                pvPairs(1) = [];
            end
            if isempty(pvPairs)
                break;
            end
            % Ensure basic format is valid
            paramName = '';
            if ~ischar(pvPairs{1})
                error('YMA:screencapture:invalidProperty','Invalid property passed to ScreenCapture');
            elseif length(pvPairs) == 1
                if isempty(paramsStruct.target)
                    paramsStruct.target = pvPairs{1};
                    break;
                else
                    error('YMA:screencapture:noPropertyValue',['No value specified for property ''' pvPairs{1} '''']);
                end
            end
            % Process parameter values
            paramName  = pvPairs{1};
            if strcmpi(paramName,'filename')  % backward compatibility
                paramName = 'target';
            end
            paramValue = pvPairs{2};
            pvPairs(1:2) = [];
            idx = find(strncmpi(paramName,supportedArgs,length(paramName)));
            if ~isempty(idx)
                %paramsStruct.(lower(supportedArgs{idx(1)})) = paramValue;  % incompatible with ML6
                paramsStruct = setfield(paramsStruct, lower(supportedArgs{idx(1)}), paramValue);  %#ok ML6
                % If 'toolbar' param specified, then it cannot be left empty - use gcf
                if strncmpi(paramName,'toolbar',length(paramName)) & isempty(paramsStruct.toolbar)  %#ok ML6
                    paramsStruct.toolbar = getCurrentFig;
                end
            elseif isempty(paramsStruct.target)
                paramsStruct.target = paramName;
                pvPairs = {paramValue, pvPairs{:}};  %#ok (more readable this way, although a bit less efficient...)
            else
                supportedArgsStr = sprintf('''%s'',',supportedArgs{:});
                error('YMA:screencapture:invalidProperty','%s \n%s', ...
                      'Invalid property passed to ScreenCapture', ...
                      ['Supported property names are: ' supportedArgsStr(1:end-1)]);
            end
        end  % loop pvPairs
    catch
        if ~isempty(paramName),  paramName = [' ''' paramName ''''];  end
        error('YMA:screencapture:invalidProperty','Error setting ScreenCapture property %s:\n%s',paramName,lasterr); %#ok<LERR>
    end
end  % processArgs
%% Convert position from handle-relative to desktop Java-based pixels
function [paramsStruct, msgStr] = convertPos(paramsStruct)
    msgStr = '';
    try
        % Get the screen-size for later use
        screenSize = get(0,'ScreenSize');
        % Get the containing figure's handle
        hParent = paramsStruct.handle;
        if isempty(paramsStruct.handle)
            paramsStruct.hFigure = getCurrentFig;
            hParent = paramsStruct.hFigure;
        else
            paramsStruct.hFigure = ancestor(paramsStruct.handle,'figure');
        end
        % To get the acurate pixel position, the figure window must be undocked
        try
            if strcmpi(get(paramsStruct.hFigure,'WindowStyle'),'docked')
                set(paramsStruct.hFigure,'WindowStyle','normal');
                drawnow; pause(0.25);
                paramsStruct.wasDocked = 1;  % no true available in ML6
            end
        catch
            % never mind - ignore...
        end
        % The figure (if specified) must be in focus
        if ~isempty(paramsStruct.hFigure) & ishandle(paramsStruct.hFigure)  %#ok ML6
            isFigureValid = 1;  % no true available in ML6
            figure(paramsStruct.hFigure);
        else
            isFigureValid = 0;  % no false available in ML6
        end
        % Flush all graphic events to ensure correct rendering
        drawnow; pause(0.01);
        % No handle specified
        wasPositionGiven = 1;  % no true available in ML6
        if isempty(paramsStruct.handle)
            
            % Set default handle, if not supplied
            paramsStruct.handle = paramsStruct.hFigure;
            
            % If position was not specified, get it interactively using RBBOX
            if isempty(paramsStruct.position)
                [paramsStruct.position, jFrameUsed, msgStr] = getInteractivePosition(paramsStruct.hFigure); %#ok<ASGLU> jFrameUsed is unused
                paramsStruct.wasInteractive = 1;  % no true available in ML6
                wasPositionGiven = 0;  % no false available in ML6
            end
            
        elseif ~ishandle(paramsStruct.handle)
            % Handle was supplied - ensure it is a valid handle
            error('YMA:screencapture:invalidHandle','Invalid handle passed to ScreenCapture');
            
        elseif isempty(paramsStruct.position)
            % Handle was supplied but position was not, so use the handle's position
            paramsStruct.position = getPixelPos(paramsStruct.handle);
            paramsStruct.position(1:2) = 0;
            wasPositionGiven = 0;  % no false available in ML6
            
        elseif ~isnumeric(paramsStruct.position) | (length(paramsStruct.position) ~= 4)  %#ok ML6
            % Both handle & position were supplied - ensure a valid pixel position vector
            error('YMA:screencapture:invalidPosition','Invalid position vector passed to ScreenCapture: \nMust be a [x,y,w,h] numeric pixel array');
        end
        
        % Capture current object (uicontrol/axes/figure) if w=h=0 (single-click in interactive mode)
        if paramsStruct.position(3)<=0 | paramsStruct.position(4)<=0  %#ok ML6
            %TODO - find a way to single-click another Matlab figure (the following does not work)
            %paramsStruct.position = getPixelPos(ancestor(hittest,'figure'));
            paramsStruct.position = getPixelPos(paramsStruct.handle);
            paramsStruct.position(1:2) = 0;
            paramsStruct.wasInteractive = 0;  % no false available in ML6
            wasPositionGiven = 0;  % no false available in ML6
        end
        % First get the parent handle's desktop-based Matlab pixel position
        parentPos = [0,0,0,0];
        dX = 0;
        dY = 0;
        dW = 0;
        dH = 0;
        if ~isFigure(hParent)
            % Get the reguested component's pixel position
            parentPos = getPixelPos(hParent, 1);  % no true available in ML6
            % Axes position inaccuracy estimation
            deltaX = 3;
            deltaY = -1;
            
            % Fix for images
            if isImage(hParent)  % | (isAxes(hParent) & strcmpi(get(hParent,'YDir'),'reverse'))  %#ok ML6
                % Compensate for resized image axes
                hAxes = get(hParent,'Parent');
                if all(get(hAxes,'DataAspectRatio')==1)  % sanity check: this is the normal behavior
                    % Note 18/4/2013: the following fails for non-square images
                    %actualImgSize = min(parentPos(3:4));
                    %dX = (parentPos(3) - actualImgSize) / 2;
                    %dY = (parentPos(4) - actualImgSize) / 2;
                    %parentPos(3:4) = actualImgSize;
                    % The following should work for all types of images
                    actualImgSize = size(get(hParent,'CData'));
                    dX = (parentPos(3) - min(parentPos(3),actualImgSize(2))) / 2;
                    dY = (parentPos(4) - min(parentPos(4),actualImgSize(1))) / 2;
                    parentPos(3:4) = actualImgSize([2,1]);
                    %parentPos(3) = max(parentPos(3),actualImgSize(2));
                    %parentPos(4) = max(parentPos(4),actualImgSize(1));
                end
                % Fix user-specified img positions (but not auto-inferred ones)
                if wasPositionGiven
                    % In images, use data units rather than pixel units
                    % Reverse the YDir
                    ymax = max(get(hParent,'YData'));
                    paramsStruct.position(2) = ymax - paramsStruct.position(2) - paramsStruct.position(4);
                    % Note: it would be best to use hgconvertunits, but:
                    % ^^^^  (1) it fails on Matlab 6, and (2) it doesn't accept Data units
                    %paramsStruct.position = hgconvertunits(hFig, paramsStruct.position, 'Data', 'pixel', hParent);  % fails!
                    xLims = get(hParent,'XData');
                    yLims = get(hParent,'YData');
                    xPixelsPerData = parentPos(3) / (diff(xLims) + 1);
                    yPixelsPerData = parentPos(4) / (diff(yLims) + 1);
                    paramsStruct.position(1) = round((paramsStruct.position(1)-xLims(1)) * xPixelsPerData);
                    paramsStruct.position(2) = round((paramsStruct.position(2)-yLims(1)) * yPixelsPerData + 2*dY);
                    paramsStruct.position(3) = round( paramsStruct.position(3) * xPixelsPerData);
                    paramsStruct.position(4) = round( paramsStruct.position(4) * yPixelsPerData);
                    % Axes position inaccuracy estimation
                    if strcmpi(computer('arch'),'win64')
                        deltaX = 7;
                        deltaY = -7;
                    else
                        deltaX = 3;
                        deltaY = -3;
                    end
                    
                else  % axes/image position was auto-infered (entire image)
                    % Axes position inaccuracy estimation
                    if strcmpi(computer('arch'),'win64')
                        deltaX = 6;
                        deltaY = -6;
                    else
                        deltaX = 2;
                        deltaY = -2;
                    end
                    dW = -2*dX;
                    dH = -2*dY;
                end
            end
            %hFig = ancestor(hParent,'figure');
            hParent = paramsStruct.hFigure;
        elseif paramsStruct.wasInteractive  % interactive figure rectangle
            % Compensate for 1px rbbox inaccuracies
            deltaX = 2;
            deltaY = -2;
        else  % non-interactive figure
            % Compensate 4px figure boundaries = difference betweeen OuterPosition and Position
            deltaX = -1;
            deltaY = 1;
        end
        %disp(paramsStruct.position)  % for debugging
        
        % Now get the pixel position relative to the monitor
        figurePos = getPixelPos(hParent);
        desktopPos = figurePos + parentPos;
        % Now convert to Java-based pixels based on screen size
        % Note: multiple monitors are automatically handled correctly, since all
        % ^^^^  Java positions are relative to the main monitor's top-left corner
        javaX  = desktopPos(1) + paramsStruct.position(1) + deltaX + dX;
        javaY  = screenSize(4) - desktopPos(2) - paramsStruct.position(2) - paramsStruct.position(4) + deltaY + dY;
        width  = paramsStruct.position(3) + dW;
        height = paramsStruct.position(4) + dH;
        paramsStruct.position = round([javaX, javaY, width, height]);
        %paramsStruct.position
        % Ensure the figure is at the front so it can be screen-captured
        if isFigureValid
            figure(hParent);
            drawnow;
            pause(0.02);
        end
    catch
        % Maybe root/desktop handle (root does not have a 'Position' prop so getPixelPos croaks
        if isequal(double(hParent),0)  % =root/desktop handle;  handles case of hParent=[]
            javaX = paramsStruct.position(1) - 1;
            javaY = screenSize(4) - paramsStruct.position(2) - paramsStruct.position(4) - 1;
            paramsStruct.position = [javaX, javaY, paramsStruct.position(3:4)];
        end
    end
end  % convertPos
%% Interactively get the requested capture rectangle
function [positionRect, jFrameUsed, msgStr] = getInteractivePosition(hFig)
    msgStr = '';
    try
        % First try the invisible-figure approach, in order to
        % enable rbbox outside any existing figure boundaries
        f = figure('units','pixel','pos',[-100,-100,10,10],'HitTest','off');
        drawnow; pause(0.01);
        oldWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
        jf = get(handle(f),'JavaFrame');
        warning(oldWarn);
        try
            jWindow = jf.fFigureClient.getWindow;
        catch
            try
                jWindow = jf.fHG1Client.getWindow;
            catch
                jWindow = jf.getFigurePanelContainer.getParent.getTopLevelAncestor;
            end
        end
        com.sun.awt.AWTUtilities.setWindowOpacity(jWindow,0.05);  %=nearly transparent (not fully so that mouse clicks are captured)
        jWindow.setMaximized(1);  % no true available in ML6
        jFrameUsed = 1;  % no true available in ML6
        msg = {'Mouse-click and drag a bounding rectangle for screen-capture ' ...
               ... %'or single-click any Matlab figure to capture the entire figure.' ...
               };
    catch
        % Something failed, so revert to a simple rbbox on a visible figure
        try delete(f); drawnow; catch, end  %Cleanup...
        jFrameUsed = 0;  % no false available in ML6
        msg = {'Mouse-click within any Matlab figure and then', ...
               'drag a bounding rectangle for screen-capture,', ...
               'or single-click to capture the entire figure'};
    end
    uiwait(msgbox(msg,'ScreenCapture'));
    
    k = waitforbuttonpress;  %#ok k is unused
    %hFig = getCurrentFig;
    %p1 = get(hFig,'CurrentPoint');
    positionRect = rbbox;
    %p2 = get(hFig,'CurrentPoint');
    if jFrameUsed
        jFrameOrigin = getPixelPos(f);
        delete(f); drawnow;
        try
            figOrigin = getPixelPos(hFig);
        catch  % empty/invalid hFig handle
            figOrigin = [0,0,0,0];
        end
    else
        if isempty(hFig)
            jFrameOrigin = getPixelPos(gcf);
        else
            jFrameOrigin = [0,0,0,0];
        end
        figOrigin = [0,0,0,0];
    end
    positionRect(1:2) = positionRect(1:2) + jFrameOrigin(1:2) - figOrigin(1:2);
    if prod(positionRect(3:4)) > 0
        msgStr = sprintf('%dx%d area captured',positionRect(3),positionRect(4));
    end
end  % getInteractivePosition
%% Get current figure (even if its handle is hidden)
function hFig = getCurrentFig
    oldState = get(0,'showHiddenHandles');
    set(0,'showHiddenHandles','on');
    hFig = get(0,'CurrentFigure');
    set(0,'showHiddenHandles',oldState);
end  % getCurrentFig
%% Get ancestor figure - used for old Matlab versions that don't have a built-in ancestor()
function hObj = ancestor(hObj,type)
    if ~isempty(hObj) & ishandle(hObj)  %#ok for Matlab 6 compatibility
        try
            hObj = get(hObj,'Ancestor');
        catch
            % never mind...
        end
        try
            %if ~isa(handle(hObj),type)  % this is best but always returns 0 in Matlab 6!
            %if ~isprop(hObj,'type') | ~strcmpi(get(hObj,'type'),type)  % no isprop() in ML6!
            try
                objType = get(hObj,'type');
            catch
                objType = '';
            end
            if ~strcmpi(objType,type)
                try
                    parent = get(handle(hObj),'parent');
                catch
                    parent = hObj.getParent;  % some objs have no 'Parent' prop, just this method...
                end
                if ~isempty(parent)  % empty parent means root ancestor, so exit
                    hObj = ancestor(parent,type);
                end
            end
        catch
            % never mind...
        end
    end
end  % ancestor
%% Get position of an HG object in specified units
function pos = getPos(hObj,field,units)
    % Matlab 6 did not have hgconvertunits so use the old way...
    oldUnits = get(hObj,'units');
    if strcmpi(oldUnits,units)  % don't modify units unless we must!
        pos = get(hObj,field);
    else
        set(hObj,'units',units);
        pos = get(hObj,field);
        set(hObj,'units',oldUnits);
    end
end  % getPos
%% Get pixel position of an HG object - for Matlab 6 compatibility
function pos = getPixelPos(hObj,varargin)
    persistent originalObj
    try
        stk = dbstack;
        if ~strcmp(stk(2).name,'getPixelPos')
            originalObj = hObj;
        end
        if isFigure(hObj) %| isAxes(hObj)
        %try
            pos = getPos(hObj,'OuterPosition','pixels');
        else  %catch
            % getpixelposition is unvectorized unfortunately!
            pos = getpixelposition(hObj,varargin{:});
            % add the axes labels/ticks if relevant (plus a tiny margin to fix 2px label/title inconsistencies)
            if isAxes(hObj) & ~isImage(originalObj)  %#ok ML6
                tightInsets = getPos(hObj,'TightInset','pixel');
                pos = pos + tightInsets.*[-1,-1,1,1] + [-1,1,1+tightInsets(1:2)];
            end
        end
    catch
        try
            % Matlab 6 did not have getpixelposition nor hgconvertunits so use the old way...
            pos = getPos(hObj,'Position','pixels');
        catch
            % Maybe the handle does not have a 'Position' prop (e.g., text/line/plot) - use its parent
            pos = getPixelPos(get(hObj,'parent'),varargin{:});
        end
    end
    % Handle the case of missing/invalid/empty HG handle
    if isempty(pos)
        pos = [0,0,0,0];
    end
end  % getPixelPos
%% Adds a ScreenCapture toolbar button
function addToolbarButton(paramsStruct)
    % Ensure we have a valid toolbar handle
    hFig = ancestor(paramsStruct.toolbar,'figure');
    if isempty(hFig)
        error('YMA:screencapture:badToolbar','the ''Toolbar'' parameter must contain a valid GUI handle');
    end
    set(hFig,'ToolBar','figure');
    hToolbar = findall(hFig,'type','uitoolbar');
    if isempty(hToolbar)
        error('YMA:screencapture:noToolbar','the ''Toolbar'' parameter must contain a figure handle possessing a valid toolbar');
    end
    hToolbar = hToolbar(1);  % just in case there are several toolbars... - use only the first
    % Prepare the camera icon
    icon = ['3333333333333333'; ...
            '3333333333333333'; ...
            '3333300000333333'; ...
            '3333065556033333'; ...
            '3000000000000033'; ...
            '3022222222222033'; ...
            '3022220002222033'; ...
            '3022203110222033'; ...
            '3022201110222033'; ...
            '3022204440222033'; ...
            '3022220002222033'; ...
            '3022222222222033'; ...
            '3000000000000033'; ...
            '3333333333333333'; ...
            '3333333333333333'; ...
            '3333333333333333'];
    cm = [   0      0      0; ...  % black
             0   0.60      1; ...  % light blue
          0.53   0.53   0.53; ...  % light gray
           NaN    NaN    NaN; ...  % transparent
             0   0.73      0; ...  % light green
          0.27   0.27   0.27; ...  % gray
          0.13   0.13   0.13];     % dark gray
    cdata = ind2rgb(uint8(icon-'0'),cm);
    % If the button does not already exit
    hButton = findall(hToolbar,'Tag','ScreenCaptureButton');
    tooltip = 'Screen capture';
    if ~isempty(paramsStruct.target)
        tooltip = [tooltip ' to ' paramsStruct.target];
    end
    if isempty(hButton)
        % Add the button with the icon to the figure's toolbar
        hButton = uipushtool(hToolbar, 'CData',cdata, 'Tag','ScreenCaptureButton', 'TooltipString',tooltip, 'ClickedCallback',['screencapture(''' paramsStruct.target ''')']);  %#ok unused
    else
        % Otherwise, simply update the existing button
        set(hButton, 'CData',cdata, 'Tag','ScreenCaptureButton', 'TooltipString',tooltip, 'ClickedCallback',['screencapture(''' paramsStruct.target ''')']);
    end
end  % addToolbarButton
%% Java-get the actual screen-capture image data
function imgData = getScreenCaptureImageData(positionRect)
    if isempty(positionRect) | all(positionRect==0) | positionRect(3)<=0 | positionRect(4)<=0  %#ok ML6
        imgData = [];
    else
        % Use java.awt.Robot to take a screen-capture of the specified screen area
        rect = java.awt.Rectangle(positionRect(1), positionRect(2), positionRect(3), positionRect(4));
        robot = java.awt.Robot;
        jImage = robot.createScreenCapture(rect);
        % Convert the resulting Java image to a Matlab image
        % Adapted for a much-improved performance from:
        % http://www.mathworks.com/support/solutions/data/1-2WPAYR.html
        h = jImage.getHeight;
        w = jImage.getWidth;
        %imgData = zeros([h,w,3],'uint8');
        %pixelsData = uint8(jImage.getData.getPixels(0,0,w,h,[]));
        %for i = 1 : h
        %    base = (i-1)*w*3+1;
        %    imgData(i,1:w,:) = deal(reshape(pixelsData(base:(base+3*w-1)),3,w)');
        %end
        % Performance further improved based on feedback from Urs Schwartz:
        %pixelsData = reshape(typecast(jImage.getData.getDataStorage,'uint32'),w,h).';
        %imgData(:,:,3) = bitshift(bitand(pixelsData,256^1-1),-8*0);
        %imgData(:,:,2) = bitshift(bitand(pixelsData,256^2-1),-8*1);
        %imgData(:,:,1) = bitshift(bitand(pixelsData,256^3-1),-8*2);
        % Performance even further improved based on feedback from Jan Simon:
        pixelsData = reshape(typecast(jImage.getData.getDataStorage, 'uint8'), 4, w, h);
        imgData = cat(3, ...
            transpose(reshape(pixelsData(3, :, :), w, h)), ...
            transpose(reshape(pixelsData(2, :, :), w, h)), ...
            transpose(reshape(pixelsData(1, :, :), w, h)));
    end
end  % getInteractivePosition
%% Return the figure to its pre-undocked state (when relevant)
function redockFigureIfRelevant(paramsStruct)
  if paramsStruct.wasDocked
      try
          set(paramsStruct.hFigure,'WindowStyle','docked');
          %drawnow;
      catch
          % never mind - ignore...
      end
  end
end  % redockFigureIfRelevant
%% Copy screen-capture to the system clipboard
% Adapted from http://www.mathworks.com/matlabcentral/fileexchange/28708-imclipboard/content/imclipboard.m
function imclipboard(imgData)
    % Import necessary Java classes
    import java.awt.Toolkit.*
    import java.awt.image.BufferedImage
    import java.awt.datatransfer.DataFlavor
    % Add the necessary Java class (ImageSelection) to the Java classpath
    if ~exist('ImageSelection', 'class')
        % Obtain the directory of the executable (or of the M-file if not deployed) 
        %javaaddpath(fileparts(which(mfilename)), '-end');
        if isdeployed % Stand-alone mode. 
            [status, result] = system('path');  %#ok<ASGLU>
            MatLabFilePath = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
        else % MATLAB mode. 
            MatLabFilePath = fileparts(mfilename('fullpath')); 
        end 
        javaaddpath(MatLabFilePath, '-end'); 
    end
        
    % Get System Clipboard object (java.awt.Toolkit)
    cb = getDefaultToolkit.getSystemClipboard;  % can't use () in ML6!
    
    % Get image size
    ht = size(imgData, 1);
    wd = size(imgData, 2);
    
    % Convert to Blue-Green-Red format
    imgData = imgData(:, :, [3 2 1]);
    
    % Convert to 3xWxH format
    imgData = permute(imgData, [3, 2, 1]);
    
    % Append Alpha data (not used)
    imgData = cat(1, imgData, 255*ones(1, wd, ht, 'uint8'));
    
    % Create image buffer
    imBuffer = BufferedImage(wd, ht, BufferedImage.TYPE_INT_RGB);
    imBuffer.setRGB(0, 0, wd, ht, typecast(imgData(:), 'int32'), 0, wd);
    
    % Create ImageSelection object
    %    % custom java class
    imSelection = ImageSelection(imBuffer);
    
    % Set clipboard content to the image
    cb.setContents(imSelection, []);
end  %imclipboard
%% Is the provided handle a figure?
function flag = isFigure(hObj)
    flag = isa(handle(hObj),'figure') | isa(hObj,'matlab.ui.Figure');
end  %isFigure
%% Is the provided handle an axes?
function flag = isAxes(hObj)
    flag = isa(handle(hObj),'axes') | isa(hObj,'matlab.graphics.axis.Axes');
end  %isFigure
%% Is the provided handle an image?
function flag = isImage(hObj)
    flag = isa(handle(hObj),'image') | isa(hObj,'matlab.graphics.primitive.Image');
end  %isFigure
%%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%
% find a way in interactive-mode to single-click another Matlab figure for screen-capture
end

end