%**************************************************************************
% Goal: Generate metadata table for first HHMM manuscript
%**************************************************************************

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
% info about tag ID strings and tag types
%*************************************************************************
whales = {
    "Zica-20191012-144029";
    "Zica-20191012-145101";
    "Zica-20191111-94810";
    "Zica-20191117-195993";
    "Zica-20211112-94819";
    "Zica-20211113-195993";
    "Zica-20220112-195994";
    "Zica-20230518-233391";
    "Zica-20230519-232950";
    "Zica-20230723-233394";
    "Zica-20230723-233395";
    "Zica-20240227-233396";
    "Zica-20240227-240128"
} ;

record_start_utc = strings(size(whales));
tagon_lat = NaN(size(whales));
tagon_lon = tagon_lat;
sensor_record_dur_hr = NaN(size(whales));
data_cropped = zeros(size(whales));
afs = zeros(size(whales));
acc_fs = afs;
mag_fs = afs;
z_fs = afs;
temp_fs = afs;

% loop over whales and process coarse scale data
%*************************************************************************
for w = 1:length(whales)
    disp(append('Processing whale:', ' ', whales{w}))
    
    % Read in NC file
    %*************************************************************************
    tag = whales{w};
    load_nc(strcat(smrt_path, tag, '-cal'));

    afs(w) = SA.sampling_rate;
    acc_fs(w) = Aw.sampling_rate;
    mag_fs(w) = Mw.sampling_rate;
    z_fs(w) = depth.sampling_rate;
    temp_fs(w) = internal_temp.sampling_rate;

    
    if exist('depth_corrected', 'var')
        depth = depth_corrected;
    end
    if exist('CorrectedDepth', 'var')
        depth = CorrectedDepth;
    end
    
    % Tagon time information (for converting with UTC times)
    % tag start time UTC as matlab datenum
    stimenum = datenum(info.dephist_device_datetime_start, info.dephist_device_regset);
    % tag start time UTC as string
    record_start_utc(w) = string(datetime(stimenum, "ConvertFrom", "datenum"));

    % tagon location
    tagon_lat(w) = info.dephist_deploy_location_lat;
    tagon_lon(w) = info.dephist_deploy_location_lon;

    % record duration
    % most are cropped to tagoff.
    % in this dataset, all that were not cropped have been checked and
    % it's verified that the tag is on animal for full record duration
    if isfield(depth, "crop")
        sensor_record_dur_hr(w) = diff(depth.crop) / 60 / 60;
        data_cropped(w) = 1;
    else
        sensor_record_dur_hr(w) = size(depth.data, 1) * depth.sampling_rate / 60 / 60;
    end

    clearvars -except smrt_path whales w record_start_utc ...
        tagon_lat tagon_lon sensor_record_dur_hr data_cropped...
        afs acc_fs mag_fs z_fs temp_fs
end % loop over whales

whale_tab = table(whales, record_start_utc, tagon_lat, tagon_lon, ...
    sensor_record_dur_hr, data_cropped, afs, acc_fs, mag_fs, z_fs);


% Write output to file
%######################################################################
writetable(whale_tab, [smrt_path, 'smrt_hhmm_metadata.csv']);