%% Script for data preprocessing of the infant audiovisual eyetracking 
% experiment - pilot data collected Dec 21 - Mar 22
% the script was coded by Ana Maria Portugal

% it uses the hdf5 file and the trial information csv file

% the task script was implemented by Lowe Wilsson in Python

clear all

%% add necessary functions
addpath '/Users/k6-c02c12kvlvdl/Dropbox/Scripts/ET_Methods/pupil-size-master/code/helperFunctions'
% we need a function for dynamic offset mean from Kret, & Sjak-Shie, 2019 
% https://github.com/ElioS-S/pupil-size


%% Set general variables

% for interpolation, both pupil and gaze, the maximum gap allowed
interp_maxGap = 150; % in ms, leave empty if no interpolation is needed

% path for stimuli that define the sound condition in the csv trial file
social = 'stimuli/audio/social/';
nonsocial = 'stimuli/audio/nonsocial/';
silent = 'silent';

% Screen limits in the gaze data (distance to screen assumed to be 65 cm)
% Monitor size is 52.69x29.64. Monitor right-center edge has 
% coordinates (22.06, 0), the  center-top edge has coordinates (0, 12.9). 
% The left-bottom corner of the monitor has coordinates (-22.06, -12.9).
width_half = 22.06; % X from -22.06 to 22.06
height_half = 12.9; % Y from -12.9 to 12.9

% Define AOI size - during task presentation the stimuli was within a 
% 9x9 visual degrees size
stim_box = 11; % AOI is increased a bit

% do we want plots and data structs to be saved?
plot_figure = 0; % 1=saves heatmap and a figure with 4 different plots 
save_data = 0; % 1=saves wide and long csv files

% initiate a matrix to store pupil samples across trials and individuals
t_for_B = 1; %this is needed to collect all pupil samples in the baseline 
t_for_R = 1; %this is needed to collect all pupil samples in response 


%% paths for data
% this needs to be edited for other projects
resultsPath = 'C:\Users\giobu365\Documents\pupil_habituation\ET_DATA'
resultsPath_calib = '/Users/k6-c02c12kvlvdl/Documents/UU face pop out exp for Ana M/Calibration/';

cd(resultsPath)
mkdir('../results/') % make directory for plots and results
files = dir('**/*.csv'); % scan folder for all trial info csv files, 
% should be one csv file per participant


for f = 1: length(files)
    
    if f > length(files)
        % this is needed because we'll remove some participants below
        continue
    end
    
    [~,name,ext] = fileparts(files(f).name);
    subject = name(1:3); % experiment subject codes were in the format of
    % "P01" so we segment the first three characters. Might need adaptation
    % in other projects
    
    s_time = name(end-3:end);
    id_name = [subject, '_', s_time] % time was saved because we had duplicates
    
%     % uncomment this to inspect particular subjects
%         if subject == "P17"
%         else
%             continue
%         end
        
    
    %% skip/ remove duplicates - this needs adaptation in case of other projects
    % when specific subjects and times appear we removed their paths from 
    % the files list and reassigned the subject variables
    
    if subject == "p01" && s_time == "1102"
        warning('skip')
        files(f)=[];
        if f > length(files) % if we are at the end of the files list we
            % move on in the loop / end loop
            continue
        else % if not we re assign the subject name and paths to the loop index            
            [~,name,ext] = fileparts(files(f).name);
            subject = name(1:3);
            s_time = name(end-3:end);
            id_name = [subject, '_', s_time]
        end
    elseif subject == "p04" && s_time == "1227"
        warning('skip')
        files(f)=[];
        if f > length(files)
            continue
        else            
            [~,name,ext] = fileparts(files(f).name);
            subject = name(1:3);
            s_time = name(end-3:end);
            id_name = [subject, '_', s_time]
        end
    elseif subject == "p13" && s_time == "1116"
        warning('skip')
        files(f)=[];
        if f > length(files)
            continue
        else            
            [~,name,ext] = fileparts(files(f).name);
            subject = name(1:3);
            s_time = name(end-3:end);
            id_name = [subject, '_', s_time]
        end  
    elseif subject == "P61" && s_time == "1219"
        warning('skip')
        files(f)=[];
        if f > length(files)
            continue
        else            
            [~,name,ext] = fileparts(files(f).name);
            subject = name(1:3);
            s_time = name(end-3:end);
            id_name = [subject, '_', s_time]
        end    
    end
    
    % trial_data_all_wide collects wide format data, one row per subject
    trial_data_all_wide.subject(f,1) = convertCharsToStrings(subject);
    trial_data_all_wide.id_name(f,1) = convertCharsToStrings(id_name);
            
    
    %% get invidiual file paths
    
    % find the hdf5 file - this file has the same filename as the csv file
    % but ends in '_hdf5.hdf5'
    Filename = dir(['**/',name,'_hdf5.hdf5']); 
    if isempty(Filename)
        % if no file was found then move on to the next subject
        disp(warning('No hdf5 file found'))
        continue
    else
        Filename =  [Filename.folder,'/', Filename.name ];
    end
    
    % find the trial information csv file - this just needs to get the
    % pathname and filename already in the files struct
    Filename_csv = [files(f).folder,'/', files(f).name ];
 
    % find the calibration file - this searches based on the subject id 
    % (eg P01) because time was not the same as the eye tracker files
    % it might give us multiple files or files that are not necessarily
    % corresponding to the correct session. This was fixed after manually
    Filename_calib = dir([resultsPath_calib, '**/',subject,'*ease_et_calibration*', '_validation_data.csv' ]);
    
    if isempty(Filename_calib)
        % if no file was found, write NaN in the variables field
        disp(warning('No Calibration file found'))
        trial_data_all_wide.calib_time(f,1) = "No file";
        
        trial_data_all_wide.mean_distance_gaze_to_target(f,1) = NaN;
        trial_data_all_wide.std_distance_gaze_to_target(f,1) = NaN;
        trial_data_all_wide.proportion_of_time_gaze_on_screen(f,1) = NaN;

    elseif length(Filename_calib) > 1
        % if more than on file was found, write NaN in the variables field
        disp(warning('More than one calibration file found'))
        trial_data_all_wide.calib_time(f,1) = "More than one file, fix manually";
        
        trial_data_all_wide.mean_distance_gaze_to_target(f,1) = NaN;
        trial_data_all_wide.std_distance_gaze_to_target(f,1) = NaN;
        trial_data_all_wide.proportion_of_time_gaze_on_screen(f,1) = NaN;
        
    else   
        % if we find only one file, we will get the time on the file so we
        % can check if the corresponds to the dataset (i.e. should be some
        % minutes earlier than the eye tracking session
        trial_data_all_wide.calib_time(f,1) = convertCharsToStrings(Filename_calib.name(end-33:end-20));
        
        Filename_calib =  [Filename_calib.folder,'/', Filename_calib.name ];
        calib_info = ease_2_getcalibinfo(Filename_calib); 
        
        trial_data_all_wide.mean_distance_gaze_to_target(f,1) = calib_info.mean_distance_gaze_to_target;
        trial_data_all_wide.std_distance_gaze_to_target(f,1) = calib_info.std_distance_gaze_to_target;
        trial_data_all_wide.proportion_of_time_gaze_on_screen(f,1) = calib_info.proportion_of_time_gaze_on_screen;        
    end
    
    
    %% load gaze and event buffer from hdf5 file
    gaze_buffer  = h5read(Filename,'/data_collection/events/eyetracker/BinocularEyeSampleEvent');
    event_buffer = h5read(Filename, '/data_collection/events/experiment/MessageEvent');
   
    test  = h5read(Filename,'/data_collection/session_meta_data'); % this only has the date of the session
    experiment_date = string(test.code(3:18)');
    
    %test  = h5read(Filename,'/data_collection/events/experiment/LogEvent'); % this only has system information
    
    % convert char event information to strings for easy reading and access
    for i=1: size(event_buffer.text,2)
        idx= find(isletter(event_buffer.text(:, i)), 1, "last" );
        event_buffer.event(i,1) = convertCharsToStrings(event_buffer.text(1:idx, i)');
    end
        
    
    %% load trial information from the csv file
    % this needs to be adapted to other projects
    % The pilot data collected Dec 21 - Mar 22 included 15 initial 
    % participants (P01-P15) who only view the pilot task and from then all
    % other participants also viewed another task (naturalistic videos)
    % Format of the trial information files is a bit different between
    % these two batches of participants (the indexes of the relevant 
    % columns differ so we used two different functions to import the data)
    
    if str2double(subject(2:3)) >= 15
        % if subject id is higher than 15 
        trial_info = ease_2_popout_gettrialinfo(Filename_csv);
    else
        % if subject is one of the first 15 participants
        trial_info = ease_2_popout_gettrialinfo_old(Filename_csv);
    end
    
    % remove rows that do not have task trial information / block breaks
    toRemove = find(isnan(trial_info.trial_global_start_time))
    trial_info(toRemove,:) = [];


    %% loop across all trials
    
    % choose which trials to plot if plot_figure = 1 above 
    trials_to_plot = randi([1 length(trial_info.trial_loop_thisRepN) ],1,2); % this plots two random trials
    %trials_to_plot = [1, 2, 3]; % this plots specific trials
    %trials_to_plot = [1:length(trial_info.trial_loop_thisRepN)]; % this plots all trials
  
    
    for t= 1:length(trial_info.trial_loop_thisRepN)
        
%         % uncomment this to inspect particular trials
%             if t ~= 45   
%                 continue
%             end
        
        
        % trial_data collects long format data, one row per trial per
        % subject
        % save trial information - id and trial nr
        trial_data.id_name(t,1) = convertCharsToStrings(id_name);
        trial_data.id_date(t,1) = experiment_date;
        trial_data.id_trial(t,1) = t;
        
        %% save trial conditions - type, volume, and order
        
        if strncmp(string(trial_info.audio_filepath(t)), social, length(social))
            trial_data.AudioType(t,1) = 1; % social sound
        elseif strncmp(string(trial_info.audio_filepath(t)), nonsocial, length(nonsocial))
            trial_data.AudioType(t,1) = 2; % nonsocial sound
        elseif strncmp(string(trial_info.audio_filepath(t)), silent, length(silent))
            trial_data.AudioType(t,1) = 0; % silent
        end
        
        trial_data.AudioVolume(t,1) = trial_info.audio_volume(t);
        
        % we get the trial order using the trial number from the trial file
        % - if even is the first repetition, if odd is the second
        % repetition
        trial_data.Order(t,1) = (~mod(trial_info.trial_loop_thisRepN(t),2) == 0) + 1 ;
        
        %% create AOIs for each stimuli category 
        % we use info about position in the screen from trial info file and add the stim_box 
        social_AOI = [trial_info.visual_social_pos_x(t) - stim_box/2, trial_info.visual_social_pos_y(t)- stim_box/2, trial_info.visual_social_pos_x(t) + stim_box/2, trial_info.visual_social_pos_y(t) + stim_box/2];
        geometric_AOI = [trial_info.visual_geometric_pos_x(t) - stim_box/2, trial_info.visual_geometric_pos_y(t)- stim_box/2, trial_info.visual_geometric_pos_x(t) + stim_box/2, trial_info.visual_geometric_pos_y(t) + stim_box/2];
        manmade_AOI = [trial_info.visual_manmade_pos_x(t) - stim_box/2, trial_info.visual_manmade_pos_y(t)- stim_box/2, trial_info.visual_manmade_pos_x(t) + stim_box/2, trial_info.visual_manmade_pos_y(t) + stim_box/2];
        natural_AOI = [trial_info.visual_natural_pos_x(t) - stim_box/2, trial_info.visual_natural_pos_y(t)- stim_box/2, trial_info.visual_natural_pos_x(t) + stim_box/2, trial_info.visual_natural_pos_y(t) + stim_box/2];
        fixation_AOI = [ -stim_box/2, - stim_box/2, stim_box/2, stim_box/2]; % central stim
                       
        
        %% segment trial
        
%         % compare global start time on csv to start of trial
%         trial_global_start_time = trial_info.trial_global_start_time(t);
%         trial_start_ebidx_global = find(event_buffer.event == ['exp1 trial ', num2str(t), ' start' ] );
%         trial_start_ebtime_global = event_buffer.time(trial_start_ebidx_global);
%         trial_global_start_time - trial_start_ebtime_global

        % trial start
        % the start of the trial will be 500 ms before the gaze was
        % captured / the offset of the attention grabber 
        % trial_start_ebtime might preced trial start but we want the 
        % window where gaze was evaluated to be in the fixation AOI
        trial_start_ebidx = find(event_buffer.event == ['exp1 trial ', num2str(t), ' attention grabber end' ] );
        trial_start_event_id = event_buffer.event_id(trial_start_ebidx);
        trial_gaze_captured_ebtime = event_buffer.time(trial_start_ebidx);
        trial_start_ebtime = event_buffer.time(trial_start_ebidx) - 0.500; 

        
        % audio onset
        trial_audio_ebidx = find(event_buffer.event == ['exp1 trial ', num2str(t), ' sound onset' ] );
        trial_audio_event_id = event_buffer.event_id(trial_audio_ebidx);
        trial_audio_ebtime = event_buffer.time(trial_audio_ebidx);
        
        % visual onset
        trial_visual_ebidx = find(event_buffer.event == ['exp1 trial ', num2str(t), ' visual onset' ] );
        trial_visual_event_id = event_buffer.event_id(trial_visual_ebidx);
        trial_visual_ebtime = event_buffer.time(trial_visual_ebidx);
        
        % visual offset = the end of the trial
        trial_end_ebidx = find(event_buffer.event == ['exp1 trial ', num2str(t), ' visual offset' ] );
        trial_end_event_id = event_buffer.event_id(trial_end_ebidx);
        trial_end_ebtime = event_buffer.time(trial_end_ebidx);
        
        % save the time from audio onset to visual onset which is in
        % interstimuli interval (random between 80-400 ms)
        trial_data.isi(t,1) = trial_visual_ebtime - trial_audio_ebtime ;
      
        % Find start and end of trial in indexes on Buffer
        % trial based on 500 ms before gaze captured and visual offset
        %trial_start_gbidx = find(gaze_buffer.event_id >= trial_start_event_id, 1, "first" )
        trial_start_gbidx = find(gaze_buffer.time >= trial_start_ebtime, 1, "first" );
        trial_end_gbidx = find(gaze_buffer.time <= trial_end_ebtime, 1, "last" );
        
        % save total duration of trial - from 500 ms before gaze captured 
        % to visual offset
        trial_data.DurationTrial(t,1) = trial_end_ebtime - trial_start_ebtime;
        
        % segment time buffer
        tb = gaze_buffer.time( trial_start_gbidx:trial_end_gbidx, 1 );
        time_trial = (tb(:) - tb(1))*1000; % convert time buffer to ms from visual onset
        
        % save sampling rate in Hz
        sr = etDetermineSampleRate(tb*1000000); % see function below 
        trial_data.SamplingRate(t,1) = sr;
        
        % segment X buffer
        gaze_lx = gaze_buffer.left_gaze_x( trial_start_gbidx:trial_end_gbidx, 1);
        gaze_rx = gaze_buffer.right_gaze_x( trial_start_gbidx:trial_end_gbidx, 1);
        
        % segment Y buffer
        gaze_ly = gaze_buffer.left_gaze_y( trial_start_gbidx:trial_end_gbidx, 1);
        gaze_ry = gaze_buffer.right_gaze_y( trial_start_gbidx:trial_end_gbidx, 1);  
        
        
        %% Save mean distance to screen
        % Lowe recommended using gaze_buffer.left_eye_cam_y and right to get distance
        distance_ly = gaze_buffer.left_eye_cam_y( trial_start_gbidx:trial_end_gbidx, 1);
        distance_ry = gaze_buffer.right_eye_cam_y( trial_start_gbidx:trial_end_gbidx, 1);
        trial_data.DistanceScreen(t,1) = nanmean( mean( [distance_ly, distance_ry] ,2, 'omitnan') );
      
        
        %% Save missing data before interpolation for the entire trial
        RawMissing = (isnan(gaze_lx) & isnan(gaze_rx)) | (isnan( gaze_ly) & isnan( gaze_ry) );
        trial_data.MissingRawGazeTrial(t,1) = sum(RawMissing) / length(RawMissing) ;
        
        
        %% interpolate and average eyes
        % Data was interpolated linearly over gaps in the data shorter 
        % than 150 ms (as in Kleberg 2019, same as CBCD/Luke Mason's scripts)
        
        [mb, flags] = etInterpBuffer(gaze_lx,gaze_ly, gaze_rx,gaze_ry, tb, interp_maxGap); % see function below
        
        lx_out = mb(:, 1);
        ly_out = mb(:, 2);
        rx_out = mb(:, 3);
        ry_out = mb(:, 4);
        
        gaze_x_trial = mean( [lx_out, rx_out] ,2, 'omitnan');
        gaze_y_trial = mean( [ly_out, ry_out] ,2, 'omitnan');
        
%         % get missing data after interpolation
%         gaze_missing = isnan(gaze_x_trial) | isnan( gaze_y_trial) ;
        

        %% Validate trial
        
        % compute binary vectors for gaze on CS, on screen, outside screen
        
        inAOIfixation =...
            gaze_x_trial >= fixation_AOI(1) &...
            gaze_x_trial <= fixation_AOI(3) &...
            gaze_y_trial >= fixation_AOI(2) &...
            gaze_y_trial <= fixation_AOI(4);
        
        onScreen = ...
            gaze_x_trial <= width_half &...
            gaze_x_trial >= -width_half &...
            gaze_y_trial <= height_half &...
            gaze_y_trial >= -height_half;
        
        outScreen = ...
            gaze_x_trial > width_half |...
            gaze_x_trial < -width_half |...
            gaze_y_trial > height_half |...
            gaze_y_trial < -height_half;
        
        % 1. Trial invalid if interpolated gaze was not at the central
        % attention-grabbing area/AOI for at least 40% of the 500 ms before
        % sound onset        
        sound_start_gbidx = find(tb >= trial_audio_ebtime, 1, "first" );
        check_cs_before = 500; %ms
        check_cs_after = 0;
        crit_minInCS = .4;
        InCS_idx = sound_start_gbidx - (check_cs_before*sr/1000) : sound_start_gbidx + (check_cs_after*sr/1000);
        trial_data.inCSbeforeSound(t,1) = sum( inAOIfixation(InCS_idx) )/ length(inAOIfixation(InCS_idx)) ;
        trial_data.ValidInCSbeforeSound(t,1) = trial_data.inCSbeforeSound(t,1) >= crit_minInCS;
         
        % 2. Trial invalid if interpolated gaze was not at the central
        % attention-grabbing area during the 200 ms before visual onset       
        visual_start_gbidx = find(tb >= trial_visual_ebtime, 1, "first" );
        check_cs_before = 200; %ms
        check_cs_after = 0;
        crit_minInCS = 1;
        InCS_idx = visual_start_gbidx - (check_cs_before*sr/1000) : visual_start_gbidx + (check_cs_after*sr/1000);
        trial_data.inCSbeforeArray(t,1) = sum( inAOIfixation(InCS_idx) )/ length( inAOIfixation(InCS_idx)) ;
        trial_data.ValidInCSbeforeArray(t,1) = trial_data.inCSbeforeArray(t,1) >= crit_minInCS;
        
        % 3.1 Trial invalid if interpolated gaze was missing or outside the screen
        % (ie it was not inside the screen) during the 500 ms after visual
        % onset      
        check_OnScreen_before = 0; %ms
        check_OnScreen_after = 500;
        crit_minOnScreen = 1;
        OnScreen_idx = visual_start_gbidx - (check_OnScreen_before*sr/1000) : visual_start_gbidx + (check_OnScreen_after*sr/1000);
        trial_data.onScreenAfterArray(t,1) = sum( onScreen(OnScreen_idx)) / length( onScreen(OnScreen_idx));
        trial_data.ValidOnScreenAfterArray(t,1) = trial_data.onScreenAfterArray(t,1) >= crit_minOnScreen;
        
        % 3.2 Trial invalid if interpolated gaze was missing or outside the screen
        % (ie it was not inside the screen) for more than 25% during
        % the subsequent 500ms,
        check_OnScreen_after_start = 500; %ms
        check_OnScreen_after_end = 1000;
        crit_minOnScreen = .75;
        OnScreen_idx = visual_start_gbidx + (check_OnScreen_after_start*sr/1000) : visual_start_gbidx + (check_OnScreen_after_end*sr/1000);
        trial_data.onScreenAfterArray500to1000(t,1) = sum( onScreen(OnScreen_idx)) / length( onScreen(OnScreen_idx));
        trial_data.ValidOnScreenAfterArray500to1000(t,1) = trial_data.onScreenAfterArray500to1000(t,1) >= crit_minOnScreen;
        
%        % previous flags/criteria
        %     % 3. Trial invalid if any interpolated gaze was missing during the
        %     % 500 ms after visual onset
        %     check_Missing_before = 0; %ms
        %     check_Missing_after = 500;
        %     crit_maxMissing = 0;
        %     Missing_idx = visual_start_gbidx - (check_Missing_before*sr/1000) : visual_start_gbidx + (check_Missing_after*sr/1000);
        %     trial_data.MissingAfterArray500(t,1) = sum( gaze_missing(Missing_idx)) / length( gaze_missing(Missing_idx));
        %     trial_data.ValidMissingAfterArray500(t,1) = trial_data.MissingAfterArray500(t,1) <= crit_maxMissing;
        %
        %     % 3. Trial invalid if interpolated gaze was missing more than .25
        %     % during the 500-2000 ms after visual onset
        %     check_Missing_after_start = 500; %ms
        %     check_Missing_after_end = 2000;
        %     crit_maxMissing = .25;
        %     Missing_idx = visual_start_gbidx + (check_Missing_after_start*sr/1000) : visual_start_gbidx + (check_Missing_after_end*sr/1000);
        %     trial_data.MissingAfterArray500to2000(t,1) = sum( gaze_missing(Missing_idx)) / length( gaze_missing(Missing_idx));
        %     trial_data.ValidMissingAfterArray500to2000(t,1) = trial_data.MissingAfterArray500to2000(t,1) <= crit_maxMissing;
        %
        %
        %     % 4. Trial invalid if there is interpolated gaze data outside the
        %     % screen in the 2000 ms after visual onset
        %     check_OutScreen_before = 0; %ms
        %     check_OutScreen_after = 1500;
        %     crit_minOutScreen = 0;
        %     OutScreen_idx = visual_start_gbidx - (check_OutScreen_before*sr/1000) : visual_start_gbidx + (check_OutScreen_after*sr/1000);
        %     trial_data.outScreenAfterArray(t,1) = sum( outScreen(OutScreen_idx)) / length( outScreen(OutScreen_idx));
        %     trial_data.ValidOutScreenAfterArray(t,1) = trial_data.outScreenAfterArray(t,1) <= crit_minOutScreen;
              
        
        %% segment visual array presentation
        % from visual onset to visual offset
        
        visual_start_gbidx = find(tb >= trial_visual_ebtime, 1, "first" );
        visual_end_gbidx = find(tb <= trial_end_ebtime, 1, "last" );
        
        idx_visual = visual_start_gbidx:visual_end_gbidx;
        
        gaze_x = gaze_x_trial(idx_visual);
        gaze_y = gaze_y_trial(idx_visual);
        time = tb(idx_visual);
        
        % save total duration of presentation in seconds from number of
        % samples collected and sampling rate (= 3 seconds)
        trial_data.DurationArray(t,1) = size(gaze_x, 1) / sr;
        
%         % get missing data after interpolation during the presentation
%         missing = isnan(gaze_x) | isnan( gaze_y);
                

        %% VISUAL RESPONSE STATS
        
        % compute inAOI binary vectors
        
        inAOIsocial =...
            gaze_x >= social_AOI(1) &...
            gaze_x <= social_AOI(3) &...
            gaze_y >= social_AOI(2) &...
            gaze_y <= social_AOI(4);
        
        inAOIgeometric =...
            gaze_x >= geometric_AOI(1) &...
            gaze_x <= geometric_AOI(3) &...
            gaze_y >= geometric_AOI(2) &...
            gaze_y <= geometric_AOI(4);
        
        inAOImanmade =...
            gaze_x >= manmade_AOI(1) &...
            gaze_x <= manmade_AOI(3) &...
            gaze_y >= manmade_AOI(2) &...
            gaze_y <= manmade_AOI(4);
        
        inAOInatural =...
            gaze_x >= natural_AOI(1) &...
            gaze_x <= natural_AOI(3) &...
            gaze_y >= natural_AOI(2) &...
            gaze_y <= natural_AOI(4);
        
        onScreen = ...
            gaze_x <= width_half &...
            gaze_x >= -width_half &...
            gaze_y <= height_half &...
            gaze_y >= -height_half;
        
        
        % save total time on Screen from number of samples collected and
        % sampling rate (in seconds)
        trial_data.onScreenArray(t,1) = sum( onScreen) / sr;
             
        % get time on each AOI from number of samples collected and
        % sampling rate (in seconds)
        trial_data.onAOI(t,1) = sum(inAOIsocial) / sr;
        trial_data.onAOI(t,2) = sum(inAOIgeometric) / sr;
        trial_data.onAOI(t,3) = sum(inAOImanmade) / sr;
        trial_data.onAOI(t,4) = sum(inAOInatural) / sr;

        
        % Find sample entries to each AOI then
        % find contiguous sections of gaze (>= criteria in ms, critInAOI)
        % in AOI during the visual presentation period to identify a "look"
        
        critInAOI = 100; % minimum time (ms) 
        % this is 24 in CBCD/ Luke Mason's scripts
        % Johan Kleberg used a different criteria based on 4 out of 5
        % samples being in the AOI
        
        
        % SOCIAL AOI
        % find gaps inAOI
        inAOIsocialCont = findcontig(inAOIsocial, true);
        
        if ~isempty(inAOIsocialCont)            
            % save first sample in AOI
            trial_data.inAOIRawLatency(t,1) = time(inAOIsocialCont(1,1)) - time(1);
            
            % convert duration of contiguous sections to ms
            inAOIsocialContMs = inAOIsocialCont(:, 3) * (1000 / sr);
            % find the first section with a duration longer than the criterion
            FirstAOIsocialContFound =...
                find(inAOIsocialContMs >= critInAOI, 1, 'first');
            if ~isempty(FirstAOIsocialContFound)
                % look up the sample of the onset of the first section of gaze that
                % met criterion, convert to time
                trial_data.inAOILatency(t,1) = time(inAOIsocialCont(FirstAOIsocialContFound,1)) - time(1);
            else
                % if no sections met criterion, mark as missing
                trial_data.inAOILatency(t,1) = -1;
            end
        else
            % if no sections were found at all, mark as missing
            trial_data.inAOILatency(t,1) = -1;
            trial_data.inAOIRawLatency(t,1) = -1;
        end
        
        
        % GEOMETRIC AOI
        inAOIgeometricCont = findcontig(inAOIgeometric, true);
        
        if ~isempty(inAOIgeometricCont)
            trial_data.inAOIRawLatency(t,2) = time(inAOIgeometricCont(1,1)) - time(1);
            % convert duration of contiguous sections to ms
            inAOIgeometricContMs = inAOIgeometricCont(:, 3) * (1000 / sr);
            % find the first section with a duration longer than the criterion
            FirstAOIgeometricContFound =...
                find(inAOIgeometricContMs >= critInAOI, 1, 'first');
            if ~isempty(FirstAOIgeometricContFound)
                % look up the sample of the onset of the section of gaze that
                % met criterion, convert to time
                trial_data.inAOILatency(t,2) = time(inAOIgeometricCont(FirstAOIgeometricContFound,1)) - time(1);
            else
                % if no sections met criterion, mark as missing
                trial_data.inAOILatency(t,2) = -1;
            end
        else
            % if no sections were found at all, mark as missing
            trial_data.inAOILatency(t,2) = -1;
            trial_data.inAOIRawLatency(t,2) = -1;
        end
        
        
        % MANMADE AOI
        inAOImanmadeCont = findcontig(inAOImanmade, true);
        
        if ~isempty(inAOImanmadeCont)
            trial_data.inAOIRawLatency(t,3) = time(inAOImanmadeCont(1,1)) - time(1);
            % convert duration of contiguous sections to ms
            inAOImanmadeContMs = inAOImanmadeCont(:, 3) * (1000 / sr);
            % find the first section with a duration longer than the criterion
            FirstAOImanmadeContFound =...
                find(inAOImanmadeContMs >= critInAOI, 1, 'first');
            if ~isempty(FirstAOImanmadeContFound)
                % look up the sample of the onset of the section of gaze that
                % met criterion, convert to time
                trial_data.inAOILatency(t,3) = time(inAOImanmadeCont(FirstAOImanmadeContFound,1)) - time(1);
            else
                % if no sections met criterion, mark as missing
                trial_data.inAOILatency(t,3) = -1;
            end
        else
            % if no sections were found at all, mark as missing
            trial_data.inAOILatency(t,3) = -1;
            trial_data.inAOIRawLatency(t,3) = -1;
        end
        
        
        % NATURAL AOI
        inAOInaturalCont = findcontig(inAOInatural, true);
        
        if ~isempty(inAOInaturalCont)
            trial_data.inAOIRawLatency(t,4) = time(inAOInaturalCont(1,1)) - time(1);
            % convert duration of contiguous sections to ms
            inAOInaturalContMs = inAOInaturalCont(:, 3) * (1000 / sr);
            % find the first section with a duration longer than the criterion
            FirstAOInaturalContFound =...
                find(inAOInaturalContMs >= critInAOI, 1, 'first');
            if ~isempty(FirstAOInaturalContFound)
                % look up the sample of the onset of the section of gaze that
                % met criterion, convert to time
                trial_data.inAOILatency(t,4) = time(inAOInaturalCont(FirstAOInaturalContFound,1)) - time(1);
            else
                % if no sections met criterion, mark as missing
                trial_data.inAOILatency(t,4) = -1;
            end
        else
            % if no sections were found at all, mark as missing
            trial_data.inAOILatency(t,4) = -1;
            trial_data.inAOIRawLatency(t,4) = -1;
        end
        

        % Was there a look in any AOI?
        % if so, save where the first latency happened (which stimuli 
        % category), wwhether it was the face, and save the latency
        if any(trial_data.inAOILatency(t,:) > 0) % if there was a look
            
            % which index is the first look?
            target = find( trial_data.inAOILatency(t,:) > 0 &...
                trial_data.inAOILatency(t,:) == min (trial_data.inAOILatency(t, trial_data.inAOILatency(t,:) > 0)) );
            switch target
                case 1
                    trial_data.FirstLook(t,1) = 'F'; % Face
                    trial_data.FirstLookFace(t,1) = 1; 
                case 2
                    trial_data.FirstLook(t,1) = 'G'; % Geometric
                    trial_data.FirstLookFace(t,1) = 0;
                case 3
                    trial_data.FirstLook(t,1) = 'M'; % Manmade
                    trial_data.FirstLookFace(t,1) = 0;
                case 4
                    trial_data.FirstLook(t,1) = 'N'; % Natural
                    trial_data.FirstLookFace(t,1) = 0;
            end
            
            trial_data.MinLatency(t,:) = min (trial_data.inAOILatency(t, trial_data.inAOILatency(t,:) > 0)) ;
        else
            % if no look was found then mark as missing
            trial_data.FirstLook(t,1) = NaN;
            trial_data.FirstLookFace(t,1) = NaN;
            trial_data.MinLatency(t,1) = NaN;
        end
        
        
        %% continue to validate trial
        
        % 4. Trial invalid if "look" was not detected at any AOI
        % (entry which lasted 100 ms or more) during the entire visual
        % presentation
        trial_data.ValidGazeInAOI(t,1) = any(trial_data.inAOILatency(t,:) > 0);
        
        % 5. Trial invalid if Latency to first AOI is below 200 ms or above 1 s
        trial_data.ValidGazeLatency(t,1) = trial_data.MinLatency(t,1) >= .2 & trial_data.MinLatency(t,1) <= 1;
        
        %% Validate trial based on all flags
        trial_data.ValidGaze(t,1) = trial_data.ValidInCSbeforeSound(t,1) == 1 &&...
            trial_data.ValidInCSbeforeArray(t,1) == 1 &&...
            trial_data.ValidOnScreenAfterArray(t,1) == 1 &&...
            trial_data.ValidOnScreenAfterArray500to1000(t,1) == 1 &&...
            trial_data.ValidGazeInAOI(t,1) == 1 &&...
            trial_data.ValidGazeLatency(t,1) == 1;       
        
        
        %% PUPIL RESPONSE STATS
        
        % segment pupil data
        pupilL = gaze_buffer.left_pupil_measure1( trial_start_gbidx:trial_end_gbidx, 1);
        pupilR = gaze_buffer.right_pupil_measure1( trial_start_gbidx:trial_end_gbidx, 1);
        
        % missing pupil
        pupilL_miss = isnan(pupilL);
        pupilR_miss = isnan(pupilR);
        
        % smooth pupil data with a moving window average (100 ms window)
        % using matlab function movmean(, 'omitnan') but replacing all 
        % missing data with NaN again after        
        mawindow_ms = 100; % 100 ms 
        window = (sr*mawindow_ms) / 1000; % in samples
        pupilL_smooth = movmean(pupilL,window, 'omitnan');
        pupilL_smooth(pupilL_miss) = NaN;
        pupilR_smooth = movmean(pupilR,window, 'omitnan');
        pupilR_smooth(pupilR_miss) = NaN;
        
        % average the pupil signal 
        % if we have data for both eyes
        if sum(~pupilL_miss) > 0 &&  sum(~pupilR_miss) > 0 && sum(~pupilL_miss & ~pupilR_miss) > 0
            
            % average the pupil signal considering a dynamic offset
            % using a function from pupil-size-master (Kret, & Sjak-Shie, 2019)
            pupil_mean = genMeanDiaSamples(time_trial, pupilL_smooth, pupilR_smooth, ~pupilL_miss, ~pupilR_miss);

            % Ana found some errors in this function when the begining of 
            % the signal was missing / only one eye was present so in some
            % cases the normal mean was computed instead
            
            % check which data stream has more data
            streams = ['L', 'R', 'Mean'];
            streams_miss = [sum(pupilL_miss), sum(pupilR_miss), sum(isnan(pupil_mean)) ];
            idx_stream = find( streams_miss == min(streams_miss) );
            
            if any(idx_stream == 3) && ~isempty(pupil_mean)
                trial_data.PupilEye(t,1) = 1; % we used mean pupil with dynamic offset
            else
                pupil_mean =  mean( [pupilL_smooth, pupilR_smooth] ,2, 'omitnan');
                trial_data.PupilEye(t,1) = 2; % we used normal mean
            end
            
        else
            % if one eye is completely missing do normal average which is going
            % to take only one eye.
            pupil_mean =  mean( [pupilL_smooth, pupilR_smooth] ,2, 'omitnan');
            trial_data.PupilEye(t,1) = 3; % we used only one eye
            
            %         % check which data stream as more data
            %         streams = ['L', 'R'];
            %         streams_miss = [sum(pupilL_miss) ,  sum(pupilR_miss)];
            %
            %         idx_stream = find( streams_miss == min(streams_miss) );
            %
            %         if idx_stream == 1
            %             trial_data.Pupil_eye(t,1) = 'Only L     '; % we used left pupil
            %         elseif idx_stream == 2
            %             trial_data.Pupil_eye(t,1) = 'Only R     '; % we used right pupil
            %         end
        end
        
        % save how much missing data before interpolation we had
        trial_data.MissingRawPupil(t,1) = sum(isnan(pupil_mean)) / length(pupil_mean) ;
        
        % Exclude invalid samples based on reasonable mm range (Kret, & Sjak-Shie, 2019)
        pupil_mean(pupil_mean < 1.5) = NaN;
        pupil_mean(pupil_mean > 9) = NaN;
        
        % Exclude invalid samples based on outliers (Mathôt & Vilotijević, 2022)
        Pmean = double( nanmean(pupil_mean));
        P3STD = double( 3*std(pupil_mean, 'omitnan'));
        pupil_mean(pupil_mean > Pmean + P3STD) = NaN;
        pupil_mean(pupil_mean < Pmean - P3STD) = NaN;
        
        % Save how much invalid data we excluded
        trial_data.PupilOutside3STD(t,1) = sum( pupil_mean > Pmean + P3STD | pupil_mean < Pmean - P3STD ) / length(pupil_mean);
        trial_data.PupilOutsideRange(t,1) = sum( pupil_mean < 1.5 | pupil_mean > 9 ) / length(pupil_mean);
        
        % Exclude outside screen samples
        pupil_mean(outScreen) = NaN;
        
        
        % Interpolate missing/invalid pupil data
        % linearly over gaps in the data shorter 
        % than 150 ms (Kleberg 2019, same as CBCD/Luke Mason's scripts)
        pupil_miss = isnan(pupil_mean);
        
        gaps = findcontig(pupil_miss, true);
        
        % check there are missing gaps and also some pupil data
        if size(gaps,1) > 0 && size(time_trial(~pupil_miss), 1) > 2 && ~isempty(interp_maxGap)
            
            pupil_interp = interp1(time_trial(~pupil_miss) ...
                ,pupil_mean(~pupil_miss)...
                ,time_trial,'linear');
            
            % select those runs that are higher than the maximum length that we
            % want to interpolate over
            large_gaps = gaps( gaps(:,3)* (1000 / sr) > interp_maxGap, : );
            
            % loop through all invalid interpolated sections (i.e. missing
            % data was more than criterion) and get where they happened
            indexes_large_gaps = [];
            for i = 1:size(large_gaps,1)
                indexes_large_gaps = horzcat(indexes_large_gaps, linspace(large_gaps(i,1),large_gaps(i,2),large_gaps(i,3)));
            end 
            pupil_interp(indexes_large_gaps) = NaN;  % Set all values in large gaps to NaN again
        else
            
            % if there are no missing gaps just keep the mean pupil
            pupil_interp = pupil_mean;
        end
        
        
        % BASELINE AND RESPONSE STATS
        
        % save mean pupil size in the baseline
        % as well as missing after interpolation
        % baseline is defined as the 200 ms before audio onset
        sound_start_gbidx = find(tb >= trial_audio_ebtime, 1, "first" );
        duration_baseline_before = 200; %ms
        duration_baseline_after = 0;
        baseline_idx = sound_start_gbidx - (duration_baseline_before*sr/1000) : sound_start_gbidx + (duration_baseline_after*sr/1000);
        trial_data.MeanBaseline(t,1) = mean( pupil_interp(baseline_idx)  , 'omitnan');
        trial_data.MissingBaseline(t,1) = sum( isnan( pupil_interp(baseline_idx) )) / size(baseline_idx,2);
      
        % save mean pupil size in the response
        % as well as missing after interpolation
        % response is defined as the 2 seconds after visual onset
        visual_start_gbidx = find(tb >= trial_visual_ebtime, 1, "first" );
        duration_response_before = 0; %ms
        duration_response_after = 2000;
        response_idx = visual_start_gbidx - (duration_response_before*sr/1000) : visual_start_gbidx + (duration_response_after*sr/1000);
        trial_data.MeanResponse(t,1) = mean( pupil_interp(response_idx)  , 'omitnan');
        trial_data.MissingResponse(t,1) = sum( isnan( pupil_interp(response_idx) )) / size(response_idx,2);
        
        % save pupil dilation indexed as mean pupil response - mean pupil baseline
        if ~isnan(trial_data.MeanBaseline(t,1)) && ~isnan(trial_data.MeanResponse(t,1))
            trial_data.PupilDilation(t,1) = trial_data.MeanResponse(t,1) - trial_data.MeanBaseline(t,1);
        else
            % if we are missing baseline or response mark as missing
            trial_data.PupilDilation(t,1) = NaN;
        end
        
        % Collect normalized pupil during baseline and during response for
        % later average tracing / plots
        pupil_normal = (pupil_interp - trial_data.MeanBaseline(t,1)) / trial_data.MeanBaseline(t,1);
        pupil_all.response(:, t_for_B) = pupil_normal(response_idx);
        pupil_all.baseline(:, t_for_R) = pupil_normal(baseline_idx);
        t_for_B = t_for_B +1;
        t_for_R = t_for_R +1;
 
        %% plot trial data
        if plot_figure &&...
                any(t == trials_to_plot)

            % plot gaze data
            spx = subplot(2, 2, 1);
            hold(spx, 'on')
            title('Gaze plot')
            plot(tb, gaze_lx, '-r', 'linewidth', 2)
            plot(tb, gaze_ly, '-g', 'linewidth', 2)
            
            plot(tb, gaze_rx, '-r', 'linewidth', 2)
            plot(tb, gaze_ry, '-g', 'linewidth', 2)
            
            plot(tb, gaze_x_trial, '-k', 'linewidth', 1)
            plot(tb, gaze_y_trial, '-k', 'linewidth', 1)
            
            % plot inCS validation period before sound onset
            sound_start_gbidx = find(tb >= trial_audio_ebtime, 1, "first" );
            check_cs_before = 500; %ms
            check_cs_after = 0;
            pos_sound = [trial_audio_ebtime - (check_cs_before/1000),...
                -width_half,...
                (check_cs_before/1000 + check_cs_after/1000),...
                width_half*2];
            rectangle('Position',pos_sound, 'FaceColor',[0.4940 0.1840 0.5560 0.20])
            
            % plot inCS/ onScreen validation period around visual onset
            check_cs_before = 200; %ms
            check_cs_after = 500;
            pos_visual = [trial_visual_ebtime - (check_cs_before/1000),...
                -width_half,...
                (check_cs_before/1000 + check_cs_after/1000),...
                width_half*2];
            rectangle('Position',pos_visual, 'FaceColor',[0.3010 0.7450 0.9330 0.20])
            
            % plot sound onset
            line([trial_audio_ebtime, trial_audio_ebtime],[-width_half, width_half], 'color', [0, 0, 0],...
                'linewidth', 2);
            % plot visual onset
            line([trial_visual_ebtime, trial_visual_ebtime],[-width_half, width_half], 'color', [0, 0, 0],...
                'linewidth', 2);
            
            legend({'X gaze','Y gaze'},'Location','southwest')
            
            ylim([-width_half, width_half]) % coordinates are X screen ones
            xlim([tb(1), tb(end)])
            
            % another plot: plot inAOI vectors
            spx = subplot(2, 2, 3);
            hold(spx, 'on')
            title('AOI Hit plot')
            
            scatter(time, inAOInatural*5, 100, 'm', 'filled')
            scatter(time, inAOImanmade*4, 100, 'c', 'filled')
            scatter(time, inAOIgeometric*3, 100, 'b', 'filled')
            scatter(time, inAOIsocial*2, 100, 'g', 'filled')
            scatter(tb, inAOIfixation*1, 100, 'r', 'filled')
            
            line([trial_audio_ebtime, trial_audio_ebtime], [0.5,5.5], 'color', [0, 0, 0],...
                'linewidth', 2);
            line([trial_visual_ebtime, trial_visual_ebtime], [0.5,5.5], 'color', [0, 0, 0],...
                'linewidth', 2);
            
            ylim([0.5,5.5])
            xlim([tb(1), tb(end)])
            
            legend({'Natural','Manmade', 'Geometric', 'Face', 'Fixation'},'Location','northwest')

            % another plot: plot pupil tracing
            spx = subplot(2, 2, 2);
            hold(spx, 'on')
            title('Pupil size plot')
            
            plot(tb, pupilL, '-b', 'linewidth', 2)
            plot(tb, pupilR, '-m', 'linewidth', 2)
            plot(tb, pupil_interp, '-k', 'linewidth', 1)
            
            min_pupil = min([pupilL;pupilR; pupil_interp]);
            max_pupil = max([pupilL;pupilR; pupil_interp]);
            
            if ~isnan(min_pupil)
                % pos = [x y w h]
                % plot baseline period
                pos_baseline = [trial_audio_ebtime - (duration_baseline_before/1000),...
                    min_pupil,...
                    (duration_baseline_before/1000 + duration_baseline_after/1000),...
                    max_pupil - min_pupil];
                rectangle('Position',pos_baseline, 'FaceColor',[0.9290 0.6940 0.1250 0.20])
                
                % plot response period
                pos_response = [trial_visual_ebtime - (duration_response_before/1000),...
                    min_pupil,...
                    (duration_response_before/1000 + duration_response_after/1000),...
                    max_pupil - min_pupil];
                rectangle('Position',pos_response, 'FaceColor',[0.8500 0.3250 0.0980 0.20])
                
                line([trial_audio_ebtime, trial_audio_ebtime], [min_pupil-.1, max_pupil+.1], 'color', [0, 0, 0],...
                    'linewidth', 2);
                line([trial_visual_ebtime, trial_visual_ebtime], [min_pupil-.1, max_pupil+.1], 'color', [0, 0, 0],...
                    'linewidth', 2);
                xlim([tb(1), tb(end)])
                ylim([min_pupil-.1, max_pupil+.1])
            end
            
            % another plot: plot some stats
            spx = subplot(2, 2, 4);
            % write down stats
            hold(spx, 'on')
            title('Stats')

            text(0.1, 0.9, strcat('ISI = ', num2str( round(trial_data.isi(t,1)*1000)), 'ms '), 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.9, ['SR =', num2str( trial_data.SamplingRate(t,1) ), 'Hz'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            
            text(0.1, 0.8, ['% In Fixation before Sound onset =', num2str( round( trial_data.inCSbeforeSound(t,1)*100 )), '%'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.8, ['% In Fixation before Visual onset =', num2str( round( trial_data.inCSbeforeArray(t,1)*100 )), '%'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.1, 0.7, ['% On Screen after Visual onset =', num2str( round( trial_data.onScreenAfterArray(t,1)*100 )), '%'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.7, ['% On Screen 500-1000 ms =', num2str( round( trial_data.onScreenAfterArray500to1000(t,1)*100 )), '%'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            
            text(0.1, 0.6, ['Min Latency Raw =', num2str( min( trial_data.inAOIRawLatency(t, trial_data.inAOIRawLatency(t,:) >=0 )) )], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.6, ['Min Latency =', num2str( min( trial_data.inAOILatency(t, trial_data.inAOILatency(t,:) >=0 )) )], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            
            text(0.1, 0.5, ['% Miss Pupil during Baseline =', num2str( round( trial_data.MissingBaseline(t,1)*100 )), '%'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.5, ['% Miss Pupil during Response =', num2str( round( trial_data.MissingResponse(t,1)*100 )), '%'], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            
            
            % mean and STD of pupil size, max and min
            text(0.1, 0.4, ['Pupil Mean =', num2str(mean(pupil_interp, 'omitnan'))], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.4, ['Pupil STD =', num2str(std(pupil_interp, 'omitnan'))], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.1, 0.3, ['Pupil Minimum =', num2str(min(pupil_interp))], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.3, ['Pupil Maximum =', num2str(max(pupil_interp))], 'linewidth', 1, 'edgecolor', [0, 0, 0]);
            
            text(0.1, 0.2, ['Pupil outside 3STD =', num2str( round( trial_data.PupilOutside3STD(t,1) *100 )), '%'],...
                'linewidth', 1, 'edgecolor', [0, 0, 0]);
            text(0.5, 0.2, ['Pupil outside Range (1.5-9) =', num2str( round( trial_data.PupilOutsideRange(t,1)  *100 )), '%'],...
                'linewidth', 1, 'edgecolor', [0, 0, 0]);
            
            text(0.3, 0.1, ['Valid Trial? ', num2str(trial_data.ValidGaze(t,1)), ' (1=Yes, 0=No)'],...
                'linewidth', 1, 'edgecolor', [0, 0, 0]);
            
            
            set(gcf, 'WindowState', 'fullscreen')
            saveas(gcf,['../results/', id_name, ' trial ', num2str(t),' plots.png'])
            %   pause;
            close(gcf)
            
            
            % another figure: plot "heatmap"
            figure, scatter(gaze_x, gaze_y)
            hold on
            rectangle('Position',[-stim_box/2 -stim_box/2 stim_box stim_box])
            rectangle('Position',[12-stim_box/2 6-stim_box/2 stim_box stim_box])
            rectangle('Position',[-12-stim_box/2 6-stim_box/2 stim_box stim_box])
            rectangle('Position',[12-stim_box/2 -6-stim_box/2 stim_box stim_box])
            rectangle('Position',[-12-stim_box/2 -6-stim_box/2 stim_box stim_box])
            % face
            rectangle('Position',[trial_info.visual_social_pos_x(t)-stim_box/2 trial_info.visual_social_pos_y(t)-stim_box/2 stim_box stim_box], 'EdgeColor','r')
            
            xlim([-width_half, width_half])
            ylim([-height_half, height_half])
            
            set(gcf, 'WindowState', 'fullscreen')
            saveas(gcf,['../results/', id_name, ' trial ', num2str(t),' heatmap.png'])
            %pause;
            close(gcf)
            close all            
        end        
    end
    
    % % uncomment to save trial data in an invidiual file
    % if save_data == 1
    %   writetable(struct2table(trial_data), strcat(id_name, '_', experiment_date, '_Dataset_long.csv') )
    % end
    
    %% collect all participants data in a struct
    if ~exist('trial_data_all') % on first subject
        trial_data_all = trial_data;
    else
        % merge trial data for subject with the rest of the trial data
        trial_fields = fieldnames(trial_data);
        for i = 1:length(trial_fields)
            trial_data_all.(trial_fields{i}) = [trial_data_all.(trial_fields{i}); trial_data.(trial_fields{i})];
        end
    end
    
    %% collect average stats for each individual (only valid trials!)
    % see variables saved described below
    
    % select valid trials
    idx_Valid = trial_data.ValidGaze == 1;
    
    % save total number of trials
    trial_data_all_wide.NumberTrials(f,1) = length(trial_data.id_name);
    
    trial_data_all_wide.MissingRawGazeTrial(f,1) = mean(trial_data.MissingRawGazeTrial(idx_Valid));
    trial_data_all_wide.onScreenArray(f,1) = mean(trial_data.onScreenArray(idx_Valid));
    
    outcomes = {'MinLatency';...
        'FirstLookFace';...
        'PupilDilation'};
    
    % select trials for each condition 
    clear idx
    idx.Silent1 = trial_data.AudioType == 0  & trial_data.Order == 1 ; %idx_Silent1
    idx.Silent2 = trial_data.AudioType == 0  & trial_data.Order == 2 ; %idx_Silent2
    
    idx.SocialLow1 = trial_data.AudioType == 1  &...
        trial_data.AudioVolume == .7 & trial_data.Order == 1 ; % idx_SocialLow1
    idx.SocialLow2 = trial_data.AudioType == 1  &...
        trial_data.AudioVolume == .7 & trial_data.Order == 2 ; % idx_SocialLow2
    idx.SocialHigh1 = trial_data.AudioType == 1  &...
        trial_data.AudioVolume == 1 & trial_data.Order == 1 ; % idx_SocialHigh1
    idx.SocialHigh2 = trial_data.AudioType == 1  &...
        trial_data.AudioVolume == 1 & trial_data.Order == 2 ; % idx_SocialHigh2
    
    idx.NonSocialLow1 = trial_data.AudioType == 2  &...
        trial_data.AudioVolume == .7 & trial_data.Order == 1 ;
    idx.NonSocialLow2 = trial_data.AudioType == 2  &...
        trial_data.AudioVolume == .7 & trial_data.Order == 2 ;
    idx.NonSocialHigh1 = trial_data.AudioType == 2  &...
        trial_data.AudioVolume == 1 & trial_data.Order == 1 ;
    idx.NonSocialHigh2 = trial_data.AudioType == 2  &...
        trial_data.AudioVolume == 1 & trial_data.Order == 2 ;
    
    f_idx = fieldnames(idx); % each possible combination of conditions

    % loop through all possible combination of conditions and save number
    % of valid trials
    for c=1:size(f_idx,1) % 10 Ns 
        trial_data_all_wide.(['NumberValidTrials', f_idx{c}])(f,1) = sum(idx.(f_idx{c}) & idx_Valid);
    end
    
    % loop through all measures and possible combination of conditions and 
    % save mean across valid trials
    for o = 1:size(outcomes,1) % 3 outcomes
        for c=1:size(f_idx,1) % 10 conditions = 30 means
            if trial_data_all_wide.(['NumberValidTrials', f_idx{c}])(f,1) > 0
                trial_data_all_wide.(['Mean',outcomes{o}, f_idx{c}])(f,1) = nanmean( trial_data.(outcomes{o})(idx.(f_idx{c}) &...
                    idx_Valid) );
            else
                % if there are no valid trials mark as missing
                trial_data_all_wide.(['Mean',outcomes{o}, f_idx{c}])(f,1) = NaN;
            end
        end
    end

    % loop through all measures except FirstLookFace and possible 
    % combination of conditions and save std across valid trials
    for o = 1:size(outcomes,1) % 3 outcomes
        if outcomes{o} == "FirstLookFace" 
            % we do not compute std for FirstLookFace because this is a
            % binary variable per trial
            continue
        end
        for c=1:size(f_idx,1) % 10 conditions = 20 stds
            if trial_data_all_wide.(['NumberValidTrials', f_idx{c}])(f,1) > 0
                trial_data_all_wide.(['STD',outcomes{o}, f_idx{c}])(f,1) = nanstd( trial_data.(outcomes{o})(idx.(f_idx{c}) &...
                    idx_Valid) );
            else
                % if there are no valid trials mark as missing
                trial_data_all_wide.(['STD',outcomes{o}, f_idx{c}])(f,1) = NaN;
            end
        end
    end
    
    % clear vars
    clear trial_data
    clearvars -except trial_data_all trial_data_all_wide...
        id interp_maxGap social nonsocial silent width_half height_half plot_figure save_data stim_box...
        resultsPath resultsPath_calib files...
        pupil_all t_for_B t_for_R
    
end

if save_data == 1
    
    %% write long file
    writetable(struct2table(trial_data_all), ['../results/', date, ' All_Datasets_long.csv'] )
    
    % % header are:
    % id_name	Participant ID
    % id_date	Experiment date
    % id_trial	Trial nr
    % AudioType	0 = silent, 1 = social sound, 2 = nonsocial sound
    % AudioVolume	0 = silent, 0.7 = low volume, 1 = high volume
    % Order	1 = first repetition, 2 = second repetition
    % isi	Time (seconds) between sound onset and visual array onset (80-400 ms)
    % DurationTrial	Time (seconds) from 500 ms before sound onset was triggered (gaze was captured) to visual array offset
    % SamplingRate	Number of samples per second, 600 Hz
    % MissingRawGazeTrial	Proportion raw (before interpolation) missing gaze data in the trial (from 500 ms before sound onset to visual array offset)
    % inCSbeforeSound	Trial validity flag 1: proportion of interpolated gaze that was at the central fixation point during the 500 ms before audio onset
    % ValidInCSbeforeSound	Trial validity flag 1: 0 = invalid = inCSbeforeSound < .4, 1 = valid
    % inCSbeforeArray	Trial validity flag 2: proportion of interpolated gaze that was at the central fixation point during the 200 ms before visual onset
    % ValidInCSbeforeArray	Trial validity flag 2: 0 = invalid = inCSbeforeArray < 1, 1 = valid
    % onScreenAfterArray	Trial validity flag 3: proportion of interpolated gaze that was inside the screen during the 500 ms after visual onset
    % ValidOnScreenAfterArray	Trial validity flag 3: 0 = invalid = onScreenAfterArray < 1, 1 = valid
    % onScreenAfterArray500to1000	Trial validity flag 3: proportion of interpolated gaze that was inside the screen during the 500-1000 ms after visual onset
    % ValidOnScreenAfterArray500to1000	Trial validity flag 3: 0 = invalid = onScreenAfterArray < .75, 1 = valid
    % DurationArray	Time (seconds) from visual onset to visual array offset
    % onScreenArray	Cumulative time (seconds) in the screen during visual presentation (visual onset to visual offset)
    % onAOI_1	Cumulative time (seconds) in each AOI (1 is Face, 2 is Geometric, 3 is Manmade, 4 is Natural)
    % onAOI_2	Cumulative time (seconds) in each AOI (1 is Face, 2 is Geometric, 3 is Manmade, 4 is Natural)
    % onAOI_3	Cumulative time (seconds) in each AOI (1 is Face, 2 is Geometric, 3 is Manmade, 4 is Natural)
    % onAOI_4	Cumulative time (seconds) in each AOI (1 is Face, 2 is Geometric, 3 is Manmade, 4 is Natural)
    % inAOIRawLatency_1	Latency (seconds) of first sample detected in AOI (1 is Face, 2 is Geometric, 3 is Manmade, 4 is Natural). -1 means gaze sample not detected.
    % inAOIRawLatency_2	Latency (seconds) of first sample detected in AOI (1 is Face, 2 is Geometric, 3 is Manmade, 4 is Natural). -1 means gaze sample not detected.
    % inAOIRawLatency_3	Latency (seconds) of first sample detected in AOI (1 is Face, 2 is Geometric, 3 is Manmade, 4 is Natural). -1 means gaze sample not detected.
    % inAOIRawLatency_4	Latency (seconds) of first sample detected in AOI (1 is Face, 2 is Geometric, 3 is Manmade, 4 is Natural). -1 means gaze sample not detected.
    % inAOILatency_1	Latency (seconds) of first look detected in AOI (1 is Face, 2 is Geometric, 3 is Manmade, 4 is Natural). -1 means look not detected. Look is an entry which lasted 100 ms or more.
    % inAOILatency_2	Latency (seconds) of first look detected in AOI (1 is Face, 2 is Geometric, 3 is Manmade, 4 is Natural). -1 means look not detected. Look is an entry which lasted 100 ms or more.
    % inAOILatency_3	Latency (seconds) of first look detected in AOI (1 is Face, 2 is Geometric, 3 is Manmade, 4 is Natural). -1 means look not detected. Look is an entry which lasted 100 ms or more.
    % inAOILatency_4	Latency (seconds) of first look detected in AOI (1 is Face, 2 is Geometric, 3 is Manmade, 4 is Natural). -1 means look not detected. Look is an entry which lasted 100 ms or more.
    % FirstLook	Target of first look. Empty if no look was detected.
    % FirstLookFace	Whether the target of first look was the face (1 = TRUE, 0 = FALSE). NaN if no look was detected.
    % MinLatency	Latency (seconds) of first look detected.
    % ValidGazeInAOI	Trial validity flag 4: 0 = invalid = no look detected, 1 = valid
    % ValidGazeLatency	Trial validity flag 5: 0 = invalid = latency of first look detected is not > .2 and < 1 seconds, 1 = valid
    % ValidGaze	Whether the trial is valid: 0 = invalid = any of the above flags are 0, 1 = valid
    % PupilEye	1 = mean pupil with dynamic offset, 2 = normal mean pupil, 2 = one eye pupil
    % MissingRawPupil	Proportion raw (before interpolation) missing pupil data in the trial (from 500 ms before sound onset to visual array offset)
    % PupilOutside3STD	Proportion of pupil samples that were excluded because they were outliers (above or below 3STD of trial mean)
    % PupilOutsideRange	Proportion of pupil samples that were excluded because they were invalid (below or above 1.5-9 mm)
    % MeanBaseline	Mean pupil size (mm) during baseline period
    % MissingBaseline	Proportion missing data after interpolation during baseline
    % MeanResponse	Mean pupil size during (mm) response period
    % MissingResponse	Proportion missing data after interpolation during response
    % PupilDilation	Mean Baseline from Mean Response (mm)
    
    %% write wide file
    writetable(struct2table(trial_data_all_wide), ['../results/', date, ' All_Datasets_wide.csv'] )
    
    % % headers are:
    % subject	Participant ID
    % id_name	Participant ID with time
    % calib_time	Calibration date and time
    % mean_distance_gaze_to_target	Mean distance (pixels) from gaze to target during calibration validation
    % std_distance_gaze_to_target	STD distance (pixels) from gaze to target during calibration validation
    % proportion_of_time_gaze_on_screen	Proportion gaze on the screen during calibration validation
    % NumberTrials	Number of trials presented
    % MissingRawGazeTrial	Mean proportion raw (before interpolation) missing gaze data across valid trials (from 500 ms before sound onset to visual array offset)
    % onScreenArray	Mean cumulative time (seconds) in the screen during visual presentation (visual onset to visual offset) across valid trials
    % NumberValidTrialsSilent1	Number of valid trials in condition
    % NumberValidTrialsSilent2	Number of valid trials in condition
    % NumberValidTrialsSocialLow1	Number of valid trials in condition
    % NumberValidTrialsSocialLow2	Number of valid trials in condition
    % NumberValidTrialsSocialHigh1	Number of valid trials in condition
    % NumberValidTrialsSocialHigh2	Number of valid trials in condition
    % NumberValidTrialsNonSocialLow1	Number of valid trials in condition
    % NumberValidTrialsNonSocialLow2	Number of valid trials in condition
    % NumberValidTrialsNonSocialHigh1	Number of valid trials in condition
    % NumberValidTrialsNonSocialHigh2	Number of valid trials in condition
    % MeanMinLatencySilent1	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencySilent2	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencySocialLow1	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencySocialLow2	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencySocialHigh1	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencySocialHigh2	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencyNonSocialLow1	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencyNonSocialLow2	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencyNonSocialHigh1	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanMinLatencyNonSocialHigh2	Mean latency (seconds) of first look detected across valid trials in condition
    % MeanFirstLookFaceSilent1	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceSilent2	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceSocialLow1	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceSocialLow2	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceSocialHigh1	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceSocialHigh2	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceNonSocialLow1	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceNonSocialLow2	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceNonSocialHigh1	Mean proportion of valid trials in condition where the first look was on the face
    % MeanFirstLookFaceNonSocialHigh2	Mean proportion of valid trials in condition where the first look was on the face
    % MeanPupilDilationSilent1	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationSilent2	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationSocialLow1	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationSocialLow2	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationSocialHigh1	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationSocialHigh2	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationNonSocialLow1	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationNonSocialLow2	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationNonSocialHigh1	Mean pupil dilation (response-baseline) across valid trials in condition
    % MeanPupilDilationNonSocialHigh2	Mean pupil dilation (response-baseline) across valid trials in condition
    % STDMinLatencySilent1	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencySilent2	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencySocialLow1	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencySocialLow2	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencySocialHigh1	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencySocialHigh2	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencyNonSocialLow1	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencyNonSocialLow2	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencyNonSocialHigh1	STD latency (seconds) of first look detected across valid trials in condition
    % STDMinLatencyNonSocialHigh2	STD latency (seconds) of first look detected across valid trials in condition
    % STDPupilDilationSilent1	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationSilent2	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationSocialLow1	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationSocialLow2	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationSocialHigh1	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationSocialHigh2	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationNonSocialLow1	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationNonSocialLow2	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationNonSocialHigh1	STD pupil dilation (response-baseline) across valid trials in condition
    % STDPupilDilationNonSocialHigh2	STD pupil dilation (response-baseline) across valid trials in condition

end

%% plot average pupil trace

if 1
    
    pupil_baseline = mean( pupil_all.baseline ,2, 'omitnan');
    tb_baseline = linspace(-280,-80,size(pupil_baseline, 1)); % edit if baseline timing changes
    
    pupil_response = mean( pupil_all.response ,2, 'omitnan');
    tb_response = linspace(0,2000,size(pupil_response, 1)); % edit if response timing changes
    
    min_pupil = min([pupil_baseline;pupil_response]);
    max_pupil = max([pupil_baseline;pupil_response]);
    
    % this needs to be adapted depending on the conditions that are
    % relevant
    % Elin/ Matilda requested social vs non-social and repetition effect        
    idx_social1 = trial_data_all.ValidGaze == 1  & trial_data_all.AudioType == 1  & trial_data_all.Order == 1 ; 
    idx_social2 = trial_data_all.ValidGaze == 1  & trial_data_all.AudioType == 1  & trial_data_all.Order == 2 ; 
    idx_nonsocial1 = trial_data_all.ValidGaze == 1  & trial_data_all.AudioType == 2  & trial_data_all.Order == 1 ; 
    idx_nonsocial2 = trial_data_all.ValidGaze == 1  &trial_data_all.AudioType == 2  & trial_data_all.Order == 2 ;
 
    
    idx = [idx_social1, idx_social2, idx_nonsocial1, idx_nonsocial2];
    idx_color = ['-r'; ':r'; '-b'; ':b'];

%% loop through conditions defined in idx 
for c = 1:4
    pupil_baseline_mean = mean( pupil_all.baseline(:, idx(:,c)) ,2, 'omitnan');
    pupil_baseline_std =  std( pupil_all.baseline(:, idx(:,c)),0, 2 , 'omitnan') / sqrt(length(pupil_all.baseline(:, idx(:,c))));
    
    % plot the average trace for the specific condition
    plot(tb_baseline, pupil_baseline_mean,idx_color(c,:), 'linewidth', 2)
    hold on 
end

for c = 1:4
    pupil_response_mean = mean( pupil_all.response(:, idx(:,c)) ,2, 'omitnan');
    pupil_response_std =  std( pupil_all.response(:, idx(:,c)),0, 2 , 'omitnan') / sqrt(length(pupil_all.response(:, idx(:,c))));
    
    plot(tb_response, pupil_response_mean, idx_color(c,:), 'linewidth', 2)
    hold on
end

    line([-80, -80], [-.1, 0.1], 'color', [0, 0, 0],...
        'linewidth', 2);
    line([0, 0], [-.1, 0.1], 'color', [0, 0, 0],...
        'linewidth', 2);
    
    pos_baseline = [-280,-.1, 200,.2];
    rectangle('Position',pos_baseline, 'FaceColor',[0.8 0.8 0.8 0.25])
    
    legend('Social 1', 'Social 2', 'Non-social 1', 'Non-social 2')
    ylim([-.02, 0.08])
    xlabel('Time (ms)')
    ylabel('Pupil dilation magnitude')
    ax = gca; 
    ax.FontSize = 16; 
    
end

%% this script was first drafted by using data from 4 subjects

% id_name="P29"

% id_name="P32"

% id_name="P49"

% id_name="P50"


%% Up-sample to 1000Hz (not implemented)
%     % Copied from pupil-size-master code (Kret, & Sjak-Shie, 2019)
%     interp_upsamplingFreq = 1000;
%     % Generate the upsampled time vector (seconds):
%     t_upsampled = ...
%         (time(1)/1000 ...
%         :(1/interp_upsamplingFreq)...
%         :time(end)/1000)';
%     diaInterp = interp1(time./1000 ...
%         ,pupilL...
%         ,t_upsampled,'linear'); % this does not interpolate missing data
%     % the resulting vector diaInterp has a frequency of 1000 Hz but
%     % has missing data

%% Filter the data (not implemented)

%    % Low pass filter
%    % Copied from pupil-size-master code (Kret, & Sjak-Shie, 2019)
%    % Calculate the low pass filter specs using the cutoff
%    % frequency [Hz], filter order, and the upsample frequency
%    % specified above:
%    LpFilt_cutoffFreq         = 4;
%    LpFilt_order              = 4;
%    [LpFilt_B,LpFilt_A] ...
%        = butter(LpFilt_order,2*LpFilt_cutoffFreq...
%        /interp_upsamplingFreq );
%    diaInterp_filtered = filtfilt(LpFilt_B...
%        ,LpFilt_A...
%        ,diaInterp);
%     % does not work: error in filtfilt function
%
%     % Multi-step filter by (Kret, & Sjak-Shie, 2019)
%     [valOut,speedFiltData,devFiltData] = rawDataFilter(time, pupilL); #
%     % function from pupil-size-master (Kret, & Sjak-Shie, 2019)
%     % does not work: error in filtfilt function

%% functions needed

function [sRate, msPerS] = etDetermineSampleRate(timeBuffer)
%AMP copied this function from CBCD Luke Mason's database (used in TABLET,
% Face pop-out, probably other Basis tasks)

if isempty(timeBuffer)
    error('Time buffer is empty.')
end

totTimeMs = double(timeBuffer(end, 1) - timeBuffer(1, 1)) / 1000000;
sRate = floor(1 / (totTimeMs / size(timeBuffer, 1)));
msPerS = (totTimeMs / size(timeBuffer, 1)) * 1000;

end

function [dataOut, flags] = etInterpBuffer(lx, ly, rx, ry, timeBuffer, maxMs, dontFilter)

% AMP copied this function from CBCD Luke Mason's database (used in TABLET,
% Face pop-out, probably other Basis tasks)

if ~exist('dontFilter', 'var') || isempty(dontFilter)
    dontFilter = false;
end


flags = false(size(timeBuffer, 1), 4);

% mb = mainBuffer; clear mainBuffer;
tb = timeBuffer; clear timeBuffer;

% get rid of invalid/offscreen data
if ~dontFilter
    %mb = etFilterGazeOnscreen(mb);
end

% get timing data on buffer
[~, msPerS] = etDetermineSampleRate(tb*1000000);
maxSamp = round(maxMs / msPerS);

% % get gaze data
% lx = mb(:, 7);
% ly = mb(:, 8);
% rx = mb(:, 20);
% ry = mb(:, 21);

% find nans
lx_nan = isnan(lx);
ly_nan = isnan(ly);
rx_nan = isnan(rx);
ry_nan = isnan(ry);

% interpolate if there is more than one valid and one invalid sample
t = 1:numel(lx);
if sum(isnan(lx)) > 0 && sum(~isnan(lx)) > 1
    lx_int = interp1(t(~lx_nan), lx(~lx_nan), t, 'linear');
else
    lx_int = lx;
end

if sum(isnan(ly)) > 0 && sum(~isnan(ly)) > 1
    ly_int = interp1(t(~ly_nan), ly(~ly_nan), t, 'linear');
else
    ly_int = ly;
end

if sum(isnan(rx)) > 0 && sum(~isnan(rx)) > 1
    rx_int = interp1(t(~rx_nan), rx(~rx_nan), t, 'linear');
else
    rx_int = rx;
end

if sum(isnan(ry)) > 0 && sum(~isnan(ry)) > 1
    ry_int = interp1(t(~ry_nan), ry(~ry_nan), t, 'linear');
else
    ry_int = ry;
end

% lx_int = interp1q(t(~lx_nan), lx(~lx_nan), t);
% ly_int = interp1q(t(~rx_nan), lx(~rx_nan), t);
% rx_int = interp1q(t(~ly_nan), lx(~ly_nan), t);
% ry_int = interp1q(t(~ry_nan), lx(~ry_nan), t);

% replace interpolated data where the number of missing samples is
% greater than the maximum specified
lx_idx = etInterp_makeIdx(lx, lx_nan, maxSamp);
lx_out = lx;
lx_out(lx_idx) = lx_int(lx_idx);

ly_idx = etInterp_makeIdx(ly, ly_nan, maxSamp);
ly_out = ly;
ly_out(ly_idx) = ly_int(ly_idx);

rx_idx = etInterp_makeIdx(rx, rx_nan, maxSamp);
rx_out = rx;
rx_out(rx_idx) = rx_int(rx_idx);

ry_idx = etInterp_makeIdx(ry, ry_nan, maxSamp);
ry_out = ry;
ry_out(ry_idx) = ry_int(ry_idx);

% store samples that were interpolated in flags output var
flags = any([lx_idx; ly_idx; rx_idx; ry_idx], 1);

% store in buffer
%dataOut = mb;
dataOut = [lx_out, ly_out, rx_out, ry_out];

end

function idx = etInterp_makeIdx(x, xnan, maxSamp)
idx = false(1, length(x));

% find contigous runs of nans, measure length of each
tmp = findcontig(xnan, 1);

% if none found, return
if isempty(tmp)
    return
else
    % select those runs that are shorter than the maximum length that
    % we'll interpolate for (default is 150ms)
    tmp = tmp(tmp(:, 3) <= maxSamp, :);
    
    % if none found, return empty
    if isempty(tmp), idx = []; end
    
    % otherwise loop through and mark as valid those segments that are
    % short enough to interpolate
    for e = 1:size(tmp, 1)
        idx(tmp(e, 1):tmp(e, 2)) = true;
    end
end
end
