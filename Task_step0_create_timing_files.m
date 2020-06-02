%Task step0 create timing files. This pulls the task events files from the
%BIDS structure and creates AFNI-style timing files to use in 3dDeconvolve.
%Also creates a regressor file for stimulus size to include as a nuisance
%regressor in the model.
%
% Brock Kirwan
% kirwan@byu.edu
% 6/1/2020

parDir = '/Volumes/Yorick/RSE1_BIDS';
derDir = [parDir '/derivatives'];
rawDir = [parDir '/rawdata'];
% sub-03 has 108 trials in run 1
subjects = {'sub-02' 'sub-03' 'sub-04' 'sub-05' 'sub-06' 'sub-07' 'sub-08' 'sub-09' 'sub-10' 'sub-11' 'sub-12' 'sub-13' 'sub-14' 'sub-15' 'sub-17' 'sub-18' 'sub-19' 'sub-20' 'sub-21' 'sub-22' 'sub-23' 'sub-24'};
trialTypes = {'bus1' 'bus2' 'bus3' 'bus4' 'bus5' 'bus6' 'mal1' 'mal2' 'mal3' 'mal4' 'mal5' 'mal6' 'novel'};

%Not everyone has exactly the same number of TRs per run. This is needed to
%make the nuisance regressor the right length. This can be obtained from
%the func files in BIDS and 3dInfo. This variable is organized by subjects
%(in the same order as the subjects variable above) and scan run.
TRs = [212   228
    189   228
    228   228
    227   226
    228   228
    228   228
    228   228
    228   228
    228   228
    228   108
    228   191
    228   228
    228   228
    228   227
    227   226
    226   228
    228   228
    228   228
    228   228
    226   228
    228   225
    227   227];

%loop through subjects
for sub = 1:length(subjects)
    %read in the event files using the tab-delimited read function
    eventFileRun1 = [rawDir '/' subjects{sub} '/func/' subjects{sub} '_task-rse1_run-1_events.tsv'];
    eventFileRun2 = [rawDir '/' subjects{sub} '/func/' subjects{sub} '_task-rse1_run-2_events.tsv'];
    
    run1 = tdfread(eventFileRun1);
    run2 = tdfread(eventFileRun2);
    
    run1TrialTypesCell = cellstr(run1.trial_type);
    run2TrialTypesCell = cellstr(run2.trial_type);
    
    %create the output directory if it doesn't already exist
    if ~exist([derDir '/' subjects{sub} '/timing_files'], 'dir')
        mkdir([derDir '/' subjects{sub} '/timing_files']);
    end
    
    %loop through trial types, create AFNI onset time file for each
    for tt = 1:length(trialTypes)
        onsets1 = run1.onset(strcmp(run1TrialTypesCell,trialTypes{tt}));
        onsets2 = run2.onset(strcmp(run2TrialTypesCell,trialTypes{tt}));
      
        string1 = [];
        string2 = [];
        %make strings with onsets 
        if isempty(onsets1)
            string1 = '*';
        else
            for i = 1:length(onsets1)
                string1 = [string1 ' ' num2str(onsets1(i))];
            end
        end
        
        if isempty(onsets2)
            string2 = '*';
        else
            for i = 1:length(onsets2)
                string2 = [string2 ' ' num2str(onsets2(i))];
            end
        end
        
        %write out the timing file to the subject's derivatives folder
        fid = fopen([derDir '/' subjects{sub} '/timing_files/' trialTypes{tt} '.txt'], 'w+');
        fprintf(fid, '%s\n', string1);
        fprintf(fid, '%s\n', string2);
        fclose(fid);
    end
    
    %resample the modulator (scaled stimulus area) into TR time steps to
    %use as a nuissance regressor
    mods1 = resample(run1.stim_area_scaled,7,4); %this function wants integer values, so scaling the 3.5s trial and 2s TR by 2.
    mods2 = resample(run2.stim_area_scaled,7,4); 

    vect = [mods1(1:TRs(sub,1)); mods2(1:TRs(sub,2))];
    vect = rescale(vect); %scaling between [0,1]
    
    %write out the nuisance variable
    dlmwrite([derDir '/' subjects{sub} '/timing_files/stimSizeScaled.txt'],vect);
    
end
