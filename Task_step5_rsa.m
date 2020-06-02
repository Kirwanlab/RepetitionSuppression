%Task_step5_rsa.m
%
% Script to use the functionally-defined seed region in the right precuneus
% as a seed region in a representational similarity analysis (RSA)
%
% Brock Kirwan
% kirwan@byu.edu
% 6/1/2020

parDir = '/Volumes/Yorick/RSE1_BIDS';
derDir = [parDir '/derivatives'];
outDir = [derDir '/grp-analysis-RSA'];
tempDir = [parDir '/code/templates'];
afniPath = '/usr/local/bin/afni'; %where you have afni installed

prefs.subjects = {...
'sub-02'
'sub-03'
'sub-04'
'sub-05'
'sub-06'
'sub-07'
'sub-08'
'sub-09'
'sub-10'
'sub-11'
'sub-12'
'sub-13'
'sub-14'
'sub-15'
'sub-17'
'sub-18'
'sub-19'
'sub-20'
'sub-21'
'sub-22'
'sub-23'
'sub-24'
};


% %Functionally-defined ROI mask based on main effect of repetition,
% thresholded at p<.001 and spatial extent threshold k=12 (nn=2). This
% selects the precuneus cluster as a seed region for the RSA

prefs.title = 'ME Repetition';
prefs.filename = 'MVM4_p.001_k12_nn2_fixed_mask+tlrc';
prefs.validROInum = [15]; 
% ROIs: 
prefs.validROIname = {'RPrecuneus'};

%deconv with regressor for novel
prefs.deconv = 'rse1_4_blur5+tlrc';
prefs.conds = [2:13]; %these are not 0-indexed.
prefs.condNames = {
    'bus1'
    'bus2'
    'bus3'
    'bus4'
    'bus5'
    'bus6'
    'mal1'
    'mal2'
    'mal3'
    'mal4'
    'mal5'
    'mal6'};

% command = ['rm ' prefs.filename '.txt'];
% unix(command);


%Run 3dROIstats
meanData = zeros(length(prefs.conds),length(prefs.validROInum));
outName = [outDir '/' prefs.filename '.txt'];

%generate the 3dROIstats output if haven't already
if ~exist(outName,'file')
    %put all the subjects into a string for the loop
    subjects = [];
    for i = 1:length(prefs.subjects)
        subjects = [subjects derDir '/' prefs.subjects{i} ...
            '/' prefs.deconv ' '];
    end
    
    %you'll have to change the relative path for AFNI
    %[status, path2program] = unix('which 3dROIstats');
    path2program = [afniPath '/3dROIstats'];
     
    cmd =['for i in ' subjects '; do ' path2program ' -mask '...
        prefs.filename ' -mask_f2short ${i} >> '...
        outName '; done'];
    unix(cmd);
end

fid = fopen(outName);
temp = textscan(fid, '%s', 'delimiter', '\t'); %read whole thing into temp variable
fclose(fid);

%count the number of ROIs
nroi = 0;
for i = 3:length(temp{:}) %assumes first 2 elements are 'File' and 'Sub-brik'
    if strncmpi('Mean',temp{1}(i),4) == 1
        nroi = nroi + 1;
    else
        break
    end
end

%set up the string format according to how many ROIs you have
form = ['%s %s'];
for i = 1:nroi
    form = [form ' %s'];
end

%open the file and read in the data again in correct format
%this gives you an array of cells with #ROIs + 2 columns
fid = fopen(outName);
file = textscan(fid, form, 'delimiter', '\t');
fclose(fid);

%now count how many conditions you have
ncond = 0;
for i = 2:length(file{1}) %assumes the first one is 'File'
    if strcmp('File',file{1}(i)) == 0
        ncond = ncond + 1;
    else
        break
    end
end

%count the number of subjects
nsub = 0;
for i = 1:length(file{1})
    if strcmp('File',file{1}(i)) == 1
        nsub = nsub + 1;
    end
end

%put the data into a matrix
data = zeros(ncond,nroi,nsub);
for cond = 1:ncond
    sub = 1;
    for i = 1:length(file{2})
        condName = char(file{2}(i)); %3dROIstats compile date May 27 2008 outputs more info than it used to...
        numEnd = strfind(condName,'[')-1;
        if strcmp(condName(1:numEnd),num2str(cond-1)) == 1
            for roi = 1:nroi
                data(cond,roi,sub) = str2num(file{roi+2}{i});
            end
            sub = sub+ 1;
        end
    end
end

%drop the conditions you don't care about, keep the ones you do
data = data(prefs.conds,:,:);
ncond = length(prefs.conds);

%and drop the ROI's you're not interested in
data = data(:,prefs.validROInum,:);
nroi = length(prefs.validROInum);

%Correlation analysis

%load the cortical mask using the BrikLoad function from here: https://sscc.nimh.nih.gov/afni/matlab/
[mask_err, mask, mask_info] = BrikLoad([tempDir '/Intersection_GM_mask+tlrc.BRIK']);

%loop through subjects
 for sub = 1:nsub

    %tell me who you're working on
    display(['working on subject ' prefs.subjects{sub}]);
    
    %load deconvolve data
    [err, V, Info, ErrMessage] = BrikLoad([derDir '/' prefs.subjects{sub} '/' prefs.deconv]);

    %loop through the voxels, calculate correlation with mean betas from
    %RSUB ROI (also loop through ROIs?)
    dims = size(V);
    r=zeros(dims(1:3));
    for x = 1:dims(1)
        for y = 1:dims(2)
            for z = 1:dims(3)
                if mask(x,y,z) == 1 %if it's within the cortical mask
                    r(x,y,z) = corr(squeeze(V(x,y,z,prefs.conds)),data(:,1,sub));
                end
            end
        end
    end
   
%     %write out the correlation file for this subject. Can skip since
%     %we'll do all stats on the nomalized datasets below
%     opt.Prefix = [prefs.subjects{sub} '_r'];
%     opt.View = '+tlrc';
%     Info.BRICK_LABS='Correl';
%     if ~isfile([opt.Prefix opt.View '.BRIK'])
%         WriteBrik(r,Info,opt);
%     end

    %z-transform the correlation
    z = atanh(r);
    
    %get rid of NANs
    z(isnan(z)) = 0;
    
    %write out the correlation file for this subject
    opt.Prefix = [outDir '/' prefs.subjects{sub} '_z'];
    opt.View = '+tlrc';
    Info.BRICK_LABS='zCorrel';
    if ~isfile([opt.Prefix opt.View '.BRIK'])
        WriteBrik(r,Info,opt);
    end
 end


 %do the ttest. Again, update for where your afni is installed.
 cmd = [afniDir '/3dttest++ -prefix ' outDir '/rprecun_correl_vs_0 -mask ' tempDir '/Intersection_GM_mask+tlrc -setA '];
 for sub = 1:nsub
     cmd = [cmd ' ' outDir '/' prefs.subjects{sub} '_z+tlrc''[0]'' '];
 end
 unix(cmd)
 

    