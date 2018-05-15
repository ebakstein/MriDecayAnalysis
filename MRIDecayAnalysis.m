%  MRIanalysis
%  analyze T2 maps over ROIs, fit relaxation time (exponential decay)
%  
%  REQUIREMENTS
%   image processing toolbox 
%   statistical toolbox
%   NIFTItools
% 
% 
% Reworked version by E. Bakstein, based on MRIanalysis_3May2018.m
% EDITS 2018-05-15:
% -	Settings: all code parameters editable in the top code section
% -	Loading voxel size from *.nii header directly
% -	ROIS: counts and names now editable in the settings section
% -	Slice selection: reworked + cleaned-up code
% -	Slice selection: divides slices into multiple figures if there are too many to fit on a single figure
% -	Fixed viewing of the results (one figure per slice)

%% INIT
%clear all; close all; 


%% SETTINGS

MRI_IMAGES_PER_FIGURE = [];       % how many MRI slices shall be shown per figure [rows columns] default: enough to plot all slices into one figure
MRI_SLICE_VIEW_ONLY = [];         % select slices to view in SOI selection. E.g. 3:9, set [] to view all
MRI_SLICE_ROTATE_VIEW = [270,90]; % rotate images: [270,90] - surface down, [90,90] surface up, [] - no rotation

ROI_NAMES = {'SINxx','DEXyy'};    % names for ROIS - defines also how many ROIS shall be entered for each slice(just for plotting & reference)

DECAY_TIMES = 8:8:128;            % decay times in the T2map sequence

%% START SCRIPT

%% 1 - DEFINE SOI (on anatomical / T2 scan)

%  LOAD T2 nii-file (filename_vol2 the smaller one)
[file, path] = uigetfile('*.nii','Open the anatomical source file','MultiSelect', 'off');
tic

% DEFAULT_SOI = 7; % CHANGE ::: no default ... always must to be selected by USER
filenameT2 = [path file];
if(~exist(filenameT2,'file'))
    error('file %s not found',filenameT2);
end

data=load_untouch_nii(filenameT2);
image = data.img;  % in the 4D data organization (x,y,z,time) =  x,y - image plate in pixels; z-slice; time=time profile of T2 echo time; 
%
% SLICE SELECTION - all slices at t=1
Nslices = size(image,3);

% select slices + prepare grid
if(~isempty(MRI_SLICE_VIEW_ONLY))
   if(any(sliceNums>Nslices || sliceNums<1))
        error('Some selected slice numbers out of bounds')
   end
   sliceNums = MRI_SLICE_VIEW_ONLY;
else
   sliceNums = 1:Nslices;
end

if(isempty(MRI_IMAGES_PER_FIGURE))
   nRows = round(sqrt(length(sliceNums)));
   nCols = ceil(length(sliceNums)/nRows);           
else
   nRows = MRI_IMAGES_PER_FIGURE(1);
   nCols = MRI_IMAGES_PER_FIGURE(2);           
end

perFig = nRows * nCols;
totFig = ceil(length(sliceNums)/perFig);
tempSOI=[]; 
imlim = [0 max(image(:))]; % same brightness/intensity range for all plots
curA = 1; figN = 0;
% select slices per figure
for ii=1:length(sliceNums) % plotting slices EB: automatic
    j = sliceNums(ii);
    if(mod(curA,perFig)==1) % new figure for slice 1, perFig+1 2*perFig + 1 etc
       curA = 1; figN = figN + 1;
       fh = figure('Name',sprintf('Select slices of interest (SOI): figure %d/%d',figN,totFig));
    end
    
    % do the plotting
    subplot(nRows,nCols,curA); % n_time_points = 16 time points suplot 4x4
   imshow(image(:,:,j,1),imlim) % ,'ImshowInitialMagnification',200
   if(~isempty(MRI_SLICE_ROTATE_VIEW))
    view(MRI_SLICE_ROTATE_VIEW);
   end
   title(sprintf('z=%d',j))   
   
   if(curA == perFig || ii==length(sliceNums))
       questOut = 'No';       
       while(~strcmp(questOut,'Yes'))  % add SOI from given figure (no fig redraw for each SOI as before)
            try
                [~,~]=ginput(1);
                soi = get(get(gca,'title'),'String'); % get clicked axis title
                soi = str2num(soi(3:end));
            catch e
                errordlg('Erroneous or no SOI selected');
                break;
            end
            tempSOI(end+1,1)=soi;
            questOut = questdlg(sprintf('Was this (%d) the last SOI from this set?',soi), 'Confirm', 'Yes','No - continue','Yes'); 
       end % select SOI per given figure
       if(ishandle(fh))
        close(fh)
       end
   end
   curA = curA + 1; 
end % iterate over all slices

results = struct();
results.filenameT2_w = file;
results.path = path;
results.soi = tempSOI;

clear fh j x y r questOut soi tempSOI DATAPATH FILE filenameT2;

%% 2 - ENTER ALL ROI (all ROI types for each slice)
% go through all requested roi
    % define Region of interest - ROI, T2, selected slice t=1   
    % ROI elliptical/polygonal selector - image processing toolbox required
    
    % INSTRUCTIONS: 
    %      1] create polygon by clicking on the image, 
    %      2] close region by clicking on the initial point (circle mouse pointer will be shown)
    %      NOTE: ROI shape can be modified by dragging anchor points      
    %      3] confirm ROI by double-click inside the ROI area
    %

% SOI
for j=1:length(results.soi)
    % ROI TYPE
    for si=1:length(ROI_NAMES)
        roiName = ROI_NAMES{si};
        questOut = 'No';
        
        while(~strcmp(questOut,'Yes'))
            % Show figure - selected slice 
            fh = figure; 
            imshow(image(:,:,results.soi(j,1),1),[]);  % can we turn image 90 o clockwise here ?and calculate all following steps from a such ROI ?
            title(sprintf('slice %d, t=1, select ROI for %s',results.soi(j,1),roiName))

            % rotate view - if set
            if(~isempty(MRI_SLICE_ROTATE_VIEW))
                view(MRI_SLICE_ROTATE_VIEW);
            end

            % GRAPHICAL input of ROI
            clear h
            h=impoly();          % input polygonal mask. imellipse(), imfreehand() and other can be used instead
            position = wait(h);  % wait for the user to input mask
            msk = createMask(h); % convert to logical mask - same size as orig. img.    

            % DIALOG: verify that mask is ok, or repeat
            hold on; [ex,ey]=find(edge(msk)); plot(ey,ex,'r.','MarkerSize',1);
            questOut = questdlg(sprintf('Is inserted %s ROI ok?',roiName), 'Confirm ROI', 'Yes','No - repeat','Yes');    
            close(fh);
        end        
        results.(roiName)(j).roi = msk;
    end   
      temp=[];
      [ex, ey] = find(edge(results.SINxx(j).roi));
      results.SINxx(j).msk=[ex, ey];
      temp(1,j)=sum(sum(image(ex(1:end),ey(1:end),results.soi(j,1))));
      clear ex ey; 
      [ex, ey] = find(edge(results.DEXyy(j).roi));
      results.DEXyy(j).msk=[ex, ey];
      temp(2,j)=sum(sum(image(ex(1:end),ey(1:end),results.soi(j,1))));
      results.SINxx(j).T2intens=temp(1,j);
      results.DEXyy(j).T2intens=temp(2,j);
      results.T2wRatio(j,1)=temp(2,j)/temp(1,j);
end
           
    

close all; clear data ex ey fh h image imlim j msk Nslices position questOut roiName si temp
%% load T2map image nii-file (~ 12MB for pups)
[file, path] = uigetfile('*.nii','Open the T2 map source file','MultiSelect', 'off');
tic
filenameT2map = [path file];

if(~exist(filenameT2map,'file'))
    error('file %s not found',filenameT2map);
end
results.filenameT2_map = file;
data1=load_untouch_nii(filenameT2map);
hdr = load_nii_hdr(filenameT2map);

image = data1.img;
Ntimes = size(image,4);
imlim = [0 max(image(:))];
% 
for j=1:length(results.soi)
    for si=1:length(ROI_NAMES)
        roiName = ROI_NAMES{si};

    % go through time sequence, compute mean intensity within current mask
    profile = nan(1,Ntimes);
    for ii = 1:Ntimes
       imgsub = image(:,:,results.soi(j,1),ii); 
       profile(1,ii) = mean(imgsub(results.(roiName)(j).roi)); 
    end

   % FIT EXPONENTIAL TO CURRENT TIME PROFILE, starting from index 2 "y = a + b*exp(c*x)"
    yy = profile; %(2:end); % fix to tak all 16 echo times
    xx = DECAY_TIMES; %1:Ntimes-1; % TE = 8
    
    % function to fit
   
    decay_fit = fit(xx', yy', 'exp1');
    decay_par = coeffvalues(decay_fit);
    decay_t = -1/decay_par(2);
  
      results.(roiName)(j).expT = decay_t; 
      results.(roiName)(j).expParams = decay_par;
      results.(roiName)(j).meanRoiT2 = profile;
    end
end

% volume calculation

% voxel volume - based on *.nii header
V = prod(hdr.dime.pixdim(1:3))/(hdr.dime.pixdim(4)^2);

% original - manual version
% p=1.5625/10; % voxel size (see image properties with ImageJ) 
% V=p*p*1; %here "1" is the tickness of slice in mm

for i=1:length(results.soi)
    cacca1=find(results.SINxx(i).roi);
    cacca2=find(results.DEXyy(i).roi);
    results.SINxx(i).vol=V*size(cacca1,1);
    results.DEXyy(i).vol=V*size(cacca2,1);
    clear cacca1 cacca2
end

% store results for given subjectmatPath = [DATAPATH FILE 'T2_results.mat'];
matPath = [path file '_results.mat'];
save(matPath,'results') % store the image as well?
fprintf('results saved to %s\n',matPath);
% close all; 
% clear data1 DATAPATH decay_fit  decay_par decay_t FILE filenameT2 i ii image imlim imgsub j matPath Ntimes p profile roiName roiNames si V xx yy;



%% visualize perf. profiles + exp. fit for all ROI
Nslices = size(image,4);
Nroi = length(ROI_NAMES);
imgmax = max(image(:)); % to plot all images in the same range

for si = 1:length(results.soi) % SOI_    
    figure('Name',sprintf('RESULTS: slice %d',results.soi(si)));  
    for ri=1:length(ROI_NAMES) % ROI types (for each slice)
        roiName = ROI_NAMES{ri};
        res = results.(roiName)(si);

        % one roi on one row
        splt = (ri-1)*3+1; % initial subplot

        % A) initial image at T=1
        subplot(Nroi,3,splt)       

    %     img = image(:,:,res.soi,1)/imgmax;
    %     imrgb = imoverlay(img, roiEdge,[1 0 0]);
    %     imshow(imrgb)

        img = image(:,:,results.soi(1),1);
        imshow(img,[0 imgmax]); 
       % rotate view - if set
        if(~isempty(MRI_SLICE_ROTATE_VIEW))
            view(MRI_SLICE_ROTATE_VIEW);
        end
        hold on;    
        title(sprintf('%s, T=1',roiName));
        [ex ey] = find(edge(res.roi));
        plot(ey,ex,'c.','MarkerSize',1)
        xlabel('decay time [ms]')

        
        % B) - temporal profile
        subplot(Nroi,3,splt+1)
        
        % T2 time profile
        % show exponential fit from T=2
         x = DECAY_TIMES; %x = 1:Ntimes-1;
        plot(x,res.meanRoiT2); hold on; 
         %y =  res.expParams(1) + res.expParams(2) * exp( -x/res.expParams);  %PROBLEM HERE ????
        y = res.expParams(1)*exp(res.expParams(2)*x);
        plot(x, y,'r.'); 
        
        %text(x(4)+2,y(4),sprintf('t = %.3f',res.expT));
        title(sprintf('%s time profile',roiName))
        axis square
        %xlim([2 17])

        % C)  plot last image in sequence
        subplot(Nroi,3,splt+2)
        imshow(image(:,:,si,end),[0,imgmax]); 
         if(~isempty(MRI_SLICE_ROTATE_VIEW))
            view(MRI_SLICE_ROTATE_VIEW);
         end
        hold on;
        plot(ey,ex,'c.','MarkerSize',1)

        title(sprintf('%s, T=%d',roiName,Nslices));
        xlabel('decay time [ms]')
    end
end
%  button=questdlg('Continue?','OK'); 
