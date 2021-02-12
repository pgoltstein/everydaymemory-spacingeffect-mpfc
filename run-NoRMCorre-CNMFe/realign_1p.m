function realign_1p( varargin )
% function realign_1p
% Realigns one-photon (epifluorescence) calcium imaging stacks
%
% Dependencies:
% - https://github.com/simonsfoundation/NoRMCorre.git
%
% Pieter Goltstein & Annet Glas
% October 31, 2017 - Version 0.1
%

    % Check for input parameters
    if nargin == 1
        ForceRealign = varargin{1};
    else
        ForceRealign = false;
    end

    % Parameters non-rigid alignment
    gSigma_align = 9; % width of the gaussian kernel
    gSize_align = 25; % maximum diameter of blobs

    % Or automatically load the first stack
    ImFiles = dir('*.tif');
    ImFile = ImFiles(1).name;

    % Do it
    AlignedStackName = 'AlignedStack';
    fprintf('\n---------------- realign_1p ----------------\n\n');
    fprintf('Datadir: %s\n', pwd);

    %% Settings
    meta = struct;
    meta.frameRate = 10; % Hz
    meta.xRes = 256; % pixels
    meta.yRes = 256; % pixels


    %% Load imaging stack and meta data

    fprintf('\nLoading imaging stack: %s\n', ImFile);
    warning('off'); % Tiff class gives stupid warnings because of missing fields

    TiffInfo = Tiff(ImFile);
    ImInfo = imfinfo(ImFile);
    meta.orig_nFrames = length(ImInfo);
    meta.orig_xRes = ImInfo(1).Width;
    meta.orig_yRes = ImInfo(1).Height;
    miniscopeProperties = regexp(ImInfo(1).ImageDescription,'\d+','match');
    meta.exposureTime = str2double(miniscopeProperties{2});
    meta.origFrameRate = 1000/meta.exposureTime;
    meta.LEDpower = str2double(miniscopeProperties{5});
    fprintf('Exposure time: %d ms, frame rate: %d Hz, LED power: %d%% \n', ...
        meta.exposureTime, meta.origFrameRate, meta.LEDpower);
    fprintf('Image dimensions: x=%d, y=%d, t=%d\n', ...
        meta.orig_xRes, meta.orig_yRes, meta.orig_nFrames );
    ImData = zeros( meta.orig_yRes, meta.orig_xRes, meta.orig_nFrames , 'uint16' );
    meta.TimeStamps = zeros(1,meta.orig_nFrames );
    fprintf('Loading frame: %6d',0);
    for f = 1:meta.orig_nFrames
        fprintf('\b\b\b\b\b\b%6d',f);
        TiffInfo.setDirectory(f);
        ImData(:,:,f) = TiffInfo.read();
        miniscopeProperties = regexp(ImInfo(f).ImageDescription,'\d+','match');
        meta.TimeStamps(f) = str2double(miniscopeProperties{4});
    end
    fprintf(' ... done\n');
    warning('on'); % Aaaand warnings back on..

    % Convert to single
    ImData = single(ImData);


    %% Interpolate dropped frames
    avgFrameBrightness = squeeze(mean(mean(ImData,1),2));
    medianBrightness = median(avgFrameBrightness);
    meta.droppedFrames = find(avgFrameBrightness<(medianBrightness/10));
    fprintf('\nInterpolating %d dropped frames\n', length(meta.droppedFrames));
    for f = meta.droppedFrames
        ImData(:,:,f) = ImData(:,:,f-1);
    end


    %% Spatial frequency filter stack
    psf = fspecial('gaussian', round(gSize_align), gSigma_align);
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;
    fImData = imfilter(ImData, psf, 'symmetric');


    %% Align stack (NoRMCorre)
    warning('off'); % NoRMCorre gives warnings because of temporary variable in parfor loop
    if ~exist('ImageRegistration.mat','file') || ForceRealign == true
        fprintf('\nPerforming non-rigid image registration using NoRMCorre\n');
        NoRMCorre_options = NoRMCorreSetParms(...
            'd1', meta.orig_yRes, 'd2', meta.orig_xRes, 'bin_width', 50, 'init_batch', 200,...
            'grid_size', [256,256], 'mot_uf', 4, 'correct_bidir', false, ...
            'overlap_pre', 32, 'overlap_post', 32, 'max_shift', 20, 'max_dev', [5,5]);
        [fImData,regShifts,regTemplate] = normcorre_batch( fImData, NoRMCorre_options );
        save('ImageRegistration.mat','regShifts','regTemplate','NoRMCorre_options');
        fprintf(' done\nSaved: ImageRegistration.mat\n');
    else
        fprintf('\nPerforming non-rigid image registration using NoRMCorre\n');
        load('ImageRegistration.mat');
        fprintf('Loaded: ImageRegistration.mat\n');
    end

    fprintf('\nApplying non-rigid image registration using NoRMCorre\n');
    ImData = apply_shifts(ImData,regShifts,NoRMCorre_options);
    warning('on'); % Aaaand warnings back on..
    clear fImData;


    %% Quality check
    nZerosPerFrame = zeros(1,meta.orig_nFrames);
    fprintf('Quality check, frame: %6d',0);
    for f = 1:meta.orig_nFrames
        fprintf('\b\b\b\b\b\b%6d',f);
        nZerosPerFrame(f) = sum(sum(ImData(:,:,f)==0));
    end
    fprintf(' ... done\n');
    meta.crappyFrames = find(nZerosPerFrame>100);
    meta.goodFrames = find(nZerosPerFrame<=100);
    fprintf('\nInterpolating %d crappy frames\n', length(meta.crappyFrames));
    for f = meta.crappyFrames
        if f < meta.goodFrames(1)
            ImData(:,:,f) = ImData(:,:,meta.goodFrames(1));
        else
            ImData(:,:,f) = ImData(:,:,f-1);
        end
    end


    %% Downsample data, uint16 and save to *.mat
    meta.nFrames  = ceil( meta.orig_nFrames * (meta.frameRate / meta.origFrameRate) );
    newFrameList = round(linspace(1,meta.nFrames,meta.orig_nFrames));
    Y = zeros(meta.yRes, meta.xRes, meta.nFrames, 'uint16');
    fprintf('Spatial & temporal downsampling, frame: %6d',0);
    for f = 1:meta.nFrames
        fprintf('\b\b\b\b\b\b%6d',f);
        Y(:,:,f) = uint16(imresize(mean(ImData(:,:,newFrameList==f),3),[meta.yRes,meta.xRes]));
    end
    fprintf(' ... done\n');
    clear ImData;
    Ysiz = [meta.yRes, meta.xRes, meta.nFrames]';
    save([AlignedStackName '.mat',], 'Y', 'Ysiz', 'meta', '-v7.3');

end
