function preprocessing_1p_concat_sessions
% function preprocessing_1p_concat_sessions
% Processes one-photon (epifluorescence) calcium imaging stacks, runs
%   on multiple processed stacks
%
% Dependencies:
% - https://github.com/simonsfoundation/NoRMCorre.git
% - https://github.com/zhoupc/CNMF_E.git
%   >> cnmfe_setup (loads path for CNMF_E -> run once and save path)
%
% Pieter Goltstein & Annet Glas
% October 24, 2018 - Version 0.1
%

    % Parameters non-rigid alignment
    find_subdirectories_by = '201*';
    gSigma_align = 9; % width of the gaussian kernel
    gSize_align = 25; % maximum diameter of blobs

    % Parameters CNMF_E
    gSigma_cnmfe = 2;        % width of the gaussian kernel, which can approximates the average neuron shape
    gSize_cnmfe = 12;        % maximum diameter of neurons in the image plane. larger values are preferred.
    min_corr_cnmfe = 0.75;  % minimum local correlation for a seeding pixel
    min_pnr_cnmfe = 8;       % minimum peak-to-noise ratio for a seeding pixel
    edge_margin_cnmfe = 10;   % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
    AlignedStackName = 'AlignedStack';

    % Or automatically load the processed stacks
    BasePathName = pwd;

    % Get filenames of processed stacks from subdirectories
    subdirs = dir(find_subdirectories_by);

    fprintf('\n---- Running alignment on individual stacks ----\n\n');
    fprintf('Datadir: %s\n', pwd);
    fprintf('Date/time: %s\n', datestr(clock));

    % First get aligned stacks per subdir
    for dirnr = 1:length(subdirs)
        cd(BasePathName);
        cd(subdirs(dirnr).name);
        Align = dir('AlignedStack.mat');
        if ~isempty(Align)
            fprintf('\nRealignment & resampling has already been done\n');
        else
            realign_1p;
        end
    end

    fprintf('\n---- Concatenating and aligning stack ----\n\n');
    fprintf('Datadir: %s\n', pwd);
    fprintf('Date/time: %s\n', datestr(clock));

    % Loop subdirs and load data separately
    nFrames = 0;
    SubDirNames = {};
    for dirnr = 1:length(subdirs)
        cd(BasePathName);
        cd(subdirs(dirnr).name);
        SubDirNames{dirnr} = subdirs(dirnr).name;
        stackdata{dirnr} = load([AlignedStackName '.mat'],'meta','Ysiz','Y');
        fprintf('Loaded data from: %s\n',[AlignedStackName '.mat']);
        nFrames = nFrames + stackdata{dirnr}.meta.nFrames;
    end
    cd(BasePathName);

    % New meta data
    meta.yRes = stackdata{1}.meta.yRes;
    meta.xRes = stackdata{1}.meta.xRes;
    meta.nFrames = nFrames;
    meta.frameRate = stackdata{1}.meta.frameRate;

    % Create one big image data matrix and fill it with all
    Y = zeros(meta.yRes, meta.xRes, nFrames, 'uint16');
    startIx = 1;
    SessionFrames = {};
    for nr = 1:length(stackdata)
        SessionFrames{nr} = [startIx (startIx+stackdata{nr}.meta.nFrames-1)];
        Y(:,:,startIx:(startIx+stackdata{nr}.meta.nFrames-1)) = stackdata{nr}.Y;
        startIx = startIx+stackdata{nr}.meta.nFrames;
    end
    clear stackdata;

    %% Spatial frequency filter stack
    psf = fspecial('gaussian', round(gSize_align), gSigma_align);
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;
    fY = imfilter(Y, psf, 'symmetric');

    %% Align the full concatenated stack (NoRMCorre)
    warning('off'); % NoRMCorre gives warnings because of temporary variable in parfor loop
    if ~exist('ImageRegistration.mat','file')
        fprintf('\nPerforming non-rigid image registration of the concatenated stack using NoRMCorre\n');
        NoRMCorre_options = NoRMCorreSetParms(...
            'd1', meta.yRes, 'd2', meta.xRes, 'bin_width', 50, 'init_batch', 200,...
            'grid_size', [256,256], 'mot_uf', 4, 'correct_bidir', false, ...
            'overlap_pre', 32, 'overlap_post', 32, 'max_shift', 20, 'max_dev', [5,5]);
        [fY,regShifts,regTemplate] = normcorre_batch( fY, NoRMCorre_options );
        save('ImageRegistration.mat','regShifts','regTemplate','NoRMCorre_options');
        fprintf(' done\nSaved: ImageRegistration.mat\n');
    else
        fprintf('\nPerforming non-rigid image registration using NoRMCorre\n');
        load('ImageRegistration.mat');
        fprintf('Loaded: ImageRegistration.mat\n');
    end

    fprintf('\Applying non-rigid image registration using NoRMCorre\n');
    Y = apply_shifts(Y,regShifts,NoRMCorre_options);
    warning('on'); % Warnings back on..
    clear fY;

    % Save the new concatenated processed stack
    Ysiz = [meta.yRes, meta.xRes, meta.nFrames]';
    save([AlignedStackName '.mat',], 'Y', 'Ysiz', 'meta', '-v7.3');
    clear Y;


    %% Source extraction (CNMF_E)

    fprintf('\n---- Source extraction ----\n\n');
    fprintf('Datadir: %s\n', pwd);
    fprintf('Date/time: %s\n', datestr(clock));

    % Select data
    data = matfile([AlignedStackName '.mat',]);

    % Options and settings
    neuron = Sources2D('d1',meta.yRes,'d2',meta.xRes, ...
        'ssub', 1, 'tsub', 1, 'gSig', gSigma_cnmfe, 'gSiz', gSize_cnmfe);
    neuron.Fs = meta.frameRate;
    neuron.options.search_method = 'ellipse';
    neuron.options.dist = 3;
    neuron.options.deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
        'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
        'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
        'optimize_pars', true, ...  % optimize AR coefficients
        'optimize_b', true, ...% optimize the baseline);
        'max_tau', 100);    % maximum decay time (unit: frame);

    % Load data
    Y = double(data.Y);
    Y = neuron.reshape(Y, 1);       % convert a 3D video into a 2D matrix
    fprintf('Data re-loaded, image dimensions: x=%d, y=%d, t=%d\n', ...
        meta.xRes, meta.yRes, meta.nFrames);


    %% compute correlation image and peak-to-noise ratio image.
    [Cn, pnr] = neuron.correlation_pnr(Y(:, round(linspace(1, meta.nFrames, min(meta.nFrames, 1000)))));

    % show correlation image
    figure('position', [10, 500, 1776, 400]);
    subplot(131);
    imagesc(Cn, [0, 1]); colorbar;
    axis equal off tight;
    title('correlation image');

    % show peak-to-noise ratio
    subplot(132);
    imagesc(pnr,[0,max(pnr(:))*0.98]); colorbar;
    axis equal off tight;
    title('peak-to-noise ratio');

    % show pointwise product of correlation image and peak-to-noise ratio
    subplot(133);
    imagesc(Cn.*pnr, [0,max(pnr(:))*0.98]); colorbar;
    axis equal off tight;
    title('Cn*PNR');
    drawnow;


    %% initialization of A, C
    debug_on = false;     % visualize the initialization procedue.
    save_avi = false;    % save the initialization procedure as an avi movie.
    patch_par = [1,1]*1; % divide the optical field into m X n patches and do initialization patch by patch. It can be used when the data is too large
    K = [];              % maximum number of neurons to search within each patch. you can use [] to search the number automatically

    min_pixel = gSigma_cnmfe^2;  % minimum number of nonzero pixels for each neuron
    neuron.updateParams('min_corr', min_corr_cnmfe, 'min_pnr', min_pnr_cnmfe, ...
        'min_pixel', min_pixel, 'bd', edge_margin_cnmfe);
    neuron.options.nk = 5;  % number of knots for detrending
    merge_thr = [.075, 0.75, 0];     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
    dmin = 1;

    % greedy method for initialization
    tic;
    [center, Cn, pnr] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi);
    fprintf('Time cost in initializing neurons: %.2f seconds\n', toc);

    % show results
    figure;
    imagesc(Cn, [0.1, 0.95]);
    hold on; plot(center(:, 2), center(:, 1), 'or');
    colormap; axis off tight equal;
    drawnow;

    % sort neurons
    neuron.orderROIs('snr');
    neuron_init = neuron.copy();


    %% iteratively update A, C and B
    % parameters, merge neurons
    display_merge = false;          % visually check the merged neurons
    view_neurons = false;           % view all neurons

    % parameters, estimate the background
    spatial_ds_factor = 1;      % spatial downsampling factor. it's for faster estimation
    thresh = 10;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)

    bg_neuron_ratio = 1.5;  % spatial range / diameter of neurons

    % parameters, estimate the spatial components
    update_spatial_method = 'hals';  % the method for updating spatial components {'hals', 'hals_thresh', 'nnls', 'lars'}
    Nspatial = 5;       % this variable has different meanings:
                        %1) udpate_spatial_method=='hals' or 'hals_thresh',
                        %then Nspatial is the maximum iteration
                        %2) update_spatial_method== 'nnls', it is the maximum
                        %number of neurons overlapping at one pixel

    % parameters for running iteratiosn
    nC = size(neuron.C, 1);    % number of neurons

    maxIter = 2;        % maximum number of iterations
    miter = 1;
    while miter <= maxIter

        % merge neurons
%         cnmfe_quick_merge;              % run neuron merges
        neuron_bk = neuron.copy();
        [merged_ROI, newIDs] = neuron.quickMerge(merge_thr);  % merge neurons based on the correlation computed with {'A', 'S', 'C'}
        % A: spatial shapes; S: spike counts; C: calcium traces
        cols = [0, 0, 1; 0, 1, 0; 0, 1, 1; 1, 0, 0; 1, 0, 1; 1, 1, 0; 1, 1, 1];
        if display_merge && ~isempty(merged_ROI)
            figure('position', [1,1, 1200, 600]);
            ind_before = false(size(neuron_bk.A, 2), 1);
            ind_after = false(size(neuron.A, 2), 1);
            m = 1;
            while m<=length(merged_ROI)
                subplot(221);
                [tmp_img, col, ~] = neuron_bk.overlapA(merged_ROI{m});
                imagesc(tmp_img);
                axis equal off tight;
                subplot(222);
                imagesc(tmp_img);
                axis equal off tight;
                [tmp_r, tmp_c, ~] = find(sum(tmp_img, 3)>0);
                xlim([min(tmp_c)-10, max(tmp_c)+10]);
                ylim([min(tmp_r)-10, max(tmp_r)+10]);
                %         neuron.image(sum(Ain(:, merged_ROI{m}), 2));
                axis off;
                subplot(2,2,3:4); cla;
                tmp_C = neuron_bk.C_raw(merged_ROI{m}, :);
                %         tmp_C = bsxfun(@times, tmp_C, 1./max(tmp_C, [], 1));
                for mm=1:size(tmp_C, 1)
                    hold on;
                    plot(tmp_C(mm,:), 'color', cols(col(mm), :),  'linewidth', 2);
                end
                temp = input('keep this merge? (y(default)/n(cancel)/b(back))/e(end)   ', 's');
                if strcmpi(temp, 'n')
                    ind_after(newIDs(m)) = true;
                    ind_before(merged_ROI{m}) = true;
                    m = m+1;
                elseif strcmpi(temp, 'b')
                    m = m-1;
                elseif strcmpi(temp, 'e')
                    break;
                else
                    m = m+1;
                end
            end

            neuron.A = [neuron.A(:, ~ind_after), neuron_bk.A(:, ind_before)];
            neuron.C = [neuron.C(~ind_after, :); neuron_bk.C(ind_before, :)];
            neuron.C_raw = [neuron.C_raw(~ind_after, :); neuron_bk.C_raw(ind_before, :)];
            neuron.S = [neuron.S(~ind_after, :); neuron_bk.S(ind_before, :)];
            neuron.P.kernel_pars = [neuron.P.kernel_pars(~ind_after, :); neuron_bk.P.kernel_pars(ind_before, :)];
            close;
        end

        % view neurons
        if view_neurons
            neuron.viewNeurons([], neuron.C_raw);
        end


%         cnmfe_merge_neighbors;          % merge neurons if two neurons' peak pixels are too close
        neuron_bk = neuron.copy();
        if ~exist('center_method', 'var')
            center_method = 'max';
        end
        [merged_ROI, newIDs] = neuron.MergeNeighbors(dmin, center_method);
        cols = [0, 0, 1; 0, 1, 0; 0, 1, 1; 1, 0, 0; 1, 0, 1; 1, 1, 0; 1, 1, 1];
        if display_merge && ~isempty(merged_ROI)
            figure('position', [1,1, 1200, 600]);
            ind_before = false(size(neuron_bk.A, 2), 1);
            ind_after = false(size(neuron.A, 2), 1);
            m = 1;
            while m<=length(merged_ROI)
                subplot(221);
                [tmp_img, col, ~] = neuron_bk.overlapA(merged_ROI{m});
                imagesc(tmp_img);
                axis equal off tight;
                subplot(222);
                imagesc(tmp_img);
                axis equal off tight;
                [tmp_r, tmp_c, ~] = find(sum(tmp_img, 3)>0);
                xlim([min(tmp_c)-10, max(tmp_c)+10]);
                ylim([min(tmp_r)-10, max(tmp_r)+10]);
                %         neuron.image(sum(Ain(:, merged_ROI{m}), 2));
                axis off;
                subplot(2,2,3:4); cla;
                tmp_C = neuron_bk.C_raw(merged_ROI{m}, :);
                %         tmp_C = bsxfun(@times, tmp_C, 1./max(tmp_C, [], 1));
                for mm=1:size(tmp_C, 1)
                    hold on;
                    plot(tmp_C(mm,:), 'color', cols(col(mm), :),  'linewidth', 2);
                end
                temp = input('keep this merge? (y(default)/n(cancel)/b(back))/e(end)   ', 's');
                if strcmpi(temp, 'n')
                    ind_after(newIDs(m)) = true;
                    ind_before(merged_ROI{m}) = true;
                    m = m+1;
                elseif strcmpi(temp, 'b')
                    m = m-1;
                elseif strcmpi(temp, 'e')
                    break;
                else
                    m = m+1;
                end
            end

            neuron.A = [neuron.A(:, ~ind_after), neuron_bk.A(:, ind_before)];
            neuron.C = [neuron.C(~ind_after, :); neuron_bk.C(ind_before, :)];
            neuron.C_raw = [neuron.C_raw(~ind_after, :); neuron_bk.C_raw(ind_before, :)];
            neuron.S = [neuron.S(~ind_after, :); neuron_bk.S(ind_before, :)];
            neuron.P.kernel_pars = [neuron.P.kernel_pars(~ind_after, :); neuron_bk.P.kernel_pars(ind_before, :)];
            close;
        end

        % view neurons
        if view_neurons
            neuron.viewNeurons([], neuron.C_raw);
        end



        % update background
        % estimate the background
        tic;
%         cnmfe_update_BG;
        clear Ysignal;
        tic;
        Ybg = Y-neuron.A*neuron.C;
        rr = ceil(neuron.options.gSiz * bg_neuron_ratio);
        active_px = []; %(sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
        [Ybg, Ybg_weights] = neuron.localBG(Ybg, spatial_ds_factor, rr, active_px, neuron.P.sn, thresh); % estiamte local background.
        % subtract the background from the raw data.
        Ysignal = Y - Ybg;

        % estimate noise
        if ~isfield(neuron.P, 'sn') || isempty(neuron.P.sn)
            % estimate the noise for all pixels
            b0 =zeros(size(Ysignal,1), 1);
            sn = b0;
            parfor m=1:size(neuron.A,1)
                [b0(m), sn(m)] = estimate_baseline_noise(Ysignal(m, :));
            end
            Ysignal = bsxfun(@minus, Ysignal, b0);
            neuron.P.sn = sn;
        end
        fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
        % neuron.playMovie(Ysignal); % play the video data after subtracting the background components.

        % update spatial & temporal components
        tic;
        for m=1:2
            %temporal
            neuron.updateTemporal_endoscope(Ysignal);

%             cnmfe_quick_merge;              % run neuron merges
            neuron_bk = neuron.copy();
            [merged_ROI, newIDs] = neuron.quickMerge(merge_thr);  % merge neurons based on the correlation computed with {'A', 'S', 'C'}
            % A: spatial shapes; S: spike counts; C: calcium traces
            cols = [0, 0, 1; 0, 1, 0; 0, 1, 1; 1, 0, 0; 1, 0, 1; 1, 1, 0; 1, 1, 1];
            if display_merge && ~isempty(merged_ROI)
                figure('position', [1,1, 1200, 600]);
                ind_before = false(size(neuron_bk.A, 2), 1);
                ind_after = false(size(neuron.A, 2), 1);
                m = 1;
                while m<=length(merged_ROI)
                    subplot(221);
                    [tmp_img, col, ~] = neuron_bk.overlapA(merged_ROI{m});
                    imagesc(tmp_img);
                    axis equal off tight;
                    subplot(222);
                    imagesc(tmp_img);
                    axis equal off tight;
                    [tmp_r, tmp_c, ~] = find(sum(tmp_img, 3)>0);
                    xlim([min(tmp_c)-10, max(tmp_c)+10]);
                    ylim([min(tmp_r)-10, max(tmp_r)+10]);
                    %         neuron.image(sum(Ain(:, merged_ROI{m}), 2));
                    axis off;
                    subplot(2,2,3:4); cla;
                    tmp_C = neuron_bk.C_raw(merged_ROI{m}, :);
                    %         tmp_C = bsxfun(@times, tmp_C, 1./max(tmp_C, [], 1));
                    for mm=1:size(tmp_C, 1)
                        hold on;
                        plot(tmp_C(mm,:), 'color', cols(col(mm), :),  'linewidth', 2);
                    end
                    temp = input('keep this merge? (y(default)/n(cancel)/b(back))/e(end)   ', 's');
                    if strcmpi(temp, 'n')
                        ind_after(newIDs(m)) = true;
                        ind_before(merged_ROI{m}) = true;
                        m = m+1;
                    elseif strcmpi(temp, 'b')
                        m = m-1;
                    elseif strcmpi(temp, 'e')
                        break;
                    else
                        m = m+1;
                    end
                end

                neuron.A = [neuron.A(:, ~ind_after), neuron_bk.A(:, ind_before)];
                neuron.C = [neuron.C(~ind_after, :); neuron_bk.C(ind_before, :)];
                neuron.C_raw = [neuron.C_raw(~ind_after, :); neuron_bk.C_raw(ind_before, :)];
                neuron.S = [neuron.S(~ind_after, :); neuron_bk.S(ind_before, :)];
                neuron.P.kernel_pars = [neuron.P.kernel_pars(~ind_after, :); neuron_bk.P.kernel_pars(ind_before, :)];
                close;
            end

            % view neurons
            if view_neurons
                neuron.viewNeurons([], neuron.C_raw);
            end

            %spatial
            neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);

            % post process the spatial components (you can run either of these two operations, or both of them)
            neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
            neuron.compactSpatial();    % run this line if neuron shapes are circular

%             cnmfe_merge_neighbors;
            neuron_bk = neuron.copy();
            if ~exist('center_method', 'var')
                center_method = 'max';
            end
            [merged_ROI, newIDs] = neuron.MergeNeighbors(dmin, center_method);
            cols = [0, 0, 1; 0, 1, 0; 0, 1, 1; 1, 0, 0; 1, 0, 1; 1, 1, 0; 1, 1, 1];
            if display_merge && ~isempty(merged_ROI)
                figure('position', [1,1, 1200, 600]);
                ind_before = false(size(neuron_bk.A, 2), 1);
                ind_after = false(size(neuron.A, 2), 1);
                m = 1;
                while m<=length(merged_ROI)
                    subplot(221);
                    [tmp_img, col, ~] = neuron_bk.overlapA(merged_ROI{m});
                    imagesc(tmp_img);
                    axis equal off tight;
                    subplot(222);
                    imagesc(tmp_img);
                    axis equal off tight;
                    [tmp_r, tmp_c, ~] = find(sum(tmp_img, 3)>0);
                    xlim([min(tmp_c)-10, max(tmp_c)+10]);
                    ylim([min(tmp_r)-10, max(tmp_r)+10]);
                    %         neuron.image(sum(Ain(:, merged_ROI{m}), 2));
                    axis off;
                    subplot(2,2,3:4); cla;
                    tmp_C = neuron_bk.C_raw(merged_ROI{m}, :);
                    %         tmp_C = bsxfun(@times, tmp_C, 1./max(tmp_C, [], 1));
                    for mm=1:size(tmp_C, 1)
                        hold on;
                        plot(tmp_C(mm,:), 'color', cols(col(mm), :),  'linewidth', 2);
                    end
                    temp = input('keep this merge? (y(default)/n(cancel)/b(back))/e(end)   ', 's');
                    if strcmpi(temp, 'n')
                        ind_after(newIDs(m)) = true;
                        ind_before(merged_ROI{m}) = true;
                        m = m+1;
                    elseif strcmpi(temp, 'b')
                        m = m-1;
                    elseif strcmpi(temp, 'e')
                        break;
                    else
                        m = m+1;
                    end
                end

                neuron.A = [neuron.A(:, ~ind_after), neuron_bk.A(:, ind_before)];
                neuron.C = [neuron.C(~ind_after, :); neuron_bk.C(ind_before, :)];
                neuron.C_raw = [neuron.C_raw(~ind_after, :); neuron_bk.C_raw(ind_before, :)];
                neuron.S = [neuron.S(~ind_after, :); neuron_bk.S(ind_before, :)];
                neuron.P.kernel_pars = [neuron.P.kernel_pars(~ind_after, :); neuron_bk.P.kernel_pars(ind_before, :)];
                close;
            end

            % view neurons
            if view_neurons
                neuron.viewNeurons([], neuron.C_raw);
            end

            % stop the iteration when neuron numbers are unchanged.
            if isempty(merged_ROI)
                break;
            end
        end
        fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);

        % pick neurons from the residual (cell 4).
        if miter==1
            seed_method = 'auto';
            [center_new, Cn_res, pnr_res] = neuron.pickNeurons(Ysignal - neuron.A*neuron.C, patch_par, seed_method, debug_on); % method can be either 'auto' or 'manual'
        end

        % stop the iteration
        temp = size(neuron.C, 1);
        if or(nC==temp, miter==maxIter)
            break;
        else
            miter = miter+1;
            nC = temp;
        end
    end


    %% delete some neurons and run CNMF-E iteration
    neuron.orderROIs('decay_time');  % you can also use {'snr', 'mean', 'decay_time'}

    % Manual control over which neurons to keep
%     neuron.viewNeurons([], neuron.C_raw);

    tic;
%     cnmfe_update_BG;
    clear Ysignal;
    tic;
    Ybg = Y-neuron.A*neuron.C;
    rr = ceil(neuron.options.gSiz * bg_neuron_ratio);
    active_px = []; %(sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
    [Ybg, Ybg_weights] = neuron.localBG(Ybg, spatial_ds_factor, rr, active_px, neuron.P.sn, thresh); % estiamte local background.
    % subtract the background from the raw data.
    Ysignal = Y - Ybg;

    % estimate noise
    if ~isfield(neuron.P, 'sn') || isempty(neuron.P.sn)
        % estimate the noise for all pixels
        b0 =zeros(size(Ysignal,1), 1);
        sn = b0;
        parfor m=1:size(neuron.A,1)
            [b0(m), sn(m)] = estimate_baseline_noise(Ysignal(m, :));
        end
        Ysignal = bsxfun(@minus, Ysignal, b0);
        neuron.P.sn = sn;
    end


    fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
    %update spatial & temporal components
    tic;
    for m=1:2
        %temporal
        neuron.updateTemporal_endoscope(Ysignal);
%         cnmfe_quick_merge;              % run neuron merges
        neuron_bk = neuron.copy();
        [merged_ROI, newIDs] = neuron.quickMerge(merge_thr);  % merge neurons based on the correlation computed with {'A', 'S', 'C'}
        % A: spatial shapes; S: spike counts; C: calcium traces
        cols = [0, 0, 1; 0, 1, 0; 0, 1, 1; 1, 0, 0; 1, 0, 1; 1, 1, 0; 1, 1, 1];
        if display_merge && ~isempty(merged_ROI)
            figure('position', [1,1, 1200, 600]);
            ind_before = false(size(neuron_bk.A, 2), 1);
            ind_after = false(size(neuron.A, 2), 1);
            m = 1;
            while m<=length(merged_ROI)
                subplot(221);
                [tmp_img, col, ~] = neuron_bk.overlapA(merged_ROI{m});
                imagesc(tmp_img);
                axis equal off tight;
                subplot(222);
                imagesc(tmp_img);
                axis equal off tight;
                [tmp_r, tmp_c, ~] = find(sum(tmp_img, 3)>0);
                xlim([min(tmp_c)-10, max(tmp_c)+10]);
                ylim([min(tmp_r)-10, max(tmp_r)+10]);
                %         neuron.image(sum(Ain(:, merged_ROI{m}), 2));
                axis off;
                subplot(2,2,3:4); cla;
                tmp_C = neuron_bk.C_raw(merged_ROI{m}, :);
                %         tmp_C = bsxfun(@times, tmp_C, 1./max(tmp_C, [], 1));
                for mm=1:size(tmp_C, 1)
                    hold on;
                    plot(tmp_C(mm,:), 'color', cols(col(mm), :),  'linewidth', 2);
                end
                temp = input('keep this merge? (y(default)/n(cancel)/b(back))/e(end)   ', 's');
                if strcmpi(temp, 'n')
                    ind_after(newIDs(m)) = true;
                    ind_before(merged_ROI{m}) = true;
                    m = m+1;
                elseif strcmpi(temp, 'b')
                    m = m-1;
                elseif strcmpi(temp, 'e')
                    break;
                else
                    m = m+1;
                end
            end

            neuron.A = [neuron.A(:, ~ind_after), neuron_bk.A(:, ind_before)];
            neuron.C = [neuron.C(~ind_after, :); neuron_bk.C(ind_before, :)];
            neuron.C_raw = [neuron.C_raw(~ind_after, :); neuron_bk.C_raw(ind_before, :)];
            neuron.S = [neuron.S(~ind_after, :); neuron_bk.S(ind_before, :)];
            neuron.P.kernel_pars = [neuron.P.kernel_pars(~ind_after, :); neuron_bk.P.kernel_pars(ind_before, :)];
            close;
        end

        % view neurons
        if view_neurons
            neuron.viewNeurons([], neuron.C_raw);
        end


        %spatial
        neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
        neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
        neuron.compactSpatial();

%         cnmfe_merge_neighbors;
        neuron_bk = neuron.copy();
        if ~exist('center_method', 'var')
            center_method = 'max';
        end;
        [merged_ROI, newIDs] = neuron.MergeNeighbors(dmin, center_method);
        cols = [0, 0, 1; 0, 1, 0; 0, 1, 1; 1, 0, 0; 1, 0, 1; 1, 1, 0; 1, 1, 1];
        if display_merge && ~isempty(merged_ROI)
            figure('position', [1,1, 1200, 600]);
            ind_before = false(size(neuron_bk.A, 2), 1);
            ind_after = false(size(neuron.A, 2), 1);
            m = 1;
            while m<=length(merged_ROI)
                subplot(221);
                [tmp_img, col, ~] = neuron_bk.overlapA(merged_ROI{m});
                imagesc(tmp_img);
                axis equal off tight;
                subplot(222);
                imagesc(tmp_img);
                axis equal off tight;
                [tmp_r, tmp_c, ~] = find(sum(tmp_img, 3)>0);
                xlim([min(tmp_c)-10, max(tmp_c)+10]);
                ylim([min(tmp_r)-10, max(tmp_r)+10]);
                %         neuron.image(sum(Ain(:, merged_ROI{m}), 2));
                axis off;
                subplot(2,2,3:4); cla;
                tmp_C = neuron_bk.C_raw(merged_ROI{m}, :);
                %         tmp_C = bsxfun(@times, tmp_C, 1./max(tmp_C, [], 1));
                for mm=1:size(tmp_C, 1)
                    hold on;
                    plot(tmp_C(mm,:), 'color', cols(col(mm), :),  'linewidth', 2);
                end
                temp = input('keep this merge? (y(default)/n(cancel)/b(back))/e(end)   ', 's');
                if strcmpi(temp, 'n')
                    ind_after(newIDs(m)) = true;
                    ind_before(merged_ROI{m}) = true;
                    m = m+1;
                elseif strcmpi(temp, 'b')
                    m = m-1;
                elseif strcmpi(temp, 'e')
                    break;
                else
                    m = m+1;
                end
            end

            neuron.A = [neuron.A(:, ~ind_after), neuron_bk.A(:, ind_before)];
            neuron.C = [neuron.C(~ind_after, :); neuron_bk.C(ind_before, :)];
            neuron.C_raw = [neuron.C_raw(~ind_after, :); neuron_bk.C_raw(ind_before, :)];
            neuron.S = [neuron.S(~ind_after, :); neuron_bk.S(ind_before, :)];
            neuron.P.kernel_pars = [neuron.P.kernel_pars(~ind_after, :); neuron_bk.P.kernel_pars(ind_before, :)];
            close;
        end

        % view neurons
        if view_neurons
            neuron.viewNeurons([], neuron.C_raw);
        end

    end
    fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);

    b0 = mean(Y,2)-neuron.A*mean(neuron.C,2);
    Ybg = bsxfun(@plus, Ybg, b0-mean(Ybg, 2));
    neuron.orderROIs('snr');

    % first save
    results = neuron.obj2struct();
    save( sprintf('%s%s%s_results.mat',BasePathName,filesep,AlignedStackName), ...
        'results', 'SessionFrames', 'BasePathName', 'SubDirNames' );

    %% display neurons
    dir_neurons = sprintf('%s%s%s_neurons%s', BasePathName, filesep, AlignedStackName, filesep);
    neuron.save_neurons(dir_neurons);


    %% display contours of the neurons
    figure;
    Cnn = correlation_image(neuron.reshape(Ysignal(:, 1:5:end), 2), 4);
    neuron.show_contours(0.6);
    colormap gray;
    axis equal; axis off;
    title('contours of estimated neurons');

    % plot contours with IDs
    % [Cn, pnr] = neuron.correlation_pnr(Y(:, round(linspace(1, T, min(T, 1000)))));
    figure;
    Cn = imresize(Cn, [meta.yRes, meta.xRes]);
    neuron.show_contours(0.6);
    colormap gray;
    title('contours of estimated neurons');
    drawnow;

    fprintf('\n---- Done with source extraction ----\n\n');
    fprintf('Datadir: %s\n', pwd);
    fprintf('Date/time: %s\n', datestr(clock));

    %% check spatial and temporal components by playing movies
    save_avi = false;
    avi_name = 'play_movie.avi';
    neuron.Cn = Cn;
%     neuron.runMovie(Ysignal, [0, 150], save_avi, avi_name);


    %% save video
    kt = 1;     % play one frame in every kt frames
    save_avi = true;
    center_ac = median(max(neuron.A,[],1)'.*max(neuron.C,[],2)); % the denoised video are mapped to [0, 2*center_ac] of the colormap


%     cnmfe_save_video;
    if ~exist('t_begin', 'var')
        t_begin = 1;
    end
    if ~exist('t_end', 'var')
        t_end = size(neuron.C, 2);
    end
    if ~exist('kt', 'var')
        kt = 1;
    end

    %% data preparation
    Y = neuron.reshape(Y, 2);
    Yac = neuron.reshape(neuron.A*neuron.C, 2);
    Ybg = neuron.reshape(Ybg, 2);
    Ysignal = neuron.reshape(Ysignal, 2);


    if ~exist('center_ac', 'var')
        center_ac = median(max(neuron.A,[],1)'.*max(neuron.C,[],2));
    end
    range_res = [-1,1]*center_ac;
    if ~exist('range_ac', 'var')
        range_ac = center_ac*1.01+range_res;
    end
    if ~exist('range_Y', 'var')
        if ~exist('multi_factor', 'var')
            temp = quantile(Y(randi(numel(Y), 10000,1)), [0.01, 0.98]);
            multi_factor = ceil(diff(temp)/diff(range_ac));
    %     else
    %         temp = quantile(Y(randi(numel(Y), 10000,1)), 0.01);
        else
            temp = quantile(Y(randi(numel(Y), 10000,1)), 0.01);
        end
        center_Y = temp(1) + multi_factor*center_ac;
        range_Y = center_Y + range_res*multi_factor;
    end

       %% add pseudo color to denoised signals
    [K, T]=size(neuron.C);
    % draw random color for each neuron
    % tmp = mod((1:K)', 6)+1;
    Y_mixed = zeros(neuron.options.d1*neuron.options.d2, T, 3);
    temp = prism;
    % temp = bsxfun(@times, temp, 1./sum(temp,2));
    col = temp(randi(64, K,1), :);
    for m=1:3
        Y_mixed(:, :, m) = neuron.A* (diag(col(:,m))*neuron.C);
    end
    Y_mixed = uint16(Y_mixed/(1*center_ac)*65536);

    %% create avi file
    if save_avi
        if ~exist('avi_filename', 'var')
            avi_filename =[BasePathName, filesep, AlignedStackName];
        end
        avi_file = VideoWriter(avi_filename);
        if ~isnan(neuron.Fs)
            avi_file.FrameRate= 2*neuron.Fs/kt;
        end
        avi_file.open();
    end


  %% play and save
    fig = figure('position', [50,100, 600, 400]);
    fig.MenuBar = 'none';
    fig.NumberTitle = 'off';
    fig.Name = 'Just update data on existing axes...';
    tic
    for m=t_begin:kt:t_end
        if m == t_begin
            % axes(ax_y); cla;
            ax_y =   axes('position', [0.015, 0.51, 0.3, 0.42]);
            im_y = imagesc(Ybg(:, :,m)+Ysignal(:, :, m), range_Y);
            %     set(gca, 'children', flipud(get(gca, 'children')));
            title('Raw data');
            axis equal off tight;

            % axes(ax_bg); cla;
            ax_bg = axes('position', [0.015, 0.01, 0.3, 0.42]);
            im_bg = imagesc(Ybg(:, :, m),range_Y);
            %     set(gca, 'children', flipud(get(gca, 'children')));
            axis equal off tight;
            title('Background');

            % axes(ax_signal); cla;
            ax_signal = axes('position', [0.345, 0.51, 0.3, 0.42]);
            im_signal = imagesc(Ysignal(:, :, m), range_ac); hold on;
            %     set(gca, 'children', flipud(get(gca, 'children')));
            title(sprintf('(Raw-BG) X %d', multi_factor));
            axis equal off tight;

            % axes(ax_denoised); cla;
            ax_denoised = axes('position', [0.345, 0.01, 0.3, 0.42]);
            im_denoised = imagesc(Yac(:, :, m), range_ac);
            %     imagesc(Ybg(:, :, m), [-50, 50]);
            title(sprintf('Denoised X %d', multi_factor));
            axis equal off tight;

            % axes(ax_res); cla;
            ax_res = axes('position', [0.675, 0.51, 0.3, 0.42]);
            im_res = imagesc(Ysignal(:, :, m)-Yac(:, :, m), range_res);
            %     set(gca, 'children', flipud(get(gca, 'children')));
            title(sprintf('Residual X %d', multi_factor));
            axis equal off tight;
            %         subplot(4,6, [5,6,11,12]+12);

            % axes(ax_mix); cla;
            ax_mix = axes('position', [0.675, 0.01, 0.3, 0.42]);
            im_mix = imagesc(neuron.reshape(Y_mixed(:, m,:),2));  hold on;
            title('Demixed');
            tx_mix = text(1, 10, sprintf('Time: %.2f second', m/neuron.Fs), 'color', 'w', 'fontweight', 'bold');

            axis equal tight off;
            %     box on; set(gca, 'xtick', []);
            %     set(gca, 'ytick', []);
        else
            im_y.CData = Ybg(:, :,m)+Ysignal(:, :, m);

            im_bg.CData = Ybg(:, :, m);

            im_signal.CData = Ysignal(:, :, m);

            im_denoised.CData = Yac(:, :, m);

            im_res.CData = Ysignal(:, :, m)-Yac(:, :, m);

            im_mix.CData = neuron.reshape(Y_mixed(:, m,:),2);
            tx_mix.String = sprintf('Time: %.2f second', m/neuron.Fs);

            axis equal tight off;
            %     box on; set(gca, 'xtick', []);
            %     set(gca, 'ytick', []);
        end

        % drawnow();
        if save_avi
            temp = getframe(gcf);
            %temp = imresize(temp.cdata, [400, 600]);
            avi_file.writeVideo(temp);
        end
    end

    if save_avi
        avi_file.close();
    end
    toc


    %% save results
    Parameters = struct;

    Parameters.gSigma_align =  gSigma_align;
    Parameters.gSize_align  =  gSize_align ;
    Parameters.gSigma_cnmfe =  gSigma_cnmfe;
    Parameters.gSize_cnmfe  =  gSize_cnmfe ;
    Parameters.min_corr_cnmfe =  min_corr_cnmfe;
    Parameters.min_pnr_cnmfe =  min_pnr_cnmfe;
    Parameters.edge_margin_cnmfe =  edge_margin_cnmfe;
    Parameters.merge_thr = merge_thr;

    results = neuron.obj2struct();
    save( sprintf('%s%s%s_results.mat',BasePathName,filesep,AlignedStackName), ...
        'results', 'SessionFrames', 'BasePathName', 'SubDirNames','Parameters' );
    fprintf('Saved results as: %s%s%s_results.mat\n',BasePathName,filesep,AlignedStackName);

    fprintf('\n---- Done ----\n\n');
    fprintf('Datadir: %s\n', pwd);
    fprintf('Date/time: %s\n', datestr(clock));
end
