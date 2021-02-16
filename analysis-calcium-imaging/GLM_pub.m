function [GLMfits] = GLM_pub(CalciumData, Shuffle)
%% READ ME
% Created in Matlab version 2016b
% Last updated 25 January 2021
%
% Used for analysis in the paper "Glas, A., Hübener, M., Bonhoeffer, T., & Goltstein, P. M. (2021).
% Spaced training enhances memory and prefrontal ensemble stability in mice"
%
% Input variables:
% Calcium data: Struct with behavioral and calcium imaging data from one mouse on one experimental session
% Calcium data can be found on the repository https://gin.g-node.org/pgoltstein/everydaymemory-spacingeffect-mpfc
% Shuffle: 1 for shuffled control, 0 for observed data
%
% Output variables:
% GLMfits: fitted GLM
%      
% Intermediate steps of the code:
% 1) creates design matrix of regressors
% 2) bins data
% 3) convolves categorical predictors with evenly spaced gaussian kernels
% 4) splits data in train and test set
% 5) fits Bernouilli GLM to training set, predicts spiking of test set 

    
%% settings
%   Shuffle = 0;   

    save_file = 0;  % 1 to save file
    Plot = 0;       % 1 to plot data
    
    DatType = 1;    % 1 for spikes, 2 for calcium data, 
    Binarize = 0;   % 0 to use real spike counts, 1 to use binarized data
    Binning = 5;    % x/10 sec bins; 5 is 500 ms
    
    Kernels = 5;            % default: 5   // spanning 2 seconds before and after 
    KernelSpacing = 5;      % default: 5   // number of frames between kernels
    sigma = 2.5;            % default: 2.5 // SD of gaussian filter; based on full width at half maximum of 500ms
    sz = 21;                % default: 21  // length of gaussFilter vector
    
    Train = 0.7; % fraction of data used for training GLM
      
    warning('off')
    
    
%% load imaging data
    if DatType == 1; ImgDat = CalciumData.ImagingData.S;
        elseif DatType ==2;  ImgDat = CalciumData.ImagingData.C;
    end
    
    % shuffle data    
    if Shuffle == 1
        ImgDat2 = nan(size(ImgDat,1),size(ImgDat,2));
        for nrn = 1:size(ImgDat,1); ND = ImgDat(nrn,:); ImgDat2(nrn,:) = ND(randperm(length(ND))); end
        ImgDat = ImgDat2; clear ImgDat2
    end
    
    % binarize data
    if Binarize == 1
        ImgDat(ImgDat > 0) = 1;
        Distribution = 'binomial';
    else
        Distribution = 'poisson';
    end
    
%%  load behavioral data
    % get arm location, distance from center
    arm = []; X = []; Y = []; meanspeed = [];
    for x = 1 : length(CalciumData.BehavioralData)
        arm = [arm, CalciumData.BehavioralData(x).mousepos.arm];
        X = [X, CalciumData.BehavioralData(1).Positions(x).x];
        Y = [Y, CalciumData.BehavioralData(1).Positions(x).y];
        meanspeed = [meanspeed, CalciumData.BehavioralData(1).Positions(x).meanspeed'];
    end
    
    %find center location and distance from center
    Center  = [(max(X) - min(X)) / 2 + min(X), (max(Y) - min(Y)) / 2 + min(Y)];
    center_distance = nan(length(X),1);
    for x = 1:length(X)
    if ~isnan(X(x)); center_distance(x) = sqrt((X(x) - Center(1)).^2 + (Y(x) - Center(2)).^2); end
    end
    center_distance = movmean(center_distance,10); 
    
    %normalize distance to max
    center_distance = center_distance - min(center_distance);
    center_distance = center_distance / max(center_distance);
    center_delta = [nan; diff(center_distance)];
    
    %find periods moving towards center (ArmToCenter) and away from center (CenterToArm)
    ArmToCenter = center_delta < -0.05;
    CenterToArm = center_delta > 0.05;

    %acceleration; normalize
    Acceleration = diff(meanspeed) / max(abs(diff(meanspeed)));     
    Acceleration(isnan(Acceleration)) = 0;
        
    % add internal state; average activity all neurons
    AllNeuron = ImgDat; 
    
    
%% step 1: create design matrix
    
    %fill in regressor cell
    regressors = {'Speed','Acceleration','Reward','Approach',...
        'DigOnset','DigOffset','EntryCenter','TurnAround'};                  % regressor name
    [regressors{2,1:size(regressors,2)}] = deal(1);                          % regressor used yes (1) / no (0)
    [regressors{3,1:size(regressors,2)}]  = deal(1);                         % regressor categorical yes (1) / no (0)
    [regressors{3,[1 2]}] = deal(0);
    
    % fill in DesignMatrix
    TrNum = zeros(length(CalciumData.AuxData(1).Total.Trial),1);    
    DM = zeros(length(CalciumData.AuxData(1).Total.Trial),size(regressors,2));
    AllTrPeriods = [];
    
    for Tr = 1:6
        TrIDX = CalciumData.DailySetting.TrialStructure(Tr).StartBehav : ...
            CalciumData.DailySetting.TrialStructure(Tr).EndBehav;
        AllTrPeriods = [AllTrPeriods, TrIDX];
        TrNum(TrIDX) = Tr; 

        % 1 speed: normalized to within trial
        DM(TrIDX,1) = meanspeed(TrIDX) / max(meanspeed(TrIDX));    

        % 2 acceleration
        DM(TrIDX,2) = Acceleration(TrIDX);      
    end

    % 3 reward
    RewardFr = find(diff(CalciumData.AuxData(1).Total.Reward)==1) + 11;                  % 11 is ofset
    DM(RewardFr(:),3) = 1;
    
    % 4 approach
    App = CalciumData.AuxData(1).Total.Approach; App(App > 1) = 1; App = find(App == 1); 
    DM(App(1:6:end),4) = 1;
       
    % 5, 6, 8; dig on & off, turn-around point
    Digging = CalciumData.AuxData(1).Total.DigAll;
        
    ArmTr = arm(AllTrPeriods); 
    deltaArm = [nan, diff(ArmTr)]; deltaArm = find(deltaArm > 0); 
    
    for idx = 1:length(deltaArm) - 1
        % find period in arm
        Period = AllTrPeriods(deltaArm(idx) : deltaArm(idx+1));
        
        %find first and last dig period
        DigOn = find(Digging(Period)==1,1,'first') + AllTrPeriods(deltaArm(idx)) - 1 ; 
        if ~isempty(DigOn)
            DM(DigOn,5) = 1;    % dig onset
            DigOff = find(Digging(Period)==1,1,'last') + AllTrPeriods(deltaArm(idx)) - 1 ;
            DM(DigOff,6) = 1;  % dig offset
        end

        % find turns: diff(CenterToArm) == -1 & arm to center
        % find period animal stays in arm
        P = diff(arm(Period(1:end))); Period = Period(1:(find(P~=0,1,'first')));
        
        CA = find(diff(CenterToArm(Period)) == -1);
        AC = find(diff(ArmToCenter(Period)) == 1);
        if length(AC) == length(CA) && length((AC==CA))==1
           DM(Period(round((CA+AC)/2)),8) = 1;
        end
    end
    
    % remove frames outside of trials
    remove_idx = 1:length(DM); remove_idx(AllTrPeriods) = []; DM(remove_idx,:) = [];
    
    % 7 all arm changes; find every time animal enters into center platform
    CenterCross = find(ArmTr == 0); deltaCenter = diff(CenterCross);
    beginCross = [CenterCross(1), CenterCross(find(deltaCenter > 1) + 1)];
    DM(beginCross,7) = 1;
    
     
   % replace nans with median values
    for x = 1:size(DM,2)
        if regressors{2,x} == 1
            idx = find(isnan(DM(:,x)));
            if ~isempty(idx); DM(idx,x) = nanmedian(DM(:,x)); end
        end
    end
    
%% step 2: bin data
    if Binning > 1
        %reshape each element to x - by - bin, nanmean over rows, fill in DM-bin
        TrNum(remove_idx) = [];
        BinIDX = [];
        for Tr = 1:6
            idx = find(TrNum == Tr);
            idx_length = ceil(length(idx)/Binning);
            if (idx_length * Binning) ~= length(idx)
                idx((end+1):(idx_length * Binning)) = nan;                 % fill out vector with nans
            end

            % reshape trial indeces to n - by bin
            idx = reshape(idx,[Binning,idx_length]); idx = idx';
            BinIDX = [BinIDX; idx];
        end

        % recreate design matrix in bins
        DM_BIN = nan(length(BinIDX),size(DM,2));
        for x = 1:length(BinIDX)
            columns = BinIDX(x,:); columns(isnan(columns)) = [];
            sub_DM = DM(columns,:);
            if size(sub_DM,1) > 1
                DM_BIN(x,:) = nanmean(sub_DM);
            else
               DM_BIN(x,:) = sub_DM;
            end
        end

        % binarize Categorical regressors
        Cat = find([regressors{3,:}]==1);
        for x = 1:length(Cat)
            DM_BIN(:,Cat(x)) = ceil(DM_BIN(:,Cat(x)));
        end
    else
        DM_BIN = DM;
    end
    
%% step 3: create kernels for categorical variables
    
    % create gaussian filter
    gaus_X = linspace(-sz / 2, sz / 2, sz);
    gaussFilter = exp(-gaus_X .^ 2 / (2 * sigma ^ 2));
    gaussFilter = gaussFilter / max (gaussFilter); % normalize   
    
    % apply gaussian filter to create offset kernel IF categorical
    % initialize fields
    reg_idx = 1;
    idx = [1:KernelSpacing:Kernels * KernelSpacing];
    Shiftnum = idx(ceil(Kernels/2));
    DM_new = [];
    regressors_new = {};
    
    for x = 1 : length(regressors)
        if regressors{3,x} == 1                                                     %which regressor
            Kern_mat = zeros(length(DM_BIN) + Kernels * KernelSpacing, Kernels);    %create kernel matrix
            Kern = conv(DM_BIN(:,x),gaussFilter,'same');                            %create filter
            
            for y = 1: Kernels                                                      %create filter offsets
                Kern_mat(idx(y):(idx(y) + length(DM_BIN)-1),y) = Kern;
            end

            Kern_mat = Kern_mat(Shiftnum:length(DM_BIN)+Shiftnum-1,:);

            DM_new = [DM_new, Kern_mat];
            [regressors_new{1,reg_idx:reg_idx+Kernels-1}] = deal(regressors{1,x});
            [regressors_new{2,reg_idx:reg_idx+Kernels-1}] = deal(regressors{2,x});
            [regressors_new{3,reg_idx:reg_idx+Kernels-1}] = deal(0);
            [regressors_new{4,reg_idx:reg_idx+Kernels-1}] = deal(x);
            reg_idx = reg_idx + Kernels;
        else
            DM_new = [DM_new, DM_BIN(:,x)];
            regressors_new{1,reg_idx} = regressors{1,x};
            regressors_new{2,reg_idx} = regressors{2,x};
            regressors_new{3,reg_idx} = regressors{3,x};
            regressors_new{4,reg_idx} = x;
            reg_idx = reg_idx +1;
        end
    end
    
    
    % fill in struct to save
    GLMfits.DesignMatrix.Full = DM_new;
    GLMfits.Regressors.Original = regressors;
    GLMfits.Regressors.Full = regressors_new;   
    GLMfits.Settings.DataType = DatType;
    GLMfits.Settings.Binning = Binning;
    GLMfits.Settings.Kernels = Kernels;
    GLMfits.Settings.KernelSpacing = KernelSpacing;
    GLMfits.Settings.Sigma = sigma;
    GLMfits.Settings.GaussianSteps = sz;
    
  
    % initialize GLM struct 
        GLM_idx = 1;
        GLMFIT = struct;
        CatVar = logical([regressors_new{3,find([regressors_new{2,:}] == 1)}]);
              
%% step 4: create indeces of the test and training set
        TrainSet = randsample(size(DM_new,1), round(Train * size(DM_new,1)),'false');
        TrainSet = sort(TrainSet)';
        TestSet = 1:size(DM_new,1); TestSet(TrainSet) = [];
      
%% step 5: fit GLM
       
        tic
        nrn = 1; fprintf('Neuron:%3.0f',0) % neuron counter 
        for nrn = 1 : size(ImgDat,1)  % loop over neurons   

             fprintf('\b\b\b%3.0f', nrn) % display neuron counter
                
            
             NeurDat = ImgDat(nrn,AllTrPeriods);
             
             % binning
             if Binning > 1
                 NeurDat2 = nan(length(BinIDX),1);
                 for x = 1:length(BinIDX)
                     columns = BinIDX(x,:); columns(isnan(columns)) = []; 
                     NeurDat2(x) = nanmean(NeurDat(columns));
                 end
                 if Binarize == 1
                     NeurDat2(NeurDat2 > 0) = 1; 
                 end
                 NeurDat = NeurDat2; clear NeurDat2
             end
                         
             % copy Design Matrix for manipulation
             DM_new2 = DM_new;

             %max normalize design matrix
             for col = 1 : size(DM_new2,2)
                 DM_new2(:,col) =  mat2gray(DM_new2(:,col));
             end

            % lasso regularization: Construct a regularized binomial regression 
            % using 10 Lambda values and 10-fold cross validation
            [B,FitInfo] = lasso(DM_new2(TrainSet,:),NeurDat(TrainSet),'NumLambda',10,'CV',10);
            
            % Find the number of nonzero model coefficients at the Lambda value 
            % with minimum deviance plus one standard deviation point
            indx = FitInfo.IndexMinMSE; B0 = B(:,indx);
            KeepRegressors = find(B(:,indx));
            
             % variance inflation factor
             R_lasso = corrcoef(DM_new2(:,KeepRegressors)); 
             VIF_lasso = diag(inv(R_lasso))';
             % remove regressors with VIF > 4
             VIF_reject = find(abs(VIF_lasso) > 4); 
             B0(KeepRegressors(VIF_reject)) = 0;
             KeepRegressors(VIF_reject) = [];
             
            % fill in struct; lasso regression
            GLMFIT(GLM_idx).Neuron = nrn;
            GLMFIT(GLM_idx).LassoRegression.DesignMatrix = DM_new2;
            GLMFIT(GLM_idx).LassoRegression.Lambda = FitInfo.LambdaMinMSE;
            GLMFIT(GLM_idx).LassoRegression.Beta = B0;
            GLMFIT(GLM_idx).LassoRegression.VIF = VIF_lasso;        
                        
            if isempty(KeepRegressors) 
                [GLMFIT(GLM_idx).FullModel.Beta, GLMFIT(GLM_idx).FullModel.P_val] = deal(nan(1,size(regressors,2)));            
                GLM_idx = GLM_idx + 1;
                continue
            end
            
            % run second model with kernels to get modulated weights
            % create 2nd design matrix with kernels
            Predictors = [];
            RegressorWeights = B0; 
            
            DM = zeros(length(DM_new2),max([regressors_new{4,:}]));
            for reg = 1:max([regressors_new{4,:}])
                idx = find([regressors_new{4,:}] == reg);
                % if single regressor: use yes / no
                if length(idx) == 1 && RegressorWeights(idx)~=0
                    DM(:,reg) = DM_new2(:,idx);             
                    Predictors = [Predictors,reg];
                end
                % if multiple regressors and used: create kernel
                if length(idx) > 1 && sum(RegressorWeights(idx))~=0
                    vector = DM_new2(:,idx) .* RegressorWeights(idx)';
                    vector = sum(vector,2);
                    DM(:,reg) = vector;             
                    Predictors = [Predictors,reg];   
                end
            end
            
             % variance inflation factor
             R_CompactModel = corrcoef(DM(:,Predictors)); 
             VIF_CompactModel = diag(inv(R_CompactModel))';
             % remove regressors with VIF > 4
             VIF_reject = find(abs(VIF_CompactModel) > 4); 
             Predictors(VIF_reject) = [];
            
            %fit model with kernels
            Model2 = fitglm(DM(TrainSet,:),NeurDat(TrainSet),'linear','Distribution',Distribution,...
            'PredictorVars', Predictors); 

            % predict model
            Prediction = predict(Model2,DM(TestSet,:));
            
            % create summary stats;
            Coeff = nan(1,size(DM,2)); Coeff(Predictors) = Model2.Coefficients{[2:1+length(Predictors)],1};
            P_val = nan(1,size(DM,2)); P_val(Predictors) = Model2.Coefficients{[2:1+length(Predictors)],4};
            P_val(Predictors) = P_val(Predictors) < (0.05 / length(Predictors)); % bonferroni correction
            
            
            % fill in GLM fit overall struct
             GLMFIT(GLM_idx).FullModel.FullModel = Model2; 
             GLMFIT(GLM_idx).FullModel.DesignMatrix = DM(:,Predictors);
             GLMFIT(GLM_idx).FullModel.VIF = VIF_CompactModel;   
             GLMFIT(GLM_idx).FullModel.Regressors = Predictors;
             GLMFIT(GLM_idx).FullModel.Beta = Coeff;
             GLMFIT(GLM_idx).FullModel.P_val = P_val;
             GLMFIT(GLM_idx).FullModel.Rsquared = Model2.Rsquared.Adjusted;
             GLMFIT(GLM_idx).FullModel.ResponseActual = NeurDat(TestSet);
             GLMFIT(GLM_idx).FullModel.ResponseModel = Prediction;

            GLM_idx = GLM_idx + 1;

        end
        
         fprintf('\n'); % end neuron counter
         
        GLMfits.Settings.TrainIdx = TrainSet;
        GLMfits.Settings.TestIdx = TestSet;
        GLMfits.NeuronData = GLMFIT; clear GLMFIT

        toc
    
    %% Summary stats
    Coeff = []; Pval = []; Rfit = [];
    for nrn = 1:length(GLMfits.NeuronData)
        Coeff = [Coeff; GLMfits.NeuronData(nrn).FullModel.Beta];
        Pval = [Pval; GLMfits.NeuronData(nrn).FullModel.P_val];
        if isfield(GLMfits.NeuronData(nrn).FullModel,'Rsquared')
            Rfit = [Rfit, GLMfits.NeuronData(nrn).FullModel.Rsquared]; end
    end
    GLMfits.Population.Beta = Coeff;
    GLMfits.Population.Pval = Pval;
    GLMfits.Population.Rfit = Rfit;
   
    %% Plot stats
    
    % plot DesignMatrix
    if Plot == 1; disp('PlottingData'); plotDesignMatrix(DM_new, regressors_new, save_file, Shuffle); end 

    % bargraph with 99-percent trimmed beta coefficients
    if Plot == 1; plotBarBeta(Coeff, regressors, save_file, Shuffle); end

    % heatmap with individual coefficients and significance
    if Plot == 1; plotHeatMap(Coeff, Pval, regressors_new, save_file, Shuffle); end

   % plot correlation matrix
    if Plot == 1; plotCorrelationMatrix(DM_new, regressors_new, save_file, Shuffle); end
       
           

 %% saving
 if save_file == 1
    cur_date = clock; 
    cur_date = [num2str(cur_date(1)),num2str(cur_date(2)),num2str(cur_date(3))];
    save(strcat('GLM_Shuffle_',num2str(Shuffle),'_',cur_date,'_Binarize_',num2str(Binarize),'.mat'),'GLMfits');
 end
 close all
 clc

end
    
function plotDesignMatrix(DM_new2, regressors_new, save_file, Shuffle)
    figure;imagesc(DM_new2); 
    Reg_keep = find([regressors_new{2,:}]==1);
    gcf; xlim([.5 size(Reg_keep,2)+.5]); xticks([1:1:size(Reg_keep,2)]); 
    xtick_idx = [regressors_new{2,Reg_keep}];
    xticklabels({regressors_new{1,Reg_keep}});
    xtickangle(90); colorbar; set(gca,'FontSize',8);
    box off; set(gca,'TickDir','out'); ylabel('Binned Trial')
    if save_file == 1
        cur_dir = pwd; cd('Figs'), if isempty(dir('GLM')); mkdir('GLM'); end
        cd('GLM'); saveas(gcf,strcat('DesignMatrix_',num2str(Shuffle),'.fig')); cd(cur_dir)
    end
    
end

function plotBarBeta(Coeff,regressors, save_file, Shuffle)
    figure; hold on
    bar(nanmean(Coeff));
    xlim([.5 size(Coeff,2)+.5]); xticks(1:size(Coeff,2)); xticklabels({regressors{1,:}});
    xtickangle(90); set(gca,'fontsize',6); set(gca,'TickDir','out'); xlabel('Regressor'),ylabel('Beta')
    if save_file == 1
        cur_dir = pwd; cd('Figs/GLM'); saveas(gcf,strcat('Betas_',num2str(Shuffle),'.fig')); cd(cur_dir)
    end
end

function plotHeatMap(Coeff, Pval, regressors_new, save_file, Shuffle)
    figure; hold on;
    
     % color nans
    RowCol = []; [RowCol(1,:),RowCol(2,:)] = find(isnan(Pval));
    value = min(min(Coeff)) - .1 * min(min(Coeff));
    for x = 1:size(RowCol,2)
        Coeff(RowCol(1,x), RowCol(2,x)) = value;
    end
    
    map = hsv; map(1,:) = 0;
    
    imagesc(Coeff); colormap(map); box off; axis tight
    xlim([.5 size(Coeff,2)+.5]); xticks(1:size(Coeff,2));  xticklabels({regressors_new{1,:}});
    xtickangle(90); set(gca,'fontsize',6); set(gca,'TickDir','out'); set(gca,'YDir','reverse')
    xlabel('Regressor'),ylabel('Neuron'); colorbar; 

    % plot significant values
    RowCol = []; [RowCol(1,:),RowCol(2,:)] = find(Pval ==1);
    if ~isempty(RowCol)
        for x = 1:size(RowCol,2)
            plot(RowCol(2,x),RowCol(1,x),'+k')
        end
    end
    
    if save_file == 1
        cur_dir = pwd; cd('Figs/GLM'); saveas(gcf,strcat('HeatMap_',num2str(Shuffle),'.fig')); cd(cur_dir)
    end
    
end

function plotCorrelationMatrix(DM_new, regressors_new, save_file, Shuffle)
 R= corrcoef(DM_new);
        figure; imagesc(R); colormap(hsv); box off
        for x = 1:size(DM_new,2)
            for y =1:size(DM_new,2)
                text(x-.25,y,sprintf('%.2f',R(x,y)))
            end
        end
        xlim([1 size(DM_new,2)]); xticks([1:1:size(DM_new,2)]); xticklabels({regressors_new{1,:}}); xtickangle(90)
        ylim([1 size(DM_new,2)]); yticks([1:1:size(DM_new,2)]); yticklabels({regressors_new{1,:}});
        colorbar; set(gca,'TickDir','out'); title('Correlation Matrix'); axis('tight')
        
         if save_file == 1
              cur_dir = pwd; cd('Figs/GLM'); saveas(gcf,strcat('CorrelationMatrix_',num2str(Shuffle),'.fig')); cd(cur_dir)
          end
end