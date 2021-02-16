function [Observed] = pActive_pub(CalciumData, ObservedPermuted)
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
% ObservedPermuted: set to 1 to use observed data, 2 to create permuted data
%
% Output variable:
% Observed: Struct with pActive values per neurons per trial, ensemble size, and ensemble stability.


%% Settings
    savefile = 0;           % 1 to save file
    
    SubS_length = 50;       % subsample 5 seconds
    min_length = 100;       % trials of at least 10 seconds (10 Hz recording); shorter trials are not analyzed
    Shuf_num = 100;         % 100 resamples
    CutOff = 95;            % cut off of Pactive (Xth percentile)
    
    
%% shuffle traces

    if ObservedPermuted == 2
        NeuronDat = CalciumData.ImagingData.S;
        for nrn = 1: size(NeuronDat,1)
            ND = NeuronDat(nrn,:);
            NeuronDat(nrn,:) = ND(randperm(length(ND)));
        end    
        CalciumData.ImagingData.S = NeuronDat;
    end

%% find active neurons per session
    
    % set up structure
    Observed = struct;
    Observed(1).Trial = 'ET1'; Observed(2).Trial = 'ET2' ; Observed(3).Trial = 'ET3';
    Observed(4).Trial = 'RT1'; Observed(5).Trial = 'RT2' ; Observed(6).Trial = 'RT3';

    tic
for i = 1 : 6                                           % for each trial  
    
    disp(strcat('Trial...',num2str(i)));
    
    if CalciumData.DailySetting.TrialStructure(i).TrialDuration > min_length   
        
            % get random subsamples from baseline
            BaseIDX = [CalciumData.DailySetting.TrialStructure(i).Base(1)...
                :(CalciumData.DailySetting.TrialStructure(i).Base(2) - SubS_length)];
            BaseStart = randi(length(BaseIDX),Shuf_num,1);
            
            % get random subsamples from trial
            TrialIDX = CalciumData.DailySetting.TrialStructure(i).TrialDuration - SubS_length;
            if TrialIDX >= Shuf_num
                TrStart = randperm(TrialIDX); 
                TrStart = TrStart(1:Shuf_num);
            else
                TrStart = randi(TrialIDX,Shuf_num,1);
            end
   
        fprintf('Neuron:%3.0f',0) % neuron counter 
        
        for nrn = 1 : size(CalciumData.ImagingData.S,1)       % for each neuron
                fprintf('\b\b\b%3.0f', nrn)
                
                Neuron(nrn).ID = nrn;
                
                % create permutation matrix to compare activity to
                if i == 1
                    y_dat = [CalciumData.ImagingData.S(nrn,:)];
                    y_dat_permuted = zeros(1000, length(y_dat));
                    for x = 1:1000
                        y_dat_permuted(x,:) = y_dat(randperm(length(y_dat)));
                    end
                    PermutationMatrix(nrn).Neuron = y_dat_permuted;
                end
                
                for tr = 1 : Shuf_num                         % for each subsample
                Neuron(nrn).Trial(tr).Tr = tr;   
               
                % mean activity in baseline subsample
                Neuron(nrn).Trial(tr).BaseFrames = ...
                    [BaseIDX(BaseStart(tr)) : ((BaseIDX(BaseStart(tr)) + SubS_length - 1))];
                Neuron(nrn).Trial(tr).BaseCount = ...
                    sum(CalciumData.ImagingData.S(nrn,[Neuron(nrn).Trial(tr).BaseFrames]));
                Neuron(nrn).Trial(tr).Mean_BaseCount = ...
                    Neuron(nrn).Trial(tr).BaseCount / SubS_length;
                
                % mean activity in trial subsample
                Neuron(nrn).Trial(tr).Frames = CalciumData.DailySetting.TrialStructure(i).Trial(TrStart(tr):...
                    (TrStart(tr) + SubS_length - 1)); 
                Neuron(nrn).Trial(tr).TrialCount  = ...
                    sum(CalciumData.ImagingData.S(nrn,[Neuron(nrn).Trial(tr).Frames]));
                Neuron(nrn).Trial(tr).Mean_TrialCount  = Neuron(nrn).Trial(tr).TrialCount / SubS_length;
                
                % mean activity difference between baselin and subsample
                Neuron(nrn).Trial(tr).DifferenceScore = Neuron(nrn).Trial(tr).Mean_TrialCount - ...
                    Neuron(nrn).Trial(tr).Mean_BaseCount;
                
                % same as above, for permuted data
                y_dat_permuted = PermutationMatrix(nrn).Neuron;
                
                Neuron(nrn).Trial(tr).BasePermute = sum(y_dat_permuted...
                    (:,Neuron(nrn).Trial(tr).BaseFrames),2) / SubS_length;
                Neuron(nrn).Trial(tr).TrialPermute = sum(y_dat_permuted...
                    (:,Neuron(nrn).Trial(tr).Frames ),2) / SubS_length;
                Neuron(nrn).Trial(tr).DifferenceScorePermute = ...
                Neuron(nrn).Trial(tr).TrialPermute - Neuron(nrn).Trial(tr).BasePermute;

                % compare observed to permuted
                PermutationVector = sort(Neuron(nrn).Trial(tr).DifferenceScorePermute);
                Cut = prctile(PermutationVector,CutOff);

                Neuron(nrn).Trial(tr).ActiveYesNo = 0;
                if Neuron(nrn).Trial(tr).DifferenceScore > Cut
                    Neuron(nrn).Trial(tr).ActiveYesNo = 1;
                end
                    
                  %remove Trial- and DS-permutation
                  Neuron(nrn).Trial(tr).DifferenceScorePermute = [];
                  Neuron(nrn).Trial(tr).TrialPermute = []; 
                  Neuron(nrn).Trial(tr).BasePermute = []; 
                end
                
                Neuron(nrn).P_active = sum([Neuron(nrn).Trial.ActiveYesNo])/Shuf_num; 
        end
            fprintf('\n');

    Observed(i).P_active = [Neuron.P_active];
    clear Neuron
    
    else
        [Observed(i).P_active] = deal(nan);
    end
    
end
    toc
    
      %% Ensemble size
    ensemble = nan(6,length(Observed(1).P_active));
    
    for tr = 1:6
        Observed(tr).EnsembleSize = ...
        sum(Observed(tr).P_active > 0) / length(Observed(tr).P_active);
        Observed(tr).EnsembleActivity = ...
        median(Observed(tr).P_active(find(Observed(tr).P_active > 0)));
   end
    
    %% Ensemble stability matrix

    for tr = 1:6; ensemble(tr,:) = Observed(tr).P_active; end
    Observed(1).EnsembleStabilityMatrix = corrcoef(ensemble'); 

%% saving file 

    if savefile == 1
        if ObservedPermuted == 1
        filename = strcat(num2str(CalciumData.DailySetting.Date),'_', ...
            CalciumData.Performance.Animal,'_Observed.mat');
        save(filename,'Observed');
        elseif ObservedPermuted == 2
            Permuted = Observed;
               filename = strcat(num2str(CalciumData.DailySetting.Date),'_', ...
            CalciumData.Performance.Animal,'_Permuted.mat');
        save(filename,'Permuted');
        end
    end  

end

