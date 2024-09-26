%% refresh environment
clear;
rng('default');
rng('shuffle');

%% functions
inverse_gamma_correction = @(x, bit, gamma) (x./(2^bit-1)).^(gamma)*(2^bit-1);

%% record experiment context
Experiment.monitorBit = 8;
Experiment.monitorGamma = 2.089;
Experiment.matlabvers = version();
Experiment.title = 'Fovea neglect effect for visual search in white noise -- Part 3 (n91)';
Experiment.subject ='placeholder_t1';
Experiment.nTrials = 1;
Experiment.nSessions = 1;
Experiment.sessionLabels={'1','2','3', '4'};
Experiment.nConditions = 1;
Experiment.conditionLabels = {'white noise'};
Experiment.nLevels = 1;
Experiment.expCompleted = zeros(Experiment.nSessions, Experiment.nConditions, Experiment.nLevels);
Experiment.cueOff = .75;
Experiment.stimOn = .25;
Experiment.stimOff = .25;
Experiment.responseWait = 3;

%% generate stimulus parameters
Stimulus.ppd = 60;
Stimulus.totalLength = 1200;
Stimulus.targets = readmatrix('../target_calen_60cpd.csv');
Stimulus.targetLength = size(Stimulus.targets, 1);

Stimulus.nDirections = 6;
Stimulus.spotLength = 94;
Stimulus.spotResize = 2;
spotPrepLen = Stimulus.spotLength / Stimulus.spotResize;
Stimulus.spotDistance = 108;
Stimulus.spotCenters = [Stimulus.totalLength/2,Stimulus.totalLength/2];

Stimulus.spotCenters = find_spot_centers(Stimulus.spotCenters, Stimulus.spotCenters,...
    Stimulus.spotLength, Stimulus.spotDistance, Stimulus.totalLength, Stimulus.nDirections);
Stimulus.nLocations = 1 + size(Stimulus.spotCenters,1);

Stimulus.clipTolerance = 0.001;
Stimulus.bgMean = 2^(Experiment.monitorBit-1);
Stimulus.bgContrast = 0.2;

Stimulus.targetAmplitude = zeros(Experiment.nSessions, Experiment.nConditions, Experiment.nLevels); 

Stimulus.tLocation = zeros(Experiment.nTrials, Experiment.nSessions, Experiment.nConditions, Experiment.nLevels);
for S = 1:Experiment.nSessions
    for C = 1:Experiment.nConditions
        for L = 1:Experiment.nLevels
            Stimulus.targetAmplitude(S,C,L) = 7;
            for T = 1:Experiment.nTrials
                if rand < 0.5
                    Stimulus.tLocation(T, S, C, L) = 0;
                else
                    Stimulus.tLocation(T, S, C, L) = randi(Stimulus.nLocations-1);
                end
            end
        end
    end   
end

%% create patches
n_background = Experiment.nLevels*Experiment.nTrials;

Stimulus.fieldRatio = 8;
Stimulus.sceneLength = Stimulus.fieldRatio * Stimulus.totalLength;
Stimulus.scene = create_power_noise_field(Stimulus.sceneLength, Stimulus.sceneLength, 0, Stimulus.fieldRatio);
Stimulus.sceneCenters = zeros(2, Experiment.nTrials, Experiment.nSessions, Experiment.nConditions, Experiment.nLevels, Stimulus.nLocations-1);
Stimulus.backgroundStd = evaluate_patch_std(Stimulus.scene, Stimulus.spotLength, Stimulus.spotLength);

for S=1:Experiment.nSessions
    for C=1:Experiment.nConditions
        for L=1:Experiment.nLevels
            for T=1:Experiment.nTrials
                for sp=1:Stimulus.nLocations-1
                    while 1
                        center_x = randi(Stimulus.sceneLength-spotPrepLen)+spotPrepLen/2;
                        center_y = randi(Stimulus.sceneLength-spotPrepLen)+spotPrepLen/2;
                        coords_scene = [center_x-spotPrepLen/2+1:center_x+spotPrepLen/2];
                        candidate = imresize(Stimulus.scene(coords_scene,coords_scene), Stimulus.spotResize, 'nearest');
                        coords_target = [(Stimulus.spotLength-Stimulus.targetLength)/2+1:(Stimulus.spotLength+Stimulus.targetLength)/2];
                        candidate(coords_target,coords_target) = candidate(coords_target,coords_target) + Stimulus.targets;
                        
                        if sum(candidate(:)<-1/Stimulus.bgContrast)+sum(candidate(:)>((2^Experiment.monitorBit-1)/Stimulus.bgMean-1)/Stimulus.bgContrast) < Stimulus.clipTolerance * (Stimulus.spotLength^2)
                            break;
                        end
                        
                        disp('Overclipped.');
                    end
                    
                    Stimulus.sceneCenters(1, T, S, C, L, sp) = center_x;
                    Stimulus.sceneCenters(2, T, S, C, L, sp) = center_y;
                end
            end
        end
    end
end

%% store experiment and stimulus
save(['p3_n91_Experiment_', Experiment.subject, '.mat'], 'Experiment');
save(['p3_n91_Stimulus_', Experiment.subject, '.mat'], 'Stimulus', '-v7.3');

%% store results
Results.hResponse = zeros(Experiment.nTrials, Experiment.nSessions, Experiment.nConditions, Experiment.nLevels);
Results.reactionTime = zeros(Experiment.nTrials, Experiment.nSessions, Experiment.nConditions, Experiment.nLevels);
save(['p3_n91_Results_', Experiment.subject, '.mat'], 'Results');