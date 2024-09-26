%% load a prepared experiment (after manual cd)
% subject = '';
% load(['Experiment_', subject, '.mat']);
% load(['Stimulus_', subject, '.mat']);
% load(['Results_', subject, '.mat']);

%% resume experiment
[currentLevel, linearIdx] = find(permute(Experiment.expCompleted,[3,2,1]) == 0, 1);

if isempty(currentLevel)
    error('Congratulations! This experiment has been completed.');
else
    [currentCondition, currentSession] = ind2sub([Experiment.nConditions, Experiment.nLevels], linearIdx);
end

disp(['Current session: ', num2str(currentSession)]);
disp(['Current condition: ', num2str(currentCondition)]);
disp(['Current level: ', num2str(currentLevel)]);

S = currentSession;
C = currentCondition;

%% functions
gamma_correction = @(x, bit, gamma) (x./(2^bit-1)).^(1/gamma)*(2^bit-1);

%% Psychtoolbox initialization
rng('default');
rng('shuffle');
Screen('Preference', 'SkipSyncTests', 1);
sca;
PsychDefaultSetup(2);
screenNumber = max(Screen('Screens'));
bgMeanGamma = gamma_correction(Stimulus.bgMean, Experiment.monitorBit, Experiment.monitorGamma);
[win, winRect] = Screen('OpenWindow', screenNumber, bgMeanGamma);
Experiment.monitorPx = winRect(3:4);

%% Keyboard
KbName('UnifyKeyNames');
left = KbName('LeftArrow');
right = KbName('RightArrow');

%% Tones
toneH = 600;
toneL = toneH/3;

%% Luminance sheet
luminanceSheet = ones(Experiment.monitorPx(2), Experiment.monitorPx(1)) * Stimulus.bgMean;

%% Stimulus mask
shiftedSpotCenters = Stimulus.spotCenters;
shiftedSpotCenters(:,1) = shiftedSpotCenters(:,1) + (Experiment.monitorPx(2)-Stimulus.totalLength)/2;
shiftedSpotCenters(:,2) = shiftedSpotCenters(:,2) + (Experiment.monitorPx(1)-Stimulus.totalLength)/2;

[displayY, displayX] = meshgrid(1:Experiment.monitorPx(1), 1:Experiment.monitorPx(2));
displayMask = boolean(zeros(Experiment.monitorPx(2), Experiment.monitorPx(1)));
for sp = 1:Stimulus.nLocations-1
    displayMask(sqrt((displayX-shiftedSpotCenters(sp,1)).^2+(displayY-shiftedSpotCenters(sp,2)).^2)...
        <= Stimulus.spotLength/2)= true;
end

%% Main
HideCursor;
start_time = clock;
levelCorrect = 0;
for L = currentLevel:Experiment.nLevels
    if Experiment.expCompleted(S, C, L) == 0
        %% Circular cues
        circularCue = zeros(Experiment.monitorPx(2), Experiment.monitorPx(1));
        for sp = 1:Stimulus.nLocations-1
            centerRing = round(sqrt((displayX-shiftedSpotCenters(sp,1)).^2+(displayY-shiftedSpotCenters(sp,2)).^2))...
                == floor(Stimulus.spotLength/2);
            innerRing = round(sqrt((displayX-shiftedSpotCenters(sp,1)).^2+(displayY-shiftedSpotCenters(sp,2)).^2))...
                == floor(Stimulus.spotLength/2) - 1;
            outerRing = round(sqrt((displayX-shiftedSpotCenters(sp,1)).^2+(displayY-shiftedSpotCenters(sp,2)).^2))...
                == floor(Stimulus.spotLength/2) + 1;

            if sp == C
                circularCue(centerRing | innerRing | outerRing) = - Stimulus.bgMean*0.2;
            else
                circularCue(centerRing | innerRing | outerRing) = Stimulus.bgMean*0.2;
            end
        end

        %% pre-display
        sampleStimulus = luminanceSheet+circularCue;
        target = Stimulus.targets * Stimulus.targetAmplitude(S, C, L) * Stimulus.bgMean;
        for sp = 1:C
            targetCoordX = [shiftedSpotCenters(sp,1)-Stimulus.targetLength/2+1:shiftedSpotCenters(sp,1)+Stimulus.targetLength/2];
            targetCoordY = [shiftedSpotCenters(sp,2)-Stimulus.targetLength/2+1:shiftedSpotCenters(sp,2)+Stimulus.targetLength/2];
            sampleStimulus(targetCoordX, targetCoordY) = sampleStimulus(targetCoordX, targetCoordY) + target;
        end
        
        sampleStimulus_gc = gamma_correction(sampleStimulus, Experiment.monitorBit, Experiment.monitorGamma);
        tex = Screen('MakeTexture', win, sampleStimulus_gc);
        Screen('DrawTexture', win, tex);
        Screen('Flip', win);
        WaitSecs(Experiment.responseWait);

        [clicks,x,y,whichButton,clickSecs] = GetClicks(win);
        %% show uniform luminance
        texture = round(gamma_correction(luminanceSheet, Experiment.monitorBit, Experiment.monitorGamma));
        tex = Screen('MakeTexture', win, texture);
        Screen('DrawTexture', win, tex);
        Screen('Flip', win, GetSecs() + Experiment.cueOff);

        for T = 1:Experiment.nTrials
            
            %% generate stimulus = background (+ target)
            stimulus = luminanceSheet;
            for sp = 1:Stimulus.nLocations-1
                sceneCoordX = [Stimulus.sceneCenters(1, T, S, C, L, sp)-Stimulus.spotLength/(2*Stimulus.spotResize)+1:Stimulus.sceneCenters(1, T, S, C, L, sp)+Stimulus.spotLength/(2*Stimulus.spotResize)];
                sceneCoordY = [Stimulus.sceneCenters(2, T, S, C, L, sp)-Stimulus.spotLength/(2*Stimulus.spotResize)+1:Stimulus.sceneCenters(2, T, S, C, L, sp)+Stimulus.spotLength/(2*Stimulus.spotResize)];
                spotCoordX = [shiftedSpotCenters(sp,1)-Stimulus.spotLength/2+1:shiftedSpotCenters(sp,1)+Stimulus.spotLength/2];
                spotCoordY = [shiftedSpotCenters(sp,2)-Stimulus.spotLength/2+1:shiftedSpotCenters(sp,2)+Stimulus.spotLength/2];
                stimulus(spotCoordX, spotCoordY) = stimulus(spotCoordX, spotCoordY) + imresize(Stimulus.scene(sceneCoordX, sceneCoordY)-mean2(Stimulus.scene(sceneCoordX, sceneCoordY)),Stimulus.spotResize,'nearest')...
                    / Stimulus.backgroundStd * Stimulus.bgContrast * Stimulus.bgMean;
                
                if Stimulus.tLocation(T, S, C, L) == sp
                    targetCoordX = [shiftedSpotCenters(sp,1)-Stimulus.targetLength/2+1:shiftedSpotCenters(sp,1)+Stimulus.targetLength/2];
                    targetCoordY = [shiftedSpotCenters(sp,2)-Stimulus.targetLength/2+1:shiftedSpotCenters(sp,2)+Stimulus.targetLength/2];
                    stimulus(targetCoordX, targetCoordY) = stimulus(targetCoordX,targetCoordY) + target;
                end
            end
            
            stimulus(~displayMask) = Stimulus.bgMean;

            %% clipping and gamma correction
            stimulus(stimulus < 0) = 0;
            stimulus(stimulus > 2^Experiment.monitorBit-1)  = 2^Experiment.monitorBit-1;
            stimulus_gc = round(gamma_correction(stimulus, Experiment.monitorBit, Experiment.monitorGamma));
            
            %% show circular cues
            texture = luminanceSheet + circularCue;
            texture = round(gamma_correction(texture, Experiment.monitorBit, Experiment.monitorGamma));
            tex = Screen('MakeTexture', win, texture);
            Screen('DrawTexture', win, tex);
            Screen('Flip', win);

            trialStart = GetSecs();

            %% show display circle without circular cues
            texture = luminanceSheet;
            texture = round(gamma_correction(texture, Experiment.monitorBit, Experiment.monitorGamma));
            tex = Screen('MakeTexture', win, texture);
            Screen('DrawTexture', win, tex);
            Screen('Flip', win, trialStart + Experiment.cueOff);

            %% show the stimulus
            tex = Screen('MakeTexture', win, stimulus_gc);
            Screen('DrawTexture', win, tex);
            Screen('Flip', win, trialStart + Experiment.cueOff + Experiment.stimOn);

            %% show circular cues
            texture = luminanceSheet + circularCue;
            texture = round(gamma_correction(texture, Experiment.monitorBit, Experiment.monitorGamma));
            tex = Screen('MakeTexture', win, texture);
            Screen('DrawTexture', win, tex);
            Screen('Flip', win, trialStart + Experiment.cueOff + Experiment.stimOn + Experiment.stimOff);

            %% time for mouse response
            SetMouse(Experiment.monitorPx(1)/2, Experiment.monitorPx(2)/2, win);
            responseStart = GetSecs();
            [clicks,x,y,whichButton,clickSecs] = GetClicks(screenNumber, 0, 'untilTime', responseStart+Experiment.responseWait);
            if clicks == 1
                Results.reactionTime(T, S, C, L) = clickSecs - responseStart;
                
                if whichButton == 1
                    Results.hResponse(T,S,C,L) = C;
                else
                    Results.hResponse(T,S,C,L) = 0;
                end
            else
                warning('You did not make a decision in the trial.');
                Results.hResponse(T,S,C,L) = 0;
            end
                
            %% auditory feedback
            if(Stimulus.tLocation(T, S, C, L) == Results.hResponse(T,S,C,L))
                Beeper(toneH, 1.5, .2);
                levelCorrect = levelCorrect + 1;
            else
                Beeper(toneL, 1.5, .2);
            end
        end
    
        %% save results
        Experiment.expCompleted(S, C, L) = 1;
        save(['Experiment_', Experiment.subject, '.mat'], 'Experiment');
        save(['Results_', Experiment.subject, '.mat'], 'Results');


    end
end

Stimulus_no_background = rmfield(Stimulus, 'scene');
save(['p1_vars_Results_',Experiment.subject,'.mat'], '-struct', 'Results');
save(['p1_vars_Experiments_',Experiment.subject,'.mat'], '-struct', 'Experiment');
save(['p1_vars_Stimulus_',Experiment.subject,'.mat'], '-struct', 'Stimulus_no_background');

Screen('Close', win);
ShowCursor;

fprintf('\n*****Progress summary*****\n');
disp(['Average percent correct for recent levels: ', num2str(levelCorrect / ((Experiment.nLevels-currentLevel+1)*Experiment.nTrials))]);
disp(['Session finishing: ', num2str(currentSession), '/',num2str(Experiment.nSessions)]);
disp(['Condition finished: ', num2str(currentCondition), '/',num2str(Experiment.nConditions)]);
disp(['Level finished: ', num2str(Experiment.nLevels), '/',num2str(Experiment.nLevels)]);
disp(['Used time: ', num2str(etime(clock, start_time)/60.0), '(min)']);
