pf = robotics.ParticleFilter;

initialize(pf, 5000, [initialPose(1:3)', 0, 0, 0], eye(6), 'CircularVariables',[0 0 1 0 0 0]);
pf.StateEstimationMethod = 'mean';
pf.ResamplingMethod = 'systematic';

% StateTransitionFcn defines how particles evolve without measurement
pf.StateTransitionFcn = @exampleHelperCarBotStateTransition;

% MeasurementLikelihoodFcn defines how measurement affect the our estimation
pf.MeasurementLikelihoodFcn = @exampleHelperCarBotMeasurementLikelihood;

% Last best estimation for x, y and theta
lastBestGuess = [0 0 0];

% Run loop at 20 Hz for 20 seconds
% Use fixed-rate support
r = robotics.Rate(1/dt);
% Reset the fixed-rate object
reset(r);

% Reset simulation time
simulationTime = 0;

while simulationTime < 20 % if time is not up

    % Generate motion command that is to be sent to the robot
    % NOTE there will be some discrepancy between the commanded motion and the
    % motion actually executed by the robot.
    uCmd(1) = 0.7*abs(sin(simulationTime)) + 0.1;  % linear velocity
    uCmd(2) = 0.08*cos(simulationTime);            % angular velocity

    drive(carbot, uCmd);

    % Predict the carbot pose based on the motion model
    [statePred, covPred] = predict(pf, dt, uCmd);

    % Get GPS reading
    measurement = exampleHelperCarBotGetGPSReading(carbot);

    % If measurement is available, then call correct, otherwise just use
    % predicted result
    if ~isempty(measurement)
        [stateCorrected, covCorrected] = correct(pf, measurement');
    else
        stateCorrected = statePred;
        covCorrected = covPred;
    end

    lastBestGuess = stateCorrected(1:3);

    % Update plot
    if ~isempty(get(groot,'CurrentFigure')) % if figure is not prematurely killed
        updatePlot(carbot, pf, lastBestGuess, simulationTime);
    else
        break
    end

    waitfor(r);

    % Update simulation time
    simulationTime = simulationTime + dt;
end