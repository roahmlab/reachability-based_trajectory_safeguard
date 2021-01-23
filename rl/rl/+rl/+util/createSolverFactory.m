function Solver = createSolverFactory(Options, Model)
% CREATESOLVERFACTORY Create a solver to optimize learnable parameters from
% gradients.

% Copyright 2019 The MathWorks, Inc.

switch class(Model)
    case 'rl.representation.model.rlLayerModel'
        % rewrite DLT internal solver for Layer API to track learn rate
        % factor (support transfer learning workflow) as value object
        switch Options.Optimizer
            case "adam"
                Solver = rl.internal.optimizer.rlLayerADAMSolver(Options, Model);
            case "sgdm"
                Solver = rl.internal.optimizer.rlLayerSGDMSolver(Options, Model);
            case "rmsprop"
                Solver = rl.internal.optimizer.rlLayerRMSPropSolver(Options, Model);
        end
    otherwise
        % wrapper for dlarray's update methods (e.g. adamupdate)
        switch Options.Optimizer
            case "adam"
                Solver = rl.internal.optimizer.rlADAMSolver(Options);
            case "sgdm"
                Solver = rl.internal.optimizer.rlSGDMSolver(Options);
            case "rmsprop"
                Solver = rl.internal.optimizer.rlRMSPropSolver(Options);
        end
end
end