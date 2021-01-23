classdef (Abstract) rlAbstractSolver
    % RLABSTRACTSOLVER   Abstract RL solver for framework 2.0 dlnetwork.
    
    %   Copyright 2019 The MathWorks, Inc.
    
    %   This utility will be replaced with DLT FW2.0 external feature. 
    
    methods (Abstract)
        % LearnableParameters = calculateUpdate(this,Learnables,Gradients,LearnRate)
        % updates the parameters in Learnables. Input gradients is a
        % table of dlnetwork Learnables after regularization and
        % threshold. Input LearnRate specifies the learning rate to use
        % for calculating the update step. The size of Gradients must be 
        % similar to LearnableParameters.
        [this, LearnableParameters] = calculateUpdate(this,Gradients,LearnRate);
        
        % REVISIT: method to get statistic to resume training
    end
    
    methods (Static)
        function UseGPU = useGPU(Statistics)
            % Determine whether the solver uses GPU by looking at the
            % statistics data type.
            
            if isnumeric(Statistics)
                UseGPU = isa(Statistics,'gpuArray');
            elseif iscell(Statistics)
                UseGPU = any(cellfun(@(x) isa(x,'gpuArray'), Statistics));
            elseif istable(Statistics)
                UseGPU = any(cellfun(@(x) isa(x,'gpuArray'), Statistics.Learnables.Value));
            else
                error(message('RL Solver does not support this learnable parameters format.'))
            end
        end
    end
end
