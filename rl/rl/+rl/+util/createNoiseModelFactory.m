function NoiseModel = createNoiseModelFactory(ActionSize,NoiseOpts,varargin)
% CREATENOISEMODELFACTORY Create a noise object given action size, noise 
% options and optionally sample time.

% Copyright 2019 The MathWorks, Inc.

switch class(NoiseOpts)
    case 'rl.option.OrnsteinUhlenbeckActionNoise'
        NoiseModel = rl.util.OrnsteinUhlenbeckActionNoise(ActionSize,NoiseOpts,varargin{:});
    case 'rl.option.GaussianActionNoise'
        NoiseModel = rl.util.GaussianActionNoise(ActionSize,NoiseOpts);
    otherwise
        % this should never be hit
        error(message('rl:agent:NoiseModelInvalidOptions'));
end

end