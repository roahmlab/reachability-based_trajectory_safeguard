function Gradients = thresholdGradient(Gradients,GradientThreshold,ThresholdMethod)
%THRESHOLDGRADIENT clips the Gradients according to ThresholdMethod if they
%exceed the value of GradientThreshold.
%Gradients and Learnables are tables (dlnetwork's Learnables property)

%   Copyright 2019 The MathWorks, Inc.

%   This utility will be replaced with DLT FW2.0 external feature. Borrow
%   from nnet.internal.cnn.GradientThresholder

% NOTE: we want to avoid error checking at every iteration so 
% validatestring is called only on 'otherwise ' case.
switch lower(ThresholdMethod)
    case 'l2norm'
        Gradients = dlupdate( ...
            @(grad)iThresholdL2Norm(grad, GradientThreshold), ...
            Gradients);
    case 'global-l2norm'
        Gradients.Value = iThresholdGlobalL2Norm(Gradients.Value, GradientThreshold);
    case 'absolute-value'
        Gradients = dlupdate( ...
            @(grad)iThresholdAbsoluteValue(grad, GradientThreshold), ...
            Gradients);
    otherwise
        validatestring(ThresholdMethod,{'l2norm','global-l2norm','absolute-value'},'','ThresholdMethod')
end
end

%% Local functions
function gradient = iThresholdL2Norm(gradient, normThreshold)
% If the L2 norm of the gradient of a learnable parameter is larger than 
% GradientThreshold, then scale the gradient so that the L2 norm equals 
% GradientThreshold.

gradientNorm = sqrt(sum(gradient(:).^2));
if gradientNorm > normThreshold
    gradient = gradient*(normThreshold/gradientNorm);
end
end

function gradients = iThresholdGlobalL2Norm(gradients, normThreshold)
% If the global L2 norm, L, is larger than GradientThreshold, then scale 
% all gradients by a factor of GradientThreshold/L. The global L2 norm 
% considers all learnable parameters.

globalL2Norm = 0;
for i = 1:numel(gradients)
    globalL2Norm = globalL2Norm + sum(gradients{i}(:).^2);
end
globalL2Norm = sqrt(globalL2Norm);

if globalL2Norm > normThreshold
    normScale = (normThreshold/globalL2Norm);
    for i = 1:numel(gradients)
        gradients{i} = gradients{i}*normScale;
    end
end

end

function gradient = iThresholdAbsoluteValue(gradient, valueThreshold)
% If the absolute value of an individual partial derivative in the gradient
% of a learnable parameter is larger than GradientThreshold, then scale the
% partial derivative to have magnitude equal to GradientThreshold and 
% retain the sign of the partial derivative.

gradient(gradient > valueThreshold) = valueThreshold;
gradient(gradient < -valueThreshold) = -valueThreshold;
end

