function Gradients = regularizeGradient(Gradients,Learnables,L2RegularizationFactor)
%REGULARIZEGRADIENT performs L2 regularization of gradients Gradients wrt
%current learnable parameters Learnables and L2RegularizationFactor.
%Gradients and Learnables are cell array or tables (dlnetwork's Learnables)

%   Copyright 2019 The MathWorks, Inc.

%   This utility will be replaced with DLT FW2.0 external feature. Borrow
%   from nnet.internal.cnn.regularizer.RegularizerL2

Gradients = dlupdate( @(grad,learnables) grad + L2RegularizationFactor.*learnables, ...
                    Gradients, Learnables );
end

