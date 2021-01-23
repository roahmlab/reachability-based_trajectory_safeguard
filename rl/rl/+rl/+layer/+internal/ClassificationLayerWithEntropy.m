classdef ClassificationLayerWithEntropy < nnet.layer.ClassificationLayer
% CLASSIFICATIONLAYERWITHENTROPY cte loss with entropy loss to encourages
% exploration for policy gradients methods

% Copyright 2018 The MathWorks, Inc.

% Examples:
% cteEntropy = rl.layer.internal.ClassificationLayerWithEntropy
% validInputSize = [1 1 2];
% checkLayer(cteEntropy,validInputSize,'ObservationDimension',4);
    
%   Copyright 2018 The MathWorks, Inc.
    
    properties
        EntropyLossWeight
    end
    
    methods
        function this = ClassificationLayerWithEntropy(varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser,'EntropyLossWeight',0);
            addParameter(parser,'Name',"cteWithEntropy",...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','Name'));
            addParameter(parser,'Description',getString(message('rl:general:ClassificationWithEntropyLayerDescription')),...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','Description'));
            parse(parser,varargin{:});
            r = parser.Results;
            
            this.Type = "ClassificationWithEntropyLayer";
            this.Name = r.Name;
            this.Description = r.Description;
            this.EntropyLossWeight = r.EntropyLossWeight;
        end
        
        function this = set.EntropyLossWeight(this,val)
            validateattributes(val,{'numeric'},{'scalar','finite','real','nonnegative','<=',1},'','EntropyLossWeight');
            this.EntropyLossWeight = val;
        end
        
        function loss = forwardLoss(this, Y, T)
            % Return the sum of cross entropy loss and scaled entropy losss
            % between the predictions Y and the training targets T.
            %
            % Inputs:
            %   Y: Predictions made by network
            %   T: Training targets
            % Output:
            %   loss - Loss between Y and T
            
            loss = rlcte(Y, T, this.EntropyLossWeight);
        end
    end
end

function Loss = rlcte(NetworkOutput,Target,EntropyLossWeight)
% Cross entropy loss with weighted entropy loss.
% Entropy loss encourages policy exploration
%   NetworkOutput: dlarray
%   Target: dlarray
%   EntropyLossWeight: scalar

% cte loss
ObsDim = finddim(NetworkOutput,'B');
if isempty(ObsDim)
    % Handle unlabelled dlarray from Custom Layer with no backward method
    % Request automatically labelling g2011546
    % REVISIT labels when support RNN or custom layer automatically correct labels
    CrossEntropyLoss = crossentropy(NetworkOutput, Target, 'DataFormat', 'SSCB');
else
    CrossEntropyLoss = crossentropy(NetworkOutput, Target);
end

% entropy loss: braching to improve performance
if EntropyLossWeight
    EntropyLoss = entropy(NetworkOutput, EntropyLossWeight);
else
    EntropyLoss = 0;
end

% total loss 
Loss = CrossEntropyLoss + EntropyLoss;
end

function EntropyLoss = entropy(NetworkOutput,EntropyLossWeight)
% Entropy Loss to encourage policy exploration
% Entropy = -sum(Prob(x).*log(Prob(x)))
% EntropyLoss = EntropyLossWeight * sum(Entropy)/NumObs
%   NetworkOutput: dlarray
%   EntropyLossWeight: scalar

ObsDim = finddim(NetworkOutput,'B');
if isempty(ObsDim)
    % Handle unlabelled dlarray from Custom Layer with no backward method
    NumObs = size(NetworkOutput,4);
else
    NumObs = size(NetworkOutput, ObsDim);
end

EntropyLoss = EntropyLossWeight * ...
    sum(NetworkOutput .* log(rl.internal.dataTransformation.boundAwayFromZero(NetworkOutput)),'all')./NumObs;
end