classdef ClipPPOLossLayer < nnet.layer.ClassificationLayer
    % CLIPPPOLOSSLAYER clipped loss limit the change in each policy update step
    % by clipping the action probability ratio between the new and old policies.
    % Include entropy loss to encourages policy exploration.
    
    % Examples:
    % clipLoss = rl.layer.internal.ClipPPOLossLayer
    
    %   Copyright 2019 The MathWorks, Inc.
    
    properties
        % Clip factor to limit the change in each policy update step
        ClipFactor
        % Weight for entropy loss to promote policy exploration.
        EntropyLossWeight
    end
    
    methods
        function this = ClipPPOLossLayer(varargin)
            % Input parser
            parser = inputParser;
            addParameter(parser,'EntropyLossWeight',0);
            addParameter(parser,'ClipFactor',0);
            addParameter(parser,'Name',"clipPPOLoss",...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','Name'));
            addParameter(parser,'Description',getString(message('rl:general:PPOClipLossDescription')),...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','Description'));
            parse(parser,varargin{:});
            r = parser.Results;
            
            this.Type = "ClipPPOLossLayer";
            this.Name = r.Name;
            this.Description = r.Description;
            this.EntropyLossWeight = r.EntropyLossWeight;
            this.ClipFactor = r.ClipFactor;
        end
        
        function Loss = forwardLoss(this, NetworkOutput, LossVariable)
            % Return the loss wrt the network output and the loss variables
            %
            % Inputs:
            %   NetworkOuput: Predictions made by network
            %   LossVariable: Contains previous action, old action
            %                 probabilities, computed advantage
            % Output:
            %   Loss - Loss computed from NetworkOuput and LossVariable
            
            Loss = ppoClipped(NetworkOutput, LossVariable, this.ClipFactor, this.EntropyLossWeight);
        end
    end
end

function Loss = ppoClipped(NetworkOutput, LossVariable, ClipFactor, EntropyLossWeight)
% Clipped PPO with entropy loss function function for discrete action space
%   NetworkOutput: dlarray of current policy action probabilities
%   LossVariable: cell array of dlarray.
%       - LossVariable{1}: logical matrix of action index
%       - LossVariable{2}: old action probability piOld(at|st)
%       - LossVariable{3}: advantage
%   ClipFactor: scalar > 0
%   EntropyLossWeight: scalar where 0 <= EntropyLossWeight <= 1

% Extract information from input
ActionIdx = LossVariable{1};
OldActionProb = LossVariable{2};
Advantage = LossVariable{3};

% rt = pi(at|st)/piOld(at|st), avoid division by zero
Ratio = NetworkOutput(ActionIdx) ./ rl.internal.dataTransformation.boundAwayFromZero(OldActionProb(ActionIdx));

% ensure same dimension
Ratio = reshape(Ratio,[],1);
Advantage = reshape(Advantage,[],1);

% obj = rt * At
Objective = Ratio .* Advantage;
ObjectiveClip = max(min(Ratio, 1 + ClipFactor), 1 - ClipFactor) .* Advantage;

% clipped surrogate loss
SurrogateLoss = -mean(min(Objective, ObjectiveClip));

% entropy loss: braching to improve performance
if EntropyLossWeight
    EntropyLoss = entropy(NetworkOutput, EntropyLossWeight);
else
    EntropyLoss = 0;
end

% total loss
Loss = SurrogateLoss + EntropyLoss;
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