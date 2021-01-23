classdef FcnLossLayer < nnet.layer.RegressionLayer
% Generic loss layer that takes in any loss functions (for layerGraph)

    properties
        % Forward function handle
        LossFcn
        
        % Whether the whole network has any RNN layer
        % Require for dlarray labelling because this layer is not aware of
        % RNN layers in the network
        IsNetworkStateful
    end
 
    methods
        function this = FcnLossLayer(varargin)           
            % Input parser
            parser = inputParser;
            addParameter(parser,'LossFcn',[],...
                @(val)validateattributes(val,{'function_handle'},{'scalar'},'','LossFcn'));
            addParameter(parser,'IsNetworkStateful',false);
            addParameter(parser,'Name',"genericLoss",...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','Name'));
            parse(parser,varargin{:});
            r = parser.Results;
            
            this.Type = "GenericLossLayer";
            this.Name = r.Name;
            this.LossFcn = r.LossFcn;
            this.IsNetworkStateful = r.IsNetworkStateful;
        end

        function Loss = forwardLoss(this, NetworkOutput, LossVariable)
            % Return the loss wrt the network output and the loss variables
            %
            % Inputs:
            %   NetworkOuput: Predictions made by network
            %   LossVariable: Contains any variable needed for loss
            %                 computation
            % Output:
            %   Loss - Loss computed from NetworkOuput and LossVariable
            
            % REVISIT: Might improve performance if we pre-infer the label,
            % only attach labels in this step.
            NetworkOutput = inferDlarrayLabel(this, NetworkOutput);
            
            Loss = feval(this.LossFcn, NetworkOutput, LossVariable);
        end
    end
    
    methods (Access = private)
        function Data = inferDlarrayLabel(this,Data)
            % Remove this when g2011546 is resolved
            switch ndims(Data)
                case 2
                    if this.IsNetworkStateful
                        Data = dlarray(Data,'CTB');
                    else
                        Data = dlarray(Data,'CB');
                    end
                case 3
                    if this.IsNetworkStateful
                        Data = dlarray(Data,'CBT');
                    else
                        % NOTE: stateless can still get here if there is a 
                        % single observation (fully connect with 2 outputs)
                        Data = dlarray(Data,'SSCB');
                    end
                case 4
                    Data = dlarray(Data,'SSCB');
                case 5
                    if this.IsNetworkStateful
                        % sequenceInputLayer with image input
                        Data = dlarray(Data,'SSCBT');
                    else
                        % image3dInputLayer
                        Data = dlarray(Data,'SSSCB');
                    end
                case 6
                    % sequenceInputLayer with image3d input
                    Data = dlarray(Data,'SSSCBT');
                otherwise
                    % TODO MSG
                    error('Generic loss layer does not support this output size.')
            end
        end
    end
end

