%rlQValueRepresentation State-action value representation Q(o) or Q(o,a).

% Copyright 2019 The MathWorks, Inc.

classdef rlQValueRepresentation < rl.representation.rlAbstractQRepresentation
    
    properties (Dependent, SetAccess = private)
        ActionInfo
    end
    
    properties (Access = private)
        % delegator to handle different syntax and error checking of 
        % single-output or multi-output Q representation
        Delegate
    end
    
    methods
        function this = rlQValueRepresentation(Model, ObservationInfo, ActionInfo, Options)
            % Constructor
            
            this = this@rl.representation.rlAbstractQRepresentation(Model, ObservationInfo, ActionInfo, Options);
            
            % build delegator to handle single-output or multi-output
            ModelOutputSize = getSize(this.Model, 'output');
            NumOutput = prod(ModelOutputSize{1});
            if NumOutput == 1
                this.Delegate = rl.representation.qdelegate.QSingleOutputDelegate(this);
                
                % actions are the single output Q representation's inputs
                this.InputDataDimension = [this.InputDataDimension {this.ActionInfo.Dimension}];
                
                % single out Q rep not support multiple action specifications
                if numel(ActionInfo) > 1 && isa(ActionInfo,'rl.util.rlFiniteSetSpec')
                    error(message('rl:agent:errQRepNotSupportMultiDiscreteActionChannel'))
                end
            else
                this.Delegate = rl.representation.qdelegate.QMultiOutputDelegate(this);
            end
            
            % validate model input dimensions
            % - check if the model has the same number of input channel
            %   specified by obs specs (and action specs for single output
            %   Q)
            % - check if each output channel of the model has compatible
            %   size with obs specs dimension (and action specs for single
            %   output Q)
            validateModelInputDimension(this)
            
            % default to use mse loss
            this = setLoss(this,rl.util.getDefaultValueRepLoss());
        end
        
        function ActionInfo = get.ActionInfo(this)
            ActionInfo = this.ActionInfo_;
        end
    end
    
    methods (Hidden)
        function Type = getQType(this)
            if isa(this.Delegate, 'rl.representation.qdelegate.QSingleOutputDelegate')
                Type = 'singleOutput';
            else
                Type = 'multiOutput';
            end
        end
    end
    
    methods (Access = protected)
        function [QValue, State] = getValueImpl(this, varargin)
            [QValue, State] = getValue(this.Delegate, this, varargin{:});
        end
        
        function [MaxQ, MaxIndex, State] = getMaxQValueImpl(this, Observation)
            [MaxQ, MaxIndex, State] = getMaxQValue(this.Delegate, this, Observation);
        end
    end
end

