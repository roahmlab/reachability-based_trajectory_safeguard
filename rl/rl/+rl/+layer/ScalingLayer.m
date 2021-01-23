classdef ScalingLayer < nnet.layer.Layer
    %   SCALING layer
    %
    % layer = rl.layer.ScalingLayer(name,scale,bias) processes input data by
    % converting the values X according to the following formula:
    %
    % Y = X.*scale + bias
    %
    % If bias is not provided, the default values is zero.
    %
    % Examples:
    % layer = rl.layer.ScalingLayer('Name','scalingWithBias','Scale',3,'Bias',-0.5)
    % layer = rl.layer.ScalingLayer('Name','scaling','Scale',3)
    
    %   Copyright 2018 The MathWorks, Inc.
        
    properties
        Scale
        Bias = 0
    end
    
    methods
        function this = ScalingLayer(varargin)
            
            parser = inputParser;
            addParameter(parser,'Scale',1);
            addParameter(parser,'Bias',0);
            addParameter(parser,'Name',"scaling",...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','Name'));
            addParameter(parser,'Description',getString(message('rl:general:ScalingLayerDescription')),...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','Description'));
            parse(parser,varargin{:});
            r = parser.Results;
            
            this.Scale = r.Scale;
            this.Bias = r.Bias;
            
            this.Type = "ScalingLayer";
            this.Name = r.Name;
            this.Description = r.Description;
        end
        function this = set.Scale(this,val)
            validateattributes(val,{'numeric'},{'finite'},'','Scale');
            this.Scale = val;
            checkParamConsistency(this);
        end
        function this = set.Bias(this,val)
            validateattributes(val,{'numeric'},{'finite'},'','Bias');
            this.Bias = val;
            checkParamConsistency(this);
        end
        function Z = predict(this,X)
            % Forward input data through the layer at prediction time and
            % output the result
            try
                Z = cast(X.*this.Scale + this.Bias,'like',X);
            catch ex
                if strcmp(ex.identifier,'MATLAB:dimagree')
                    error(message('rl:general:ScalingLayerPredictDimsDontAgree',...
                        mat2str(size(X)),mat2str(size(this.Scale)),mat2str(size(this.Bias))));
                else
                    rethrow(ex);
                end
            end
        end
        function dLdX = backward(this,~,~,dLdZ,~)
            % Backward propagate the derivative of the loss function through 
            % the layer
            try
                dLdX = cast(dLdZ.*this.Scale,'like',dLdZ);
            catch ex
                if strcmp(ex.identifier,'MATLAB:dimagree')
                    error(message('rl:general:ScalingLayerBackwardDimsDontAgree',...
                        mat2str(size(dLdZ)),mat2str(size(this.Scale)),mat2str(size(this.Bias))));
                else
                    rethrow(ex);
                end
            end
        end
    end
    methods (Access = private)
        function checkParamConsistency(this)
            m = this.Scale;
            b = this.Bias ;
            sz1 = size(m);
            sz2 = size(b);
            if ~isscalar(m) && ~isscalar(b) && ~isequal(sz1,sz2)
                error(message('rl:general:ScalingLayerInconsistentScaleBias',mat2str(sz1),mat2str(sz2)));
            end
        end
    end
end
