classdef SoftplusLayer < nnet.layer.Layer
    %   Softplus layer
    %
    % layer = rl.layer.SoftplusLayer(name) processes input data by
    % converting the values X according to the following formula:
    %
    % Y = log(1+exp(X))
    %
    % Examples:
    % layer = rl.layer.SoftplusLayer('Name','softplus')
    % layer = rl.layer.SoftplusLayer
    
    %   Copyright 2019 The MathWorks, Inc.
    
    methods
        function this = SoftplusLayer(varargin)
            
            parser = inputParser;
            addParameter(parser,'Name',"softplus",...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','Name'));
            addParameter(parser,'Description',getString(message('rl:general:SoftplusLayerDescription')),...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','Description'));
            parse(parser,varargin{:});
            r = parser.Results;
            
            this.Type = "SoftplusLayer";
            this.Name = r.Name;
            this.Description = r.Description;
        end
        function Z = predict(~,X)
            Z = rl.layer.SoftplusLayer.evaluate(X);
        end
    end
    
    methods (Static)
        function Z = evaluate(X)
            % Numerically stable implementation of softplus: Z = log(1 + exp(X))
            % When X < 0, Z = log(1 + exp(X)) 
            % When X >=0, Z = X + log(1 + exp(-X))
            
            Z = max(X,0) + log(1 + exp(-abs(X)));
            
            % the following implementation results in NaN gradient if exp(X) = Inf
            % Z = log(1 + exp(X));
            % Z(Z == Inf) = X(Z == Inf);
        end
    end
end
