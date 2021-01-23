classdef Architecture < nnet.internal.cnn.analyzer.constraints.Constraint
    % Architecture  Cloned version of "Architecture" constraint for working
    %               with MIMO to bypass softmax/regression check
    % Borrowed from nnet.internal.cnn.analyzer.constraints.Constraint.Architecture()
    % TODO: reuse methods from nnet.internal.cnn.analyzer.constraints.Constraint.Architecture()
    
    %   Copyright 2019 The MathWorks, Inc.
    
    methods
        
        function testAtLeastOneInputLayer(test)
            % Test that the network has at least one input layer.
            %
            % At least one input layer is always a valid constraint for
            % both, DAG and series networks.
            
            isInput = [test.LayerAnalyzers.IsInputLayer];
            names = [test.LayerAnalyzers(isInput).Name]';
            
            if isempty(names)
                test.addIssue("E", "Network", [], ...
                    "Architecture:MissingInputLayer");
            end
        end
        
        function testAtLeastOneOutputLayer(test)
            % Test that the network has at least one output layer.
            %
            % At least one output layer is always a valid constraint for
            % both, DAG and series networks.
            
            isOutput = [test.LayerAnalyzers.IsOutputLayer];
            names = [test.LayerAnalyzers(isOutput).Name]';
            
            if isempty(names)
                test.addIssue("E", "Network", [], ...
                    "Architecture:MissingOutputLayer");
            end
        end
        
        function testClassificationMustBePrecededBySoftmax(test)
            % Test that all classification layers are preceded by a softmax
            % layer.
            %
            % The input to a classification layer has to be a probability
            % vector, which is what a softmax layer does.
            
            src = test.InternalConnections(:,1);
            dst = test.InternalConnections(:,3);
            
            for i=1:numel(test.LayerAnalyzers)
                if ~test.LayerAnalyzers(i).IsClassificationLayer
                    continue;
                end
                
                sources = src(dst == i);
                srcSoftmax = [test.LayerAnalyzers(sources).IsSoftmaxLayer];
                
                offending = sources(~srcSoftmax);
                offending = {test.LayerAnalyzers(offending).Name}';
                
                if ~isempty(offending)
                    test.addLayerError(i, ...
                        "Architecture:ClassificationMustBePrecededBySoftmax" );
                end
            end
        end
    end
end