classdef AbstractFigureVisualizerNew < handle
% ABSTRACTFIGUREVISUALIZER

% Copyright 2018 The MathWorks, Inc.

    properties (Access = protected)
        Figure
        Environment
        A
    end
    properties (Access = private)
        EnvUpdatedListener
        EnvDeletedListener
        FigDeletedListener
    end
    methods
        function this = AbstractFigureVisualizerNew(A)
            
            %this.Environment = env;
            this.A = A;
            
            %this.EnvDeletedListener = addlistener(...
            %    env,'ObjectBeingDestroyed',@(src,ed)delete(this));
            %this.EnvUpdatedListener = addlistener(...
            %    env,'EnvUpdated',@(src,ed)envUpdatedCB(this,src,ed));
            
            plot(this)
        end
        function delete(this)
            delete(this.EnvUpdatedListener);
            delete(this.EnvDeletedListener);
            delete(this.Figure);
        end
        function plot(this)
            if isempty(this.Figure) || ~isvalid(this.Figure)
                this.Figure = buildFigure(this);
            end
            updatePlot(this);
        end
        function bringToFront(this)
            figure(this.Figure);
        end
    end
    methods (Hidden)
        function f = qeGetFigure(this)
            f = this.Figure;
        end
    end
    methods (Access = protected,Abstract)
        f = buildFigure(this)
        updatePlot(this)
    end
    methods (Access = private)
        function envUpdatedCB(this,~,~)
            plot(this);
        end
    end
end