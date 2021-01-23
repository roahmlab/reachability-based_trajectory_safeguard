classdef DoubleIntegratorVisualizer < rl.env.viz.AbstractFigureVisualizer
% DOUBLEINTEGRATORVISUALIZER

% Copyright 2018 The MathWorks, Inc.

    methods
        function this = DoubleIntegratorVisualizer(env)
            this = this@rl.env.viz.AbstractFigureVisualizer(env);
        end
    end
    methods (Access = protected)
        function f = buildFigure(this)
            
            f = figure(...
                'Toolbar','none',...
                'Visible','on',...
                'HandleVisibility','off', ...
                'NumberTitle','off',...
                'Name',getString(message('rl:env:VizNameDoubleIntegrator')),... 
                'MenuBar','none',...
                'CloseRequestFcn',@(~,~)delete(this));
            if ~strcmp(f.WindowStyle,'docked')
                f.Position(3:4) = [600 175];
            end
            ha = gca(f);
            
            ha.XLimMode = 'manual';
            ha.YLimMode = 'manual';
            ha.ZLimMode = 'manual';
            ha.DataAspectRatioMode = 'manual';
            ha.PlotBoxAspectRatioMode = 'manual';
            ha.YTick = [];
            
            env = this.Environment;
            
            md = env.MaxDistance;
            
            ha.XLim = [-md md];
            ha.YLim = [-0.1 0.5*md-0.1];
        end
        function updatePlot(this)

            env = this.Environment;
                
            f = this.Figure;
            ha = gca(f);
            
            md = env.MaxDistance;
            
            ha.XLim = [-md md];
            ha.YLim = [-0.1 0.5*md-0.1];
            
            x = env.State(1);
            w = 1;
            h = 1;
            
            pos = [x-w/2,0,w,h];
            
            rect = findobj(ha,'Tag','Mass');
            if isempty(rect)
                rectangle(ha,'Position',pos,'Tag','Mass','FaceColor','r');
            else
                rect.Position = pos;
            end
            drawnow();
        end
    end
end