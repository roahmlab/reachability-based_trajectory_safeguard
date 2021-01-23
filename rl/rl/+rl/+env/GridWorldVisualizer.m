classdef GridWorldVisualizer < handle
    %GridWorldVisualizer Creates a grid world visualizer
    %
    
    %   Copyright 2018-2019 The MathWorks, Inc.
    
    
    properties
        Obstacles %matrix of [r_1,c_1;r_2,c_2;...;r_n,c_n]
        CurrentState % [r,c]
        TerminalStates %matrix of [r_1,c_1;r_2,c_2;...;r_n,c_n]
        ShowTrace = false
        
    end
    
    properties
        Figure
        Ax
        CurrentStateMarker
        GWSize % [m,n]
        GridLines % line object
        TracedStates
    end
    
    methods
        function this = GridWorldVisualizer(GWSize,Obstacles)
            this.GWSize = GWSize;
            createView(this)
            setObstacles(this,Obstacles)
        end
        
        function createView(this)
            if isempty(this.Figure)
                this.Figure = figure('MenuBar','none','Toolbar','none');
                this.Ax = axes('Parent',this.Figure);%,'HandleVisibility','off');
                this.Ax.Toolbar.Visible = 'off';
                hold(this.Ax,'on') 
                this.Ax.Visible = 'off';
                axis(this.Ax,'equal')
            end
            drawGrid(this)
            this.CurrentStateMarker = patch(this.Ax,nan,nan,'r');
        end
        
        function drawGrid(this)
            m = this.GWSize(1);
            n = this.GWSize(2);
            delete(this.GridLines)
            xLineData = [];
            yLineData = [];
            y0 = 0.5;
            x0 = 0.5;
            for r = 0:m 
                yLineData = [yLineData;y0+r;y0+r;nan];
                xLineData = [xLineData;x0;n+0.5;nan];
            end
            for c = 0:n
                xLineData = [xLineData;x0+c;x0+c;nan];
                yLineData = [yLineData;y0;m+0.5;nan];
            end
            this.GridLines = plot(this.Ax,xLineData,-yLineData);
            
        end
        
        function setObstacles(this,Obstacles)
            delete(this.Obstacles)
            this.Obstacles = [];
            % Revisit and check dimensions
            for ct = 1:size(Obstacles,1)
                OS = eval(Obstacles(ct));
                r = OS(1);
                c = OS(2);
                NewPatch = localCreatePatch(this.Ax,r,c,'k');
                this.Obstacles = [this.Obstacles; NewPatch];
            end
        end
        
        function setTerminalStates(this,TerminalStates)
            delete(this.TerminalStates)
            this.TerminalStates = [];
            % Revisit and check dimensions
            for ct = 1:size(TerminalStates,1)
                TS = eval(TerminalStates(ct));
                r = TS(1);
                c = TS(2);
                NewPatch = localCreatePatch(this.Ax,r,c,'c');
                this.TerminalStates = [this.TerminalStates; NewPatch];
            end
        end
        
        function clearTrace(this)
            delete(this.TracedStates)
            this.TracedStates = [];
        end
        
        function setCurrentState(this,State)
            % REVISIT: check if valid state, including obstacles
            this.CurrentState = State;
            Data = 0.3*exp(1j*[0:.2:2*pi+.2]);
            XData = State(2) + real(Data);
            YData = -State(1) + imag(Data);
            set(this.CurrentStateMarker,'XData',XData,'YData',YData,'ZData',ones(size(YData)))
            if this.ShowTrace
                NewPatch = patch(this.Ax,nan,nan,'r');
                set(NewPatch,'XData',XData,'YData',YData,'ZData',0.5*ones(size(YData)))
                this.TracedStates = [this.TracedStates; NewPatch];
                for ct = 1:numel(this.TracedStates)
                    this.TracedStates(ct).FaceAlpha = 0.9*this.TracedStates(ct).FaceAlpha;
                end
            end
            
 
        end
        
        function animate(this,States)
            for ct = 1:size(States,1)
                setCurrentState(this,States(ct,:))
                drawnow
                pause(0.2)
            end
        end
 
        function delete(this)
            if ishandle(this.Figure)
                close(this.Figure)
            end
        end
    end
    
end

function NewPatch = localCreatePatch(Ax,r,c,Color)
    X = [c-0.5,c+0.5,c+0.5,c-0.5,c-0.5];
    Y = -[r-0.5,r-0.5,r+0.5,r+0.5,r-0.5];
    NewPatch = patch(Ax,X,Y,Color);
end

