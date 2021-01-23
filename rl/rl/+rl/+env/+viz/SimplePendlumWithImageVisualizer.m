classdef SimplePendlumWithImageVisualizer < rl.env.viz.AbstractFigureVisualizer
% SIMPLEPENDLUMWITHIMAGEVISUALIZER

% Revised: 11-8-2018
% Copyright 2018 The MathWorks, Inc.

    methods
        function this = SimplePendlumWithImageVisualizer(env)
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
                'Name',getString(message('rl:env:VizNameSimplePendulum')),...
                'MenuBar','none',...
                'CloseRequestFcn',@(~,~)delete(this));
            
            env = this.Environment;
            L = env.RodLength;
            
            h1 = subplot(1,3,1:2,'Parent',f);
            h2 = subplot(1,3,3  ,'Parent',f);
            h1.Tag = 'h1';
            h2.Tag = 'h2';
            
            cla(h1);
            daspect(h1,[1 1 1])
            set(h1,'XLim',[-L,L].*1.5);
            set(h1,'YLim',[-L,L].*1.5);
            
            set(h1,'xticklabel',[]);
            set(h1,'yticklabel',[]);
            
            grid(h1,'on');
        end
        function updatePlot(this)

            env = this.Environment;
            L = env.RodLength;
            d = L/10;
            theta = env.State(1);
            c = cos(theta);
            s = sin(theta);
            
            y =  L*c;
            x = -L*s;
                
            f = this.Figure;
            h1 = findobj(f,'Tag','h1');
            h2 = findobj(f,'Tag','h2');
            
            rod  = findobj(h1,'Tag','rod' );
            mass = findobj(h1,'Tag','mass');
            img  = findobj(f ,'Tag','img' );
            
            if isempty(rod)
                line(h1,'XData',[0,x],'YData',[0,y],...
                    'Color','g','LineWidth',3,'tag','rod');
            else
                set(rod,'XData',[0,x]);
                set(rod,'YData',[0,y]);
            end
            if isempty(mass)
                rectangle(h1,'Position',[x-d/2,y-d/2,d,d],'Curvature',[1,1],...
                    'FaceColor','b','EdgeColor',[0,0,0],'LineWidth',3,'tag','mass');
            else
                set(mass,'Position',[x-d/2,y-d/2,d,d]);
            end
            if isempty(img)
                img = imshow(generateImage(env),'Parent',h2);
                img.Tag = 'img';
            else
                img.CData = generateImage(env);
            end
            drawnow();
        end
    end
end