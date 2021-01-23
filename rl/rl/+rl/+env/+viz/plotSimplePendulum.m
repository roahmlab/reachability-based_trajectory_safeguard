function plotSimplePendulum(x,y)
% plot the simple pendulum given the x and y position of the mass

% Revised: 11-7-2018
% Copyright 2018 The MathWorks, Inc.

persistent f

L = sqrt(x^2 + y^2);
d = L/10;

if isempty(f) || ~isvalid(f)
    f = figure(...
        'Toolbar','none',...
        'NumberTitle','off',...
        'Name',getString(message('rl:env:VizNameSimplePendulum')),...
        'Visible','on',...
        'MenuBar','none');
    
    ha = gca(f);
    
    cla(ha);
    set(ha,'XLim',[-L,L].*1.5);
    set(ha,'YLim',[-L,L].*1.5);
    
    set(ha,'xticklabel',[]);
    set(ha,'yticklabel',[]);
    
    grid(ha,'on');
end

ha = gca(f);

rod  = findobj(ha,'Tag','rod');
mass = findobj(ha,'Tag','mass');

if isempty(rod)
    line(ha,'XData',[0,x],'YData',[0,y],...
        'Color','g','LineWidth',3,'tag','rod');
else
    set(rod,'XData',[0,x]);
    set(rod,'YData',[0,y]);
end
if isempty(mass)
    rectangle(ha,'Position',[x-d/2,y-d/2,d,d],'Curvature',[1,1],...
        'FaceColor','b','EdgeColor',[0,0,0],'LineWidth',3,'tag','mass');
else
    set(mass,'Position',[x-d/2,y-d/2,d,d]);
end
drawnow();

