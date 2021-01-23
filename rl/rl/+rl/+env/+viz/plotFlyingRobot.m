function plotFlyingRobot(x,y,theta,tL,tR,Lx,Ly)
% plot the flying robot given the x, y, and theta

% PLOTFLYINGROBOT

% Revised: 11-7-2018
% Copyright 2018 The MathWorks, Inc.

persistent f

w = Lx/20;
h = Ly/30;

if isempty(f) || ~isvalid(f)
    f = figure(...
        'Toolbar','none',...
        'NumberTitle','off',...
        'Name',getString(message('rl:env:VizNameFlyingRobot')),...
        'Visible','on',...
        'MenuBar','none');
    
    ha = gca(f);
    localResetAxes(ha,Lx,Ly)
    
    grid(ha,'on');
    hold(ha,'on');
end

ha = gca(f);

c = cos(theta);
s = sin(theta);
R = [c -s;s c];
T = [R [x y]';zeros(1,3)];

V0 = [  -w -w  w  w w 0      0    ;
         h -h -h  h 0 h*1.5 -h*1.5;
        ones(1,7)   ];
V1 = T*V0;

vx = V1(1,1:4);
vy = V1(2,1:4);
ux = [x V1(1,5)];
uy = [y V1(2,5)];
wx = V1(1,6:7);
wy = V1(2,6:7);
tx = c*[tL tR];
ty = s*[tL tR];

body = findobj(ha,'Tag','body');
nose = findobj(ha,'Tag','nose');
quiv = findobj(ha,'Tag','quiv');

if isempty(quiv)
    quiver(ha,wx,wy,tx,ty,'Color','r','LineWidth',2,'Tag','quiv');
else
    quiv.XData = wx;
    quiv.YData = wy;
    quiv.UData = tx;
    quiv.VData = ty;
end
if isempty(body)
    patch(vx,vy,'y','Tag','body');
else
    body.XData = vx;
    body.YData = vy;
end
if isempty(nose)
    line(ux,uy,'Color','k','LineWidth',1,'Tag','nose');
else
    nose.XData = ux;
    nose.YData = uy;
end

drawnow();

function localResetAxes(ha,Lx,Ly)
cla(ha);
set(ha,'XLim',[-Lx,Lx]);
set(ha,'YLim',[-Ly,Ly]);

