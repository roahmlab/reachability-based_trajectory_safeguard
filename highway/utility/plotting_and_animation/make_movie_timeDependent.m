%% user parameters
filename = 'example_trial_timeDependent_2' ;

frame_rate = 30 ;

ti = 0 ;
tf = [] ;

display_time = true ;
save_video = false ;

%% automated from here
% get trajectory timing
t = A.time ;

if isempty(tf)
    tf = t(end) ;
end

t_total = tf - ti ;
Nframes = ceil(t_total*frame_rate) ;
tvec = linspace(ti,tf,Nframes) ;

% agent plotting setup
za = A.state ;
ta = A.time ;
ac = A.footprint_contour ;
head_x = 4*0.075*cos(linspace(0,2*pi,4)) ;
head_y = 4*0.05*sin(linspace(0,2*pi,4)) ;

% contour plotting setup
[X,Y] = makeContourAxes() ;
XY = makeContourAxes() ;
U = A.input ;
tU = A.input_time ;
if tU(end) < tf
    % if the input time doesn't last as long as the trajectory, append a
    % 'stopped' input for the last bit
    tU = [tU, tf] ;
    U = [U, [0;0]] ;
end
w = P1.FRS_sol.out.w ;
t = P1.FRS_sol.t ;
z = P1.FRS_sol.z ;
k = P1.FRS_sol.k ;
D = P1.FRS_sol.D ;
rcar = P1.FRS_sol.rcar ;
z0 = [P1.FRS_sol.x0 ; P1.FRS_sol.y0] ;

% world plotting setup
goal = W.goal ;
goal_x = W.goal_radius*cos(linspace(0,2*pi)) + goal(1) ;
goal_y = W.goal_radius*sin(linspace(0,2*pi)) + goal(2) ;
WB = boundsToContour(W.bounds) ;

% set up for video
fh = figure(312) ; hold on ; axis equal
if save_video
    vid = VideoWriter([filename,'.mp4'],'MPEG-4') ;
    vid.Quality = 100 ;
    vid.FrameRate = frame_rate ;
    
    open(vid) ;
end


for idx = 1:Nframes
    cla ;
    % get current time
    tidx = tvec(idx) ;
    if display_time
        text(W.bounds(1)+0.2,W.bounds(3)+0.2,['t = ',num2str(tidx-ti)])
    end
    
    % plot obstacles
    O = W.obstacles ;
    for oidx = 1:W.N_obstacles
        o = O(oidx) ;
        zo = o.state ;
        to = o.time ;
        
        % obstacle trajectory over time
        zto = matchTrajectories(ti+0.01:0.01:tidx,to,zo) ;
        plot(zto(1,:),zto(2,:),'--','Color',[1 0.7 0.3],'LineWidth',1.25)
        
        % obstacle at current location
        zoi = matchTrajectories(tidx,to,zo) ;
        oc = o.footprint_contour ;
        oc = oc + repmat(zoi,1,size(oc,2)) ;
        
        plot(oc(1,:),oc(2,:),'r','LineWidth',1.5) ;

    end
    
    % plot world
    plot(WB(1,:),WB(2,:),'r','LineWidth',1.5) ;
    plot(goal_x,goal_y,'k--','LineWidth',1) ;
    plot(goal(1),goal(2),'k*')
    
    % segway trajectory over time
    zidx = matchTrajectories(ti+0.01:0.01:tidx,ta,za) ;
    plot(zidx(1,:),zidx(2,:),'Color',[0.5 0.5 1],'LineWidth',1.25)
    
    % plot segway
    zidx = matchTrajectories(tidx,ta,za) ;
    hidx = zidx(3) ;
    plot(ac(1,:)+zidx(1),ac(2,:)+zidx(2),'b','LineWidth',1.5)
    
    % segway reach set slices
    uidx = matchTrajectories(tidx,tU,U,'previous') ;
    kidx = nan(2,1) ;
    kidx(1) = uidx(1)./A.wmax ;
    kidx(2) = (2/A.vmax).*(uidx(2) - (A.vmax/2)) ;
    
    wt = msubs(w,t,mod(tidx,P1.t_move)) ;
    wk = msubs(wt,k,kidx) ;
    Wc = reshape(full(msubs(wk,z,XY)),100,100) ;
    [wx,wy] = shiftScaleAndRotateContourAxes(X,Y,zidx,hidx,D,z0) ;
    contour(wx,wy,Wc,[1 1],'LineWidth',1.5,'Color',[0.2 0.8 0.4]) ;
    
    
    % plot heading
    [hxidx,hyidx] = rotxy(head_x,head_y,hidx) ;
    plot(hxidx+zidx(1),hyidx+zidx(2),'b','LineWidth',1.25)
    
    axis(W.bounds) ; 
    
    % save video or pause
    if save_video
        M(idx) = getframe(gca) ;
    else
        pause(0.001) ;
    end
end

if save_video
    writeVideo(vid,M) ;
    close(vid) ;
end
