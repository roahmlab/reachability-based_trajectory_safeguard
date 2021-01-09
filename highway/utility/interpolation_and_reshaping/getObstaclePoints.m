function P = getObstaclePoints(oa,ob,d,buf)
% Given a set of line segments defining obstacles, return points on those
% line segments spaced less than or equal to a given distance apart. If a
% buffer is required, get points on an oval around each line segment.
    if nargin < 4
        buf = 0 ;
    end
    
    [~,n] = size(oa) ;

    P = [] ;
    
    if buf == 0
        % If no buffer is required, get points evenly spaced on the lines
        % defined by oa and ob
        for idx = 1:n
            % for each obstacle, determine its length, then return points
            % spaced by distance d apart
            P = [P, getPointsOnLine(oa(:,idx),ob(:,idx),d)] ;
        end
    else
        % If a buffer is required, then, for each (oa,ob) pair, build an
        % oval around the line segment and get even spacing on that oval
        for idx = 1:n
            % extract the (x,y) coordinates from oa and ob
            oax = oa(1,idx) ;
            oay = oa(2,idx) ;
            obx = ob(1,idx) ;
            oby = ob(2,idx) ;
            
            % find the angle of the current line segment
            h = atan2(oby - oay, obx - oax) ;
            
            % create line segments parallel to the current line segment
            dx = buf*sin(h) ; dy = buf*cos(h) ;
            oax1 = oax + dx ; oax2 = oax - dx ;
            obx1 = obx + dx ; obx2 = obx - dx ;
            oay1 = oay - dy ; oay2 = oay + dy ;
            oby1 = oby - dy ; oby2 = oby + dy ;
            
            P1 = getPointsOnLine([oax1;oay1],[obx1;oby1],d) ;
            P2 = getPointsOnLine([oax2;oay2],[obx2;oby2],d) ;
            
            P = [P, P1, P2] ;
            
            % create 'caps' for the oval
            Pcirca = getPointsOnCircle(oax,oay,h+(pi/2),h+(3*pi/2),buf,d) ;
            Pcircb = getPointsOnCircle(obx,oby,h-(pi/2),h+(pi/2),buf,d) ;

            P = [P, Pcirca, Pcircb] ;
        end 
    end

    P = unique(P','rows','stable')' ;
end

function Pline = getPointsOnLine(pa,pb,dline)
% Given a line segment between points pa and pb, and a spacing distance
% dline, return points spaced evenly between pa and pb (inclusive) at a
% distance no greater than dline

    pax = pa(1) ;
    pay = pa(2) ;
    pbx = pb(1) ;
    pby = pb(2) ;

    l = sqrt((pax - pbx)^2 + (pay - pby)^2) ;

    nP = ceil(l/dline) ;

    svec = linspace(0,1,nP) ;

    Pline = [pax + svec.*(pbx - pax) ; pay + svec.*(pby - pay)] ;
end


function Pcirc = getPointsOnCircle(cx,cy,ha,hb,R,dcirc)
% Given a circular arc defined with a center (cx,cy), a radius R, and start
% and end angles ha and hb, return points spaced evenly between the
% endpoints of the arc at a distance no greater than dcirc.

    H = abs(hb - ha) ; % arc span in radians
    nP = ceil((R*H)/dcirc) ; % number of points required along arc
    dh = linspace(ha,hb,nP) ;
    
    Pcirc = [cx + R*cos(dh) ;
             cy + R*sin(dh) ] ;  
end