function [S,delta_v_vals] = make_v_peak_sphere(delta_v_max,N_v,v_0,sphere_density)
    if ~exist('sphere_density','var')
        sphere_density = 10;
    end
    if nargin < 1
        delta_v_max = 5 ;
    end
    
    if nargin < 2
        N_v = 2*delta_v_max + 1 ;
    end
    
    if nargin < 3
        v_0 = zeros(3,1) ;
    end
    
    delta_v_vals = linspace(0,delta_v_max,N_v) ;
    
    S = zeros(3,1)  ;
    for v = delta_v_vals(2:end)
        [x,y,z] = sphere(sphere_density*ceil(v)) ;
        x = v*x(:) + v_0(1) ;
        y = v*y(:) + v_0(2) ;
        z = v*z(:) + v_0(3) ;
        
        S = [S, [x,y,z]'] ;
    end
end