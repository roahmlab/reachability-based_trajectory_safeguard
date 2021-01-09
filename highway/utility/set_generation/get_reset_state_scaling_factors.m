function [z0_scale,z0_offset,Z0_reset_range] = get_reset_state_scaling_factors(f_fun,Z0_range,K_range,R,varargin)

%input 
%dynamics_fun: function handle (t,z,k) that returns the systems dynamics
%Z0_range: N_z x 2 vector with initial condition range [min,max]
%K_range: N_k x 2 vector with parameter range [min,max]
%R: function of (t,z and k) applied to the end time
%varargin:
N=1000; %number of random combos of (z,k) in Z0_range and K_range to generate
T=1; %Time to forward integrate
dt=0.1; %time discretization to check for scaling
n_scale = sqrt(2)/2*ones(size(Z0_range,1),1); %box to scale to (typically less than 1)
plot_on = false; %plotting on or off?
g = []; %no disturbance function
hZ0 = []; %no semialgebriac sets defining hZ0

%what this function does
%gives scaling factors for space such that|(z(t)+z_offset)/z_scale|<n_scale
%where z(t) = int_0^t dynamics_fun(t,z,k) dt, and t is in [0,T], and z(0)
%and k are sampled from Z0_range and K_range

%output (if not t_box_times given)
%z_offset: n_z x 1 double satisfying the above inequality
%z_scale: n_z x 1double satisfying the above inequality

%output (if t_box_times given)
%z_offset(:,1): n_z x 1 double satisfying the above inequality for all t
%z_scale(:,1): n_z x 1double satisfying the above inequality for all t
%z_offset(:,i+1): n_z x 1 double satisfying the above inequality for t =
%t_box_times(i)
%z_scale(:,i+1): n_z x 1double satisfying the above inequality for t =
%t_box_times(i)

for idx = 1:2:length(varargin)
    switch varargin{idx}
        case 'N'
            N = varargin{idx+1};
        case 'T'
            T = varargin{idx+1};
        case 'dt'
            dt = varargin{idx+1};
        case 'n_scale'
            n_scale = varargin{idx+1};
        case 'plotting'
            plot_on = varargin{idx+1};
        case 'g'
            g = varargin{idx+1};
        case 'hZ0'
            hZ0 = varargin{idx+1};
    end
end

N_z = size(Z0_range,1);
N_k = size(K_range,1);
N_param = N_k;

%if we have a disturbance function add it to the parameters and function
if ~isempty(g)
    K_range = [K_range;[-ones(N_z,1),ones(N_z,1)]];
    
    dynamics_fun = @(t,z,k) f_fun(t,z,k(1:N_k))+g(t,z,k(1:N_k)).*k((N_k+1):(N_k+N_z));
    
    N_k = N_k+N_z;
else
    dynamics_fun = f_fun;
end

%generate random combinations of z(0) and k in rectangluar sets
L = combinator(2,sum(N_k+N_z))';

combos_k = NaN(N_k,size(L,2));
for i = 1:N_k
    combos_k(i,:) = K_range(i,L(i,:));
end

combos_z0 = NaN(N_z,size(L,2));
for i = 1:N_z
    combos_z0(i,:) = Z0_range(i,L(N_k+i,:));
end

combos_k = [combos_k,randRange(K_range(:,1),K_range(:,2),[],[],1,N-size(L,2))];
combos_z0 =  [combos_z0,randRange(Z0_range(:,1),Z0_range(:,2),[],[],1,N-size(L,2))];

%if we have semi algebraic sets, check that all combos satisfy them, then
%add more until we reach at least N combos

if ~isempty(hZ0)
    
    num_combos = N-1;
    
    time_in_loop = tic;
    
    while num_combos<N
        
        combos_k = [combos_k,randRange(K_range(:,1),K_range(:,2),[],[],1,N)];
        combos_z0 =  [combos_z0,randRange(Z0_range(:,1),Z0_range(:,2),[],[],1,N)];
        
        
        for i = 1:length(hZ0)
            Lz0 = hZ0{1}(combos_z0,combos_k)>=0;
            combos_z0 = combos_z0(:,Lz0);
            combos_k = combos_k(:,Lz0);
        end
        num_combos = size(combos_z0,2);
        
        
        if toc(time_in_loop)>60
            warning('searching for points in semialg sets exceeded 1 min')
        end
        
    end
    
    %trim combo back down to N
    combos_k = combos_k(:,1:N);
    combos_z0 = combos_z0(:,1:N);
end

%forward integrate system and store results
Tvec = linspace(0,T,round(T/dt));

Z=NaN(N_z,size(combos_k,2));

parfor i = 1:N
    [~,temp] = ode45(@(t,z)dynamics_fun(t,z,combos_k(:,i)),Tvec,combos_z0(:,i));
    Z(:,i) = temp(end,:)';
end

if plot_on
    figure
    for i = 1:N_z
        subplot(N_z,1,i)
        hold on
        histogram(Z(i,:))
    end
    sgtitle('end states')
end


%apply reset map
Z_reset = NaN(size(Z));
combos_param = combos_k(1:N_param,:);

parfor i = 1:N
    Z_reset(:,i) = R(T,Z(:,i),combos_param(:,i));
end

if plot_on
    figure
    for i = 1:N_z
        subplot(N_z,1,i)
        hold on
        histogram(Z_reset(i,:))
    end
    sgtitle('end states reset')
end


z0_offset = NaN(N_z,1);
z0_scale = NaN(N_z,1);
max_z = NaN(N_z,1);
min_z = NaN(N_z,1);
    
%find max and mins and scaling factors

for i=1:N_z
    max_z(i,1) = max(Z_reset(i,:));
    min_z(i,1) = min(Z_reset(i,:));
    z0_offset(i,1) = -(max_z(i,1)+min_z(i,1))/2;
    z0_scale(i,1) = (max_z(i,1)-min_z(i,1))/(n_scale(i)*2);
end


if plot_on
    figure
    for i = 1:N_z
        subplot(N_z,1,i)
        hold on
        
        h = histogram((Z_reset(i,:)+z0_offset(i))/z0_scale(i));
        ylim = max(h.Values);
        plot(n_scale(i)*ones(size(2,1)),[0,ylim],'k--')
        hold on
        plot(-n_scale(i)*ones(size(2,1)),[0,ylim],'k--')
    end
    sgtitle('end states reset and scaled')
end


Z0_reset_range = [min_z,max_z];

end

