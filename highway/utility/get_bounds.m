function varargout = get_bounds(prob)
x = prob.x;
f = prob.f;
hX = prob.hX;
degree = deg(f);

if isfield(prob,'bound')
    bound = prob.bound;
else
    bound = 'both';
end

if isfield(prob,'gX')
    gX = prob.gX;
else
    gX = [];
end


prog = spotsosprog();
prog = prog.withIndeterminate(x);

[prog,lambda] = prog.newFree(1);

prog_ub = sosOnK(prog,lambda-f,x,hX,gX,degree);
prog_lb = sosOnK(prog,-lambda+f,x,hX,gX,degree);

options = spot_sdp_default_options();
options.verbose = 1;
options.domain_size = 1;
options.solveroptions = [];

if ~strcmp(bound,'lower')
    sol_ub = prog_ub.minimize(lambda, @spot_mosek, options);
    varargout{1} = double(sol_ub.eval(lambda));
    if ~strcmp(bound,'upper')
        sol_lb = prog_lb.minimize(-lambda, @spot_mosek, options);
        varargout{2} = double(sol_lb.eval(lambda));
    end
else
    sol_lb = prog_lb.minimize(-lambda, @spot_mosek, options);
    varargout{1} = double(sol_lb.eval(lambda));
end





end