function out = brownian1d(in)
%BROWNIAN1D is a wrapper for using myrandomwalk a simulator for a 1D
%brownian motion in phase space

dflt.number_of_particles = 5e4;
dflt.viewbox_size = [0 6 -5 5];
dflt.boundary_conditions = [1 0];
dflt.height = 10;
dflt.initial_positions=@(N,h)h*rand(N,1);
dflt.initial_velocities=@(N)zeros(N,1);
dflt.acceleration_field=@(t,x)-ones(size(x));
dflt.friction_coeff=1;
dflt.kT_over_mass=1;
dflt.random_kicks=@(t,x)randn(size(x));
dflt.time_step=1e-3;
dflt.time_span=100;
dflt.step_algorithm='rk4';
dflt.plotframe_skips = 24;
dflt.write_evolution_skips = 0;
dflt.stop_at_equilibrium = true;

% input handling and checks

if nargin == 0
    out = dflt;
    return;
end

% fill all missing fields from default
fname = fieldnames(dflt);
for jname = 1:length(fname)
    if ~isfield(in,fname{jname})
        in.(fname{jname}) = dflt.(fname{jname});
    end
end
%% input handling and checks

if nargin == 0
    out = dflt;
    return;
end

% fill all missing fields from default
fname = fieldnames(dflt);
for jname = 1:length(fname)
    if ~isfield(in,fname{jname})
        in.(fname{jname}) = dflt.(fname{jname});
    end
end

% Force myrandomwalk settings

in.save_evolution=true;
in.space_is_phase_space=true;

if isfield(in,'boundary_conditions')
    if length(in.boundary_conditions)>2
        in.boundary_conditions(3:end)=[];
    end
end

if isfield(in,'viewbox_size')
    if length(in.viewbox_size)>4 
        in.viewbox_size(5:end)=[];
    end
end

%short variables fill

N = in.number_of_particles;
gam = in.friction_coeff;
D = in.kT_over_mass/gam;
h = in.height;

% Fixed input setting
in.sigma=sqrt(in.time_step*2*in.kT_over_mass);
in.initial_positions=[feval(in.initial_positions,N,h), ...
    feval(in.initial_velocities,N)];
in.drift_field = @(t,x) -gam*x+in.acceleration_field(t,x);
in.random_jumps = @(t,x) sqrt(2*D)*gam*in.random_kicks(t,x);

out = myrandomwalk(in);

end

