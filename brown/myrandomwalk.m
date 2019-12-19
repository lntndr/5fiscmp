function out = myrandomwalk(in)
% MYRANDOMWALK is a custom version of randomwalk with added functionality
% as phase space managing and evolution saving.
% n-dimensional random walk
% in = randomwalk returns the default setup as a struct.
% out = randomwalk(in) returns the structure out with field 
% final_positions and field in containing the corresponding setup.

narginchk(0,1)

%% set defaults (2-dim case)

dflt.number_of_particles = 1e3;
dflt.viewbox_size = [-1 1 -1 1];
dflt.boundary_conditions = [0 0 0 0];
dflt.initial_positions = zeros(1e3,2);
dflt.drift_field = @(t,x)zeros(size(x));
dflt.random_jumps = @(t,x)randn(size(x));
dflt.time_step = 0.001;
dflt.time_span = inf;
dflt.step_algorithm = 'rk4';
dflt.plotframe_skips = 0;

% Custom fields
dflt.write_com_skips = 0;
dflt.space_is_phase_space = false;
dflt.stop_at_equilibrium = false;
dflt.sigma = 1; % Statistical equilibrium threshold

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

% fill short-named variables and perform some consistency check
box = in.viewbox_size;
dim = length(box)/2;
bc = zeros(1,2*dim); % default free boundary conditions
bc(1:length(in.boundary_conditions)) = in.boundary_conditions;

N = in.number_of_particles;
x = in.initial_positions;
if isscalar(x) && ~isa(x,'function_handle')
    x = x*ones(N,dim); 
end
if size(x,1) ~= N
    error(['size mismatch between number of particles ' ...
        'and initial positions'])
end
if size(x,2) ~= dim
    error('size mismatch between dimensions and initial positions')
end

B = in.drift_field;
rj = in.random_jumps;
dt = in.time_step;
nt = round(in.time_span/dt);
if dim ~= 2 && dim ~= 3
    skipf = inf; % do not plot in dimensions other that 2 or 3
else
    skipf = in.plotframe_skips;
end

% custom settings short variables

skipw = in.write_evolution_skips;
sps = in.space_is_phase_space;
sae = in.stop_at_equilibrium;
sgm = in.sigma;

% custom settings consistency checks

if skipf == Inf && nt == Inf && ~sae
    error("The function thus set would cause an infinite loop.");
end

%% main

% setting space phase
if sps
    pickbound=@xdotbound;
    stepper=@phspacestep;
    eqcheck=@phaseeq;
else
    pickbound=@xbound;
    stepper=@spacestep;
    eqcheck=@spaceeq;
end

if ~sae 
    eqcheck=@dontcheckeq;
end

% allocate vector for store center of mass and mean velocity evolution
if skipw < Inf
    register=@writeev;
    if nt < Inf
        em=zeros(floor(nt/(skipw+1)),dim); %evolution matrix
    else
        em=zeros(2,dim);
    end
else
    register=@nothing2;
end

% handle graphics
if skipf < nt
    shg
    clf
    set(gcf,'numbertitle','off','name','Random walk')
    if dim == 2
        h = plot(x(:,1),x(:,2),'.');
    else
        h = plot3(x(:,1),x(:,2),x(:,3),'.');
    end
    axis(box)
    axis square
    ttl = title(sprintf('t = %-8.1f',0));
    stop = uicontrol('style','toggle','string','stop',...
        'units','normalized','position',[.45 .01 .1 .05]);
end

% Service embedded variables
Io=0; % Logical vector for particles passing bo boundaries
Ie=0; % Logical vector for particles passing be boundaries
eq_period=10/dt;
ateq=false; %system is at equilibrium
jt = 0;

while jt < nt && ~ateq
    
    t = jt*dt;
    
    %boundary conditions
    for jdim = 1:dim
       jo = 2*jdim-1;
       je = 2*jdim;
       if mod(jdim,2)
           x=xbound(x);
       else
           x=pickbound(x);
       end
    end
    if mod(jt+1,skipf+1) == 0
        h.XData = x(:,1);
        h.YData = x(:,2);
        if dim == 3
            h.ZData = x(:,3);
        end
        ttl.String = sprintf('t = %-8.1f',jt*dt);
        drawnow
        if get(stop,'value')
            break
        end
    end
    if mod(jt+1,skipw+1) == 0
        em=register(em,x);
    end
    
    ateq=eqcheck(0);
    
    %make a step
    if ~ateq
        x=stepper(x);
    end
    
    jt = jt+1;
    
end


out.final_positions = x;

if skipw < Inf
    out.evolution = em;
end

out.stopped_at_eq = ateq;

out.in = in;

if skipf < nt
    set(stop,'string','close','value',0,'callback','close(gcf)')
end

% Nested functions
    
    function x=xbound(x)
        if bc(jo)
            Io = x(:,jdim) < box(jo);
            x(Io,jdim) = -x(Io,jdim) + 2*box(jo); 
        end
        if bc(je)
            Ie = x(:,jdim) > box(je);
            x(Ie,jdim) = -x(Ie,jdim) + 2*box(je);
        end
    end

    function x=xdotbound(x)
        if sum(Io|Ie)>0 % Some particle has bounced and velocities
            % change sign as consequence
            x(Io|Ie,jdim) = -x(Io|Ie,jdim);
        end
    end

    function x=spacestep(x)
        xnew = feval(in.step_algorithm,B,t,x,dt);
        x = xnew + sqrt(dt)*rj(t,x);
    end

    function x=phspacestep(x)
        %Evaluates analytical forces
        xb=feval(in.step_algorithm,B,t,x(:,2:2:end),dt);
        x(:,1:2:end)=x(:,1:2:end)+x(:,2:2:end)*dt;
        %Adds random component
        x(:,2:2:end)=xb + sqrt(dt)*rj(t,xb);
    end

    function em=writeev(em,x)
        em(ceil(jt+1/(skipw+1)),:)=mean(x,1);
    end

    function nothing2(~,~)
        
    end

    function eq=phaseeq(~)
        curr=ceil(jt+1/(skipw+1)); %current point
        if curr > eq_period %check if in the last eq_period
                            %speed has been under the variance of random
                            %steps for every single step and also in the
                            %mean value (a bit like autovelox + tutor)                          
            if abs(mean(em(curr-eq_period:curr,2:2:end))) < sgm && ...
                    abs(max(em(curr-eq_period:curr,2:2:end))) < sgm
                eq=true;
                % Cut allocated but unused space
                em(curr+1:end,:)=[];
            else
                eq=false;
            end
        else
            eq = false;
        end
    end

    function eq=spaceeq(~) %Not tested, basically a placeholder
        if abs(mean(diff(em(curr-eq_period:curr,:)),1)) < sgm
            eq=true;
            % Cut allocated but unused space
            em(curr+1:end,:)=[];
        else
            eq=false;
        end
    end

    function eq=dontcheckeq(~)
        eq=0;
    end

end

%% aux functions

function yp = euler(F,t,y,h) %#ok<DEFNU>
    yp = y + h*F(t,y);
end

function yp = rk4(F,t,y,h) %#ok<DEFNU>
    s1 = F(t,y);
    s2 = F(t+h/2,y+h/2*s1);
    s3 = F(t+h/2,y+h/2*s2);
    s4 = F(t+h,y+h*s3);
    yp = y + h*(s1 + 2*s2 + 2*s3 + s4)/6;
end