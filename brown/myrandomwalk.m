function out = myrandomwalk(in)
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

% my fields
dflt.save_evolution = false;
dflt.space_is_phase_space = false;

%% input handling and checks

if nargin == 0
    out = dflt;
    return;
end

% fill all missing fields from default
fname = fieldnames(dflt)
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
if isscalar(x)
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
    skip = inf; % do not plot in dimensions other that 2 or 3
else
    skip = in.plotframe_skips;
end

% custom settings short variables

sev = in.save_evolution;
sps = in.space_is_phase_space;

% custom settings consistency checks

if skip == Inf && nt == Inf
    error("The function thus set would cause an infinite loop.");
end

%% main

% setting space phase
if sps
    pickbound=@xdotbound;
    stepper=@phspacestep;
else
    pickbound=@xbound;
    stepper=@spacestep;
end

% defining variables only in nested function is deprecated in Matlab 2019a
Io=0;
Ie=0;

% allocate vector for store time evolution
if sev
    register=@writeev;
    if nt < Inf
        em=zeros(dim,nt,N); %evolution matrix
        % Puntatore ad aggiungi punti in posizione ennesima
    else
        em=zeros(dim,1,N);
        warning('WarningTAG:TagName', strcat ( 'The evolution matrix', ...
        ' is impossible to\n preallocate if the number of step is not', ...
        ' finite or predefined') );
        % Puntatore ad accoda matrice in coda
    end
else
    register=@nothing2;
end

% handle graphics
if skip < nt
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

jt = 0;

while jt < nt
    
    t = jt*dt;

    %make a step
    
    x=stepper(x);
    
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
       
    if mod(jt+1,skip+1) == 0
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
    jt = jt+1;
    % writes out on a vector
    em=register(em,x);
end

out.final_positions = x;
if sev
    out.evolution = em;
end
out.in = in;

if skip < nt
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
            x(Io|Ie,jdim) = -x(Io|Ie,jdim);
        end
    end

    function x=spacestep(x)
        xnew = feval(in.step_algorithm,B,t,x,dt);
        x = xnew + sqrt(dt)*rj(t,x);
    end

    function x=phspacestep(x)
        xpre = x(:,1:2:end);
        xnxt = feval(in.step_algorithm,B,t,xpre,dt);
        x(:,1:2:end) = xnxt + sqrt(dt)*rj(t,xpre);
        %Is a velocity
        x(:,2:2:end) = (xpre-x(:,1:2:end))./dt;
    end

    function em=writeev(em,x)
        em(:,jt,:)=x';
    end

    function nothing2(~,~)
        
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