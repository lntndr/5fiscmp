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
dflt.step_algorithm = 'euler';
dflt.plotframe_skips = 0;

%% input handling and checks

if nargin == 0
    out = dflt;
    return;
end

% fill all missing fields from default
for fname = fieldnames(dflt)
    if ~isfield(in,fname)
        in.(fname) = dflt.(fname);
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

%% main

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
    xnew = feval(in.step_algorithm,B,t,x,dt);
    x = xnew + sqrt(dt)*rj(t,x);
    %x = xnew + sqrt(dt)*(xi(t,x) + xi(t+dt,xnew))/2;
    for jdim = 1:dim
        jo = 2*jdim-1;
        je = 2*jdim;
        if bc(jo)
            I = x(:,jdim) < box(jo);
            x(I,jdim) = -x(I,jdim) + 2*box(jo);
        end
        if bc(je)
            I = x(:,jdim) > box(je);
            x(I,jdim) = -x(I,jdim) + 2*box(je);
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
end

out.final_positions = x;
out.in = in;

if skip < nt
    set(stop,'string','close','value',0,'callback','close(gcf)')
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

