function out=brownian1dbench(friction,height)
lh=length(height);
lf=length(friction);

in.plotframe_skips=inf;
in.time_span=100;
in.stop_at_equilibrium=true;
in.timestep=1e3;
in.number_of_particles=1e4;
in.kT_over_mass=1;

results=cell(lh*lf,1);
tmat=zeros(length(friction)*length(height),3);

for k=1:lf
    % Set friction
    in.friction_coeff=friction(k);
    for j=1:lh
        index=lf*(k-1)+j;
        fprintf("Tentativo g=%d h=%d\n",friction(k),height(j));
        % Set height
        in.height=height(j);
        tic;results{index}=brownian1d(in);toc
        if results{index}.stopped_at_eq
            tmat(index,:)=[friction(k),height(j), ...
                length(results{index}.evolution)/in.timestep];
        else
            tmat(index,:)=[friction(k),height(j),NaN];
        end
    end
end

out.results=results;
out.bench_in=in;
out.friction=friction;
out.height=height;
out.eq_timemat=tmat;

