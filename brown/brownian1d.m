function [x,ev,dt] = brownian1d(in)
%BROWNIAN1D is a wrapper for using myrandomwalk for simulating a 1D
%brownian motion in phase space

in.save_evolution=true;
in.space_is_phase_space=true;

if isfield(in,'boundary_conditions')
    if length(in.boundary_conditions)>4 || ...
       sum(in.boundary_conditions(3:4)>0)
        in.boundary_conditions(3:4)=0;
        in.boundary_conditions(5:end)=[];
    end
end

if isfield(in,'viewbox_size')
    if length(in.viewbox_size)>4 
        in.viewbox_size(5:end)=[];
    end
end

mrwout=myrandomwalk(in);

x=mrwout.final_positions;

ev=mrwout.evolution;

dt=mrwout.in.time_step;

end

