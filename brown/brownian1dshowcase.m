function brownian1dshowcase(b1din)
% BROWNIAN1DSHOWCASE gives a representation to the data made from
% browniand1dbench by fitting the final distributions, the evolution of
% the mean speed and of the center of mass and correlating the parameters
% gamma and height to the time needed to reach equilibrium.

%% basic check input
if ~isstruct(b1din) 
    error('Data must be a struct made by minibench.m');
else
    if ~iscell(b1din.results)
        error('Data must be a struct made by minibench.m');
    end
    %altre cose
end
%% Plot final distributions

% Allocate space for storing fit in.results
for j=1:length(b1din.results)
    figure('Name',strcat('g=', ...
        num2str(b1din.results{j}.in.friction_coeff),...
        ' h=',num2str(b1din.results{j}.in.height)))
    
    %sgtitle needs MATLAB>=2018b
    sgtitle(strcat('$\gamma=$', ...
        num2str(b1din.results{j}.in.friction_coeff),...
        ' $h=$',num2str(b1din.results{j}.in.height)),'Interpreter','latex')
    
    subplot(2,2,1)
        hold on
        title('$x$ final distribution','Interpreter','latex');
        xlabel('$x$','Interpreter','latex');
        ylabel('N of particles');
        nbin=ceil(2*length(b1din.results{j}.final_positions(:,1))^(1/3));
        histfit(b1din.results{j}.final_positions(:,1),nbin, ...
                'Exponential');
        hold off
        
    subplot(2,2,2)
        hold on
        title('$\dot{x}$ final distribution','Interpreter','latex');
        xlabel('$\dot{x}$','Interpreter','latex');
        ylabel('N of particles');
        histfit(b1din.results{j}.final_positions(:,2),nbin);
        hold off
        
    subplot(2,2,[3,4])
        hold on
        if b1din.results{j}.stopped_at_eq
            eqtoc=length(b1din.results{j}.evolution)*b1din.results{j}.in.time_step;
            title(sprintf('Evolution until eq. at t=%d',eqtoc));
        else
            neqtoc=b1din.results{j}.in.time_span;
            title(sprintf('Evolution until t=%d; Eq. not reached',neqtoc));
        end
        xlabel('Steps');
        plot(b1din.results{j}.evolution);
        yline(b1din.bench_in.kT_over_mass,':');
        legend({'$x$','$\dot{x}$','$\frac{K_BT}{m}$'},...
            'Interpreter','latex');
        hold off
end

%% Study correlation between t, gamma and height

tplot=b1din.eq_timemat;
nanindex=isnan(tplot(:,3)); 
tplot(nanindex,:)=[]; %removes data when eq is not reached
x=tplot(:,1)./tplot(:,2);tplot(:,1);
figure
hold on
title("Correlation beetween gamma/height and time to equilibrium");
scatter(x,tplot(:,3),'filled');
hold off

