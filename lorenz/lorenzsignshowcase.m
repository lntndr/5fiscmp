function lorenzsignshowcase
spv = [99.65 100.5 160 350]; %Sparrow's periodic orbits rho values
grv = [20 24 25 28 40 80 95 97]; % Given non-periodic rho values

fprintf(['The Lorenz attractor is a chaotic system. As natural \n', ...
    '  consequence of this fact even a small difference in\n', ...
    '  evaluating the signature can produce different results\n', ...
    '  for a generic value of rho. However for some values the \n', ...
    '  trajectory is periodic. Some of this values are:\n\n']);

disp(spv);

fprintf(['We can observe that for this values we obtain the same\n', ...
    '  signature regardless of the method or the lenght.\n', ...
    'If ratio is equal to one, this means that the sequence of\n', ...
    '  symbols is identical, if it is equal to zero it means that\n', ...
    '  none of the symbols is equal.\n\n']);

askforuseraction;

lorenztester(spv);

fprintf(['Conversely for other values of rho using different methods\n' ...
    '  can produce different signature because of numerical noise or\n' ...
    '  roundoff errors. For example, by picking\n\n']);

disp(grv)

fprintf(['  we can observe differences in the signature for big\n', ...
    '  values of N.\n\n']);

askforuseraction;

lorenztester(grv);

fprintf(['If rho is much greater than the other parameters, e.g.\n', ...
    '  rho==1e5, the trajectory becomes pseudoperiodical, showing a\n', ...
    '  periodical signature even if the trajectory is not.\n\n']);

askforuseraction;

[~,ppb]=lorenzsigneskere(1e5,60*14);
fprintf('\n------------------   Signature   ------------------\n');
disp(cell2mat(ppb));
    
function askforuseraction
fprintf("Press any key to continue...\n");
pause;

function lorenztester(rho)
ratio=rho;
for k=2:3
    fprintf('********************\nIf N==%d\n********************\n',10^k);
    for c=1:size(rho,2)
        fprintf('\n          Rho==%d\n\n',rho(c));
        [ratio(c),ssm,ssp]=lorenzcompare(rho(c),10^k);
        fprintf('\n------------------   Signature   ------------------\n');
        disp(cell2mat(ssm));
        fprintf('\n---------------   Given signature   ---------------\n');
        disp(ssp);
        fprintf('\nRatio==%d\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n',ratio(c));
    end
end

function [r,ssm,ssp]=lorenzcompare(rho,N)
tic;[sm,ssm]=lorenzsign(rho,N);toc
tic;[sp,ssp]=lorenzsignp(rho,N);toc
s=sum(sm==sp);
r=s/size(sm,1);
