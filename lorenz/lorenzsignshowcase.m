function lorenzsignshowcase
%LORENZSINGSHOWCASE documents the behaviour of lorenzsign by comparing it
%to lorenzsignp.p

spv = [99.65 100.5 160 350]; %Sparrow's periodic orbits rho values
grv = [20 24 25 28 40 80 95 97]; % Given non-periodic rho values

fprintf(['The Lorenz attractor is a chaotic system. As natural \n', ...
    '  consequence of this fact even a small difference in\n', ...
    '  evaluating the signature can produce different results\n', ...
    '  for a generic value of rho. However for some values the \n', ...
    '  trajectory is periodic. For Sparrow[1982] they are:\n\n']);

disp(spv);

fprintf(['We can observe that for this values we obtain the same\n', ...
    '  signature regardless of the method or the lenght. We are\n', ...
    '  also able to determine the period.\n\n']);

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

[pb,ppb]=lorenzsign(1e5,60*15);
fprintf('\n------------------   Signature   ------------------\n');
disp(cell2mat(ppb));
per=islorenzsignper(pb,60*15,40);
fprintf('\nThe signature shows a period of lenght = %d \n',per);
    
function askforuseraction
fprintf("Press any key to continue...\n");
pause;

function lorenztester(rho)

for k=2:13:15 %Ask for 2 and 15 lines of signature
    fprintf('********************\nIf N==%d\n********************\n',60*k);
    for c=1:size(rho,2)
        fprintf('\n          Rho==%d\n\n',rho(c));
        [sm,ssm,ssp,fd]=lorenzcompare(rho(c),60*k);
        fprintf('\n------------------   Signature   ------------------\n');
        disp(cell2mat(ssm));
        fprintf('\n---------------   Given signature   ---------------\n');
        disp(ssp);
        if isempty(fd)
            fprintf('\nAll the elements of the signatures are equal\n');
        else
            fprintf('\nSignatures start to differ from element #%d\n',fd);
        end
        per=islorenzsignper(sm,60*k,40);
        fprintf(' \nThe signature shows a period %d elements long\n',per);
    end
end

function [sm,ssm,ssp,fd]=lorenzcompare(rho,N)
% LORENZCOMPARE given a rho and a lenght, returns
% -the signature computed by lorenzsign.m
% -the easy-readable version of the signatures computed by lorenzsign and
%   lorenzsignp
% -If exists, the first element that differs in the signatures.
[sm,ssm]=lorenzsign(rho,N);     % lorenzsign.m signature
[sp,ssp]=lorenzsignp(rho,N);    % lorenzsignp.p signature
fd=find((sm~=sp),1); % first different element
