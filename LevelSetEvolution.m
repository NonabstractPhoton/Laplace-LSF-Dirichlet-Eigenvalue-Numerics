% nxn LSF grid size
n = 128;

% Creates an initial signed distance function LSF 
[X,Y] = meshgrid(1:n);
X = (X - n/2)*(2/n);
Y = (n/2 - Y) * (2/n);
A = min(cat(3,[sqrt((X-.2).^2 + (Y-.2).^2) - .32],[sqrt((X+.2).^2 + (Y+.2).^2) - .32]),[],3);
% A = X.^2+Y.^2 - sqrt(1/pi);

lsf = LSF(A, 1);
model = createpde();

meshGenerator = @(model) generateMesh(model,GeometricOrder="quadratic",Hmax=0.03);
lsf = lsf.updateMesh(model, meshGenerator);

% hyperparameters
% DRLSE relevant hyperparameters are specified within the method
k=1;
gamma = 1;
displayMeshPer = 5;
iterations = 0;
pastObjectiveFunctionValue = Inf;

% also change in internal cost function
p = 20;

while (1)
    eigenResults = solveIteration(model, k);
    lsf = lsf.DRLSE(eigenResults, model, k, gamma, meshGenerator, pastObjectiveFunctionValue);
    
    iterations = iterations + 1;

    objectiveFunctionValue = (1-gamma)*model.Mesh.area*eigenResults.Eigenvalues(k) ...
    + gamma*model.Mesh.area*(eigenResults.Eigenvalues(k)^p + eigenResults.Eigenvalues(k+1)^p)^(1/p);
    disp(objectiveFunctionValue);

    if (mod(iterations, displayMeshPer) == 0)
        figure(1);
        pdemesh(model);
        % figure(2);
        % scatter(lsf.contour(:,1), lsf.contour(:,2))
        title([num2str(iterations) ' Iterations']);
        pause;
    end
    
    if (abs(objectiveFunctionValue - pastObjectiveFunctionValue) < 10^-2)
        disp('small difference detected');
        % break;
    end

    pastObjectiveFunctionValue = objectiveFunctionValue;

end

% Assumes mesh already generated
function [eigenResults] = solveIteration(model, k)

    eigenResults = LDEigSolver(model, [0,20*(k+1)],k);

end