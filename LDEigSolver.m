function eigenResults = LDEigSolver(model, range, k)
    % dirichlet boundary conditions with u=0
    applyBoundaryCondition(model, "dirichlet", ...
        "Edge", 1:model.Geometry.NumEdges,"u", 0);

    % convert matlab internal equation −∇·(c∇u)+au=λdu to Δf = -λu
    specifyCoefficients(model,m=0,d=1,c=1,a=0,f=0);

    % solve, increasing eigenvalue range if kth val is not returned

    model.SolverOptions.MaxShift = realmax;

    eigenResults = solvepdeeig(model,range);
    %{
    if (k > 10^3)
        eigenvalues = sparse(eigenResults.Eigenvalues);
        eigenvectors = sparse(eigenResults.Eigenvectors);  
    else
    %}
        eigenvalues = eigenResults.Eigenvalues;
        eigenvectors = eigenResults.Eigenvectors; 
    %end
    
    while 1
        if (numel(eigenvalues) >= k)
            eigenResults = createPDEResults(model,eigenvectors,eigenvalues,"eigen");
            break;
        else
            range(1) = range(2);
            range(2) = range(2) + 20*(k);

            if (k > 500)
                [eigenvalues, eigenvectors] = parallelCalculate(model,range,eigenvalues,eigenvectors);
            else
                % directly calculate
                eigenResults = solvepdeeig(model,range);
                eigenvalues = cat(1,eigenvalues,eigenResults.Eigenvalues);
                eigenvectors = cat(2,eigenvectors,eigenResults.Eigenvectors);
            end

            % disp(append("Current size: ",num2str(numel(eigenvalues))));                            
        end
    end
end

function [eigenvalues, eigenvectors] = parallelCalculate(model,range,eigenvalues,eigenvectors)
    
    parallelTasks = int32(log(range(2)));

    rangeSpace = linspace(range(1),range(2),parallelTasks+1);
    collectedVals = cell(1,parallelTasks);
    collectedVecs = cell(1,parallelTasks);
    parfor i=1:parallelTasks
        blockRange = [rangeSpace(i),rangeSpace(i+1)];
        eigenResults = solvepdeeig(model,blockRange);
        collectedVals{i} = eigenResults.Eigenvalues;
        collectedVecs{i} = eigenResults.Eigenvectors;
        disp(append("Calculated Eigenvalues between ",num2str(blockRange)));
    end
    
    
    eigenvalues = cat(1,eigenvalues,collectedVals{:});
    eigenvectors = cat(2,eigenvectors, collectedVecs{:});
end