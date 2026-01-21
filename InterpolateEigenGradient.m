function gradpsi = InterpolateEigenGradient(model, eigenResults, kRange, queryPoints)
    %   Inputs:
    %   - model: PDE model
    %   - eigenResults: Eigenresults object from PDEtoolkit
    %   - k: eigenvalue number
    %   - queryPoints: Points where gradient should be interpolated (N x 2 matrix)
    %
    %   Output:
    %   - gx: grad x of L2 normalized eigenfunction at query points
    %   - gy: grad y of L2 normalized eigenfunction at query points
    
    % OLD CODE
    % xCoords = model.Mesh.Nodes(1, :);
    % yCoords = model.Mesh.Nodes(2, :);
    
    % xNodeCoords = xCoords(model.Mesh.Elements);
    % yNodeCoords = yCoords(model.Mesh.Elements);
    
    % centers = [mean(xNodeCoords, 1)' mean(yNodeCoords, 1)'];
    
    % intVals = eigenResults.interpolateSolution(centers',k);
    % [~,elemAreas] = area(model.Mesh);
    % L2Norm = sqrt(elemAreas*(intVals.^2));

    
    [~, M, ~] = assema(model,0,1,0);
    correctionFactor = 1.5;
    gradpsi = zeros(size(queryPoints,1),2,numel(kRange));

    for k=kRange
        eigenvector = eigenResults.Eigenvectors(:,k); 

        L2Norm = sqrt(eigenvector'*M*eigenvector)*correctionFactor;
      
        results = createPDEResults(model,eigenvector/L2Norm);
        
        [gradpsi(:,1,k), gradpsi(:,2,k)] = results.evaluateGradient(queryPoints');
        
    end
    
end