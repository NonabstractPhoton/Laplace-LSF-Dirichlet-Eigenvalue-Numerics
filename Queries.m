% Meant for use with objects from other scripts
% nxn LSF grid size
n = 128;

k = 5;
eigenvalues = vpa([18,36,56,78,78.6]);
pRange = 1:1000;
data = zeros(5,1000);
for p=pRange
    lowerSum = sum(eigenvalues(1:k).^p);
    data(:,p) = double((lowerSum^(1/p - 1)) * (eigenvalues(1:k).^(p-1)));
end

plot(pRange,data);
return;

% Creates an initial signed distance function LSF 
[X,Y] = meshgrid(1:n);
X = (X - n/2)*(2/n);
Y = (n/2 - Y) * (2/n);
A = min(cat(3,[sqrt((X-.2).^2 + (Y-.2).^2) - .32],[sqrt((X+.2).^2 + (Y+.2).^2) - .32]),[],3);
% A = X.^2+Y.^2 - sqrt(1/pi);

lsf = LSF(A, 1);
model = createpde();

meshGenerator = @(model) generateMesh(model,GeometricOrder="quadratic");
lsf = lsf.updateMesh(model, meshGenerator);

applyBoundaryCondition(model,"dirichlet", ...
                                         "Edge",1:model.Geometry.NumEdges, ...
                                         "u",0);
            specifyCoefficients(model,"m",0,...
                                      "d",0,...
                                      "c",1,...
                                      "a",0,...
                                      "f",-1);
            results = solvepde(model);
pdeplot(model,XYData=results.NodalSolution);
return;

lsf = lsf.DRLSE(0, 0, 1, 0);
return;

xNodes = model.Mesh.Nodes(1,:)';
yNodes = model.Mesh.Nodes(2,:)';
I = unique(boundary(xNodes, yNodes), "stable");
for i=1:length(I)
    [~,y]=ind2sub(size(model.Mesh.Elements), find(model.Mesh.Elements == I(i)));
    indices = model.Mesh.Elements(:,y);
    hold on;
    
    for j=1:size(indices,2)
        switch (j)
            case 1
                mkr = 'o';
            case 2
                mkr = '|';
            case 3
                mkr = '_';
            otherwise
                mkr = '*';
        end
        p = model.Mesh.Nodes(:, indices(:, j));
        scatter(p(1,:),p(2,:),2000,mkr);
    end
    nodes = model.Mesh.Nodes;
    thePoint = nodes(:, I(i));
    scatter(thePoint(1), thePoint(2), 3000,'magenta','x');
    if (i > 1)
        RNode = nodes(:,I(i+1));
        LNode = nodes(:, I(i-1));
        tangent =  LNode - RNode;
        normal = [-1*tangent(2); tangent(1)];
        quiver(thePoint(1), thePoint(2), tangent(1), tangent(2), 2, 'g');        
        quiver(thePoint(1), thePoint(2), normal(1), normal(2), .5, 'r');
        [gx, gy] = evaluateGradient(stationaryResult, thePoint(1), thePoint(2));
        quiver(thePoint(1), thePoint(2), gx, gy, .5,'b');
    end

    
    hold off;
    pause;

end