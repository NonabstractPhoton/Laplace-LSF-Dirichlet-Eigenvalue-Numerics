% define omega0 as a disk with |omega0| = 1 
r = 1/sqrt(pi); 
center = [0; 0]; 
gd = [1; center; r]; 

% or a rectangle
% gd = [3; 4; -0.25; 0.25; 0.25; -0.25; -1; -1; 1; 1];

omega0 = decsg(gd);

% initialize
visualizations = true;
displayMeshPer = 30;
k = 2;
alpha = .5;
boundarydi = Inf;
[model, eigenResults, stationaryResult, boundaryIndices] = initializeProblem(omega0, k, visualizations);


% Iterate, Stopping Condition depends on norm of boundary translations vector
iterations = 0;
while (norm(boundarydi) > .001 || iterations == 100)

    iterations = iterations+1;
    visualizations = mod(iterations,displayMeshPer) == 0;

    boundarydi = varyBoundary(model, eigenResults, stationaryResult, ...
        boundaryIndices, k, alpha, visualizations);

    disp(norm(boundarydi));
    
    if (visualizations)
        pause;
        close all;
    end

    [eigenResults, stationaryResult, boundaryIndices] = solveIteration(model, k);
end

% kth eigenvalue
function [model, eigenResults, stationaryResult, I] = initializeProblem(initialGeometry, k, visualizations)
    
    model = createpde();
    
    geometryFromEdges(model, initialGeometry);
    
    % generate initial mesh
    generateMesh(model);
    
    if (visualizations)
        pdemesh(model);
        figure;
    end

    [eigenResults, stationaryResult, I] = solveIteration(model, k);
end


function [eigenResults, stationaryResult, I] = solveIteration(model, k)
    
    % dirichlet boundary conditions with u=0
    applyBoundaryCondition(model, "dirichlet", ...
        "Edge", [1:model.Geometry.NumEdges],"u", 0);
    
    % convert matlab internal equation −∇·(c∇u)+au=λdu to Δf = -λu
    specifyCoefficients(model,m=0,d=1,c=1,a=0,f=0);
    
    % solve, increasing eigenvalue range if kth val is not returned
    range = [0,100];
    while 1
        eigenResults = solvepdeeig(model, range);
        if (size(eigenResults.Eigenvectors,2) >= k)
            break;
        else
            range(2) = range(2) + 60;
        end
    end

    % calculate norm for later normalization
    l2norm = norm(eigenResults.Eigenvectors(:,k));
    
    % get boundary node indices
    xNodes = model.Mesh.Nodes(1,:)';
    yNodes = model.Mesh.Nodes(2,:)';
    I = unique(boundary(xNodes, yNodes), "stable");
    
    % create stationary solution object for normal gradient evaluation of
    % normalized solution
    tempModel = createpde();
    tempModel.Geometry = model.Geometry;
    tempModel.Mesh = model.Mesh;
    tempModel.BoundaryConditions = model.BoundaryConditions;
    
    % use Δf = -λu with computed lambda and "l2 normalized" interpolated u to rapidly converge
    % this is unfortunately necessary due to the limitations of the
    % EigenResults object

    function interpolant = interpolator(location, state, eigenResults, l2norm)
        interpolant(1, :) = interpolateSolution(eigenResults, location.x, location.y, k)/l2norm;
    end
    
    lambdau = @(location, state) eigenResults.Eigenvalues(k)*interpolator(location, state, eigenResults, l2norm); 
    
    % m * ∂^2u/∂t^2+d*∂u/∂t−∇·(c∇u)+au=f
    specifyCoefficients(tempModel,m=0,d=0,c=1,a=0,f=lambdau);

    setInitialConditions(tempModel, @(location) interpolator(location, 0, eigenResults, l2norm));
    
    stationaryResult = solvepde(tempModel);
    
    % optional extra visualizations
    %{
    pdeplot(model, "XYData", eigenResults.Eigenvectors(:,k));
    figure;
    pdeplot(tempModel, "XYData",stationaryResult.NodalSolution);
    figure;
    pdeplot(model,"XYData",stationaryResult.NodalSolution,"Contour","on","FlowData",[cx,cy]);
    figure;
    %}
end

function [boundarydi] = varyBoundary (model, eigenResults, stationaryResult, I, k, alpha, visualizations)
    
    % derivative of cost function
    % J(Ω) = |Ω|λk(Ω)
    % d(|Ω|)*λk(Ω) + d(λk(Ω))|Ω|
    % λk(Ω)*∫∂Ω V · n dσ + |Ω| * -∫∂Ω (∂uk/∂n)^2 * V · n dσ
    % Vi = vi * n, let vi be approximated by a hat function
    % λk(Ω)*2*(dist(Ni-1,Ni)+dist(Ni,Ni+1d)∫u(Ni)*t dt + |Ω| * -∫∂Ω (∂uk/∂n)^2 * V · n dσ
    % (dist(Ni-1,Ni)+dist(Ni,Ni+1)) * λk(Ω)*u(Ni) + 
    %   |Ω| * -(dist(Ni-1,N)*∫C1 (grad(u).n)^2 * (1-t) dt + dist(N,Ni+1)*∫C2 (grad(u).n)^2 * t))
    % use tangent vectors to parameterize path and get normals

    
    nodes = model.Mesh.Nodes;

    boundarydi = [];
    
    function [translation] = computeTranslation(LNode, thisNode, RNode)
        
        distTL = norm(LNode-thisNode);
        distRT = norm(RNode-thisNode);

        lu = eigenResults.Eigenvalues(k)*interpolateSolution(stationaryResult,thisNode(1), thisNode(2));
        
        tangentTL = LNode - thisNode;
        tangentRT = thisNode - RNode;

        normalTL = [-1*tangentTL(2); tangentTL(1)];
        normalTL = normalTL / norm(normalTL);

        normalRT = [-1*tangentRT(2); tangentRT(1)];
        normalRT = normalRT / norm(normalRT);

        pathTL = @(t) thisNode + tangentTL*t;
        pathRT = @(t) RNode + tangentRT*t;

        function val = integrandTL(t)
            point = pathTL(t);
            [gx, gy] = evaluateGradient(stationaryResult, point(1), point(2));
            val = (1-t) * dot([gx;gy], normalTL)^2;
            if (~isfinite(val))
                val = zeros(size(t));
            end
        end

        function val = integrandRT(t)
            point = pathRT(t);
            [gx, gy] = evaluateGradient(stationaryResult, point(1), point(2));
            val = (t) * dot([gx;gy], normalRT)^2;
            if (~isfinite(val))
                val = zeros(size(t));
            end
        end

        intTL = integral(@integrandTL,0,1);
        intRT = integral(@integrandRT,0,1);

        di = (distTL + distRT) * lu - eigenResults.Mesh.area * ...
            (distTL*intTL + distRT*intRT);

        boundarydi = cat(1, boundarydi, di);

        tangent = LNode - RNode;
        normal = [-tangent(2); tangent(1)];
        normal = normal / norm(normal);

        translation = -1*alpha*di*normal;
        
        % optional extra visualizations
        %{
        if (visualizations)
            scatter(thisNode(1), thisNode(2), 'r');
            scatter(LNode(1), LNode(2), 'b');
            scatter(RNode(1), RNode(2), 'b');
            quiver(thisNode(1), thisNode(2), normal(1), normal(2),.1,'red');
            quiver(thisNode(1), thisNode(2), gradient(1), gradient(2), 1,'blue');
        end
        %}
    end
    
    % Invariant to use as reference while calculating new boundary node values
    boundaryNodes = nodes(:,I);
    
    % First Boundary Point
    
    if (visualizations)
        hold on;
    end

    nodes(:, I(1)) = boundaryNodes(:, 1) ...
        + computeTranslation(boundaryNodes(:, end), boundaryNodes(:, 1), boundaryNodes(:, 2));
    
    % Each Other Boundary Point
    
    for i = 2:length(I)-1
        nodes(:, I(i)) = boundaryNodes(:, i) ...
            + computeTranslation(boundaryNodes(:, i-1), boundaryNodes(:, i), boundaryNodes(:, i+1));
    end
    
    % Last Boundary Point
    
    nodes(:, I(end)) = boundaryNodes(:, end) ...
        + computeTranslation( boundaryNodes(:, end-1), boundaryNodes(:, end), boundaryNodes(:, 1));
    
    if (visualizations)
        hold off;
        figure;
    end

    % Generate new geometry and mesh
    model.Geometry = [];
    TR = delaunayTriangulation(nodes');
    geometryFromMesh(model, TR.Points', TR.ConnectivityList');
    newBoundary = unique(freeBoundary(TR));
    generateMesh(model,GeometricOrder="quadratic");
    if (visualizations)
        pdemesh(model);
    end
end


