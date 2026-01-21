%{
    This object represents a discretized Level Set Function
%}

classdef LSF_EvaluateGradientVersion
    properties
        
        % i,j -> (y, x)
        coordinateGrid {mustBeNumeric}
        
        % phi(i,j)
        values {mustBeNumeric}
        backwardGradients {mustBeNumeric}
        forwardGradients {mustBeNumeric}
        contour {mustBeNumeric}
        
        % dilation from [-1,1] x [-1, 1] to coordinate space
        % coordinateDilation = max(|x|) = max(|y|) 
        % helps define isomorphism coordinate space -> grid space
        coordinateDilation {mustBeNumeric} = 1
    end

    methods

        function obj = LSF_EvaluateGradientVersion(initialValues, dilation)
            
            n = size(initialValues,1);
            assert(n == size(initialValues,2));

            obj.values = initialValues;
            obj.coordinateDilation = dilation;
    
            [x,y] = meshgrid(1:n);
            x = (x - n/2) * 2/n * dilation;
            y = (n/2 - y) * 2/n * dilation;
            obj.coordinateGrid = cat(3,x,y);

        end

        function phi = interpolate(obj, coordinates)
            coordinateList = reshape(obj.coordinateGrid, [], 2);
            phi = griddata(coordinateList(1,:), coordinateList(2,:), obj.values(:), coordinates(1,:), coordinates(2,:));
        end

        % returns gridded points closest to 0 level set contour
        function obj = computeContour(obj)          
            obj.contour = rmoutliers( ...
                unique(contourc(obj.values, [0,0])', 'rows') ...
                , "percentiles", [20,80]);
            n = size(obj.values,1);
            obj.contour(:,1) = (obj.contour(:,1) - n/2) * 2/n;
            obj.contour(:,2) = (n/2 - obj.contour(:,2)) * 2/n;
        end
        
        % Creates a mesh for the PDEmodel with the provided function
        % Recomputes contour and interior points each time
        % Modifies the object, so capture the return
        function obj = updateMesh(obj, model, generateMeshFunction)
            n = size(obj.values,1);
            obj = obj.computeContour();
            boundedPoints = unique(cat(1,obj.contour, reshape((obj.values <= 0) .* obj.coordinateGrid, n^2,2)), "rows");           
            shp = alphaShape(boundedPoints);
            
            model.Geometry = [];
            geometryFromMesh(model, shp.Points', shp.alphaTriangulation()');
            generateMeshFunction(model);
        end

        % uses ENO
        % only run after updating LSF values, memory expensive
        function obj = computeGradients(obj)

            n = size(obj.values(),1);

            % square grid so arbitrary choice
            gridDelta = obj.coordinateGrid(1,2,1) - obj.coordinateGrid(1,1,1);
            

            forwardDiffx = [diff(obj.values,1,2)] / gridDelta;
            backwardDiffx = [zeros(n,1) forwardDiffx];
            forwardDiffx = [forwardDiffx zeros(n,1)];
            
            forwardDiffy = [-diff(obj.values,1)] / gridDelta;
            backwardDiffy = [forwardDiffy; zeros(1,n)];
            forwardDiffy = [zeros(1,n); forwardDiffy];
            
            obj.forwardGradients = zeros(n,n, 2);
            obj.backwardGradients = obj.forwardGradients;
            obj.forwardGradients = LSF_EvaluateGradientVersion.getWENOGradient(forwardDiffx, forwardDiffy);
            obj.backwardGradients = LSF_EvaluateGradientVersion.getWENOGradient(backwardDiffx, backwardDiffy);
        end

        function grads = evaluateGradientAtCoords(obj, coordinates, dir)
            
            n = size(obj.coordinateGrid(),1);
            
            % invert ij -> xy isomorphism
            gridPts = 2*obj.coordinateDilation*coordinates;
            gridPts = round([0,-1;1,0] * gridPts + n/2 + .5);
            grads = dir * obj.forwardGradients(gridPts(1,:), gridPts(2,:)) ...
                + (1-dir) * obj.backwardGradients(gridPts(1,:), gridPts(2,:));            
        end

        % The algorithms below are adapted from Distance Regularized Level Set Evolution and Its Application to Image Segmentation by Li et al.
        % Note gradients are WENO but divergence and laplace are inbuilts using central diff and finite diff respectively. 
        % These are numerically stable for DRLSE according to the cited paper.
        function obj = DRLSE(obj, eigenvalues, stationaryResults, gamma)
           
            obj = obj.reinitialize();
            obj = obj.computeContour();
            % obj = obj.NeumannBoundCond();
            obj = obj.computeGradients();
            
            [V, narrowband] = obj.getDeformationField(eigenvalues, stationaryResults, gamma);

            gradients = (V <= 0).* obj.forwardGradients + (V > 0) .* obj.backwardGradients;
            s = sqrt(gradients(:, :, 1).^2 + gradients(:, :, 2).^2);

            distRegTerm = distRegP2(obj, gradients, s);

            % add armijo wolfe linsearch for timestep

            % timestep = 1/max(abs(V),[],"all");
            timestep = LSF_EvaluateGradientVersion.armijoWolfeLinesearch();
            mu = 0;

            obj.values = obj.values + timestep * (mu*distRegTerm.*narrowband - V.*s);
        end

        function [V, narrowband] = getDeformationField(obj, eigenvalues, stationaryResults, gamma)
                    
            n = size(obj.values, 1);
            
            centralGrad = (obj.forwardGradients + obj.backwardGradients) ./ 2;
            centralGrad(:, 1, 1) = obj.forwardGradients(:, 1, 1);
            centralGrad(:, end, 1) = obj.backwardGradients(:, end, 1);
            centralGrad(1, :, 2) = obj.backwardGradients(1, :, 2);
            centralGrad(end, :, 2) = obj.forwardGradients(end, :, 2);

            magntidues = sqrt(centralGrad(:, :, 1).^2 + centralGrad(:, :, 2).^2);

            % avoid / by 0
            N = centralGrad ./ (magntidues + (magntidues == 0) * eps(0));

            % Get linear indices closest to boundary
            interior = obj.values <= 0; 
            cI = round(obj.contour * n/2);
            cI(:,1) = cI(:,1) + n/2;
            cI(:,2) =  n/2 - cI(:,2);

            % Select narrowband as every neighboring point of boundary
            cn = size(cI,1);
            col10 = [ones(cn,1), zeros(cn,1)];
            col01 = [zeros(cn,1),ones(cn,1)];
            cI = rmoutliers(unique(cat(1,cI,cI+1,cI-1,cI+col10, cI-col10,cI+col01,cI-col01),"rows"),"percentiles", [20,80]);
            cI = sub2ind([n n], cI(:,1), cI(:,2));
            
            narrowband = zeros(n);
            narrowband(cI) = 1;
            I = find(narrowband & interior);
            X = obj.coordinateGrid(:,:,1);
            Y = obj.coordinateGrid(:,:,2);

            % Evaluate eigenfunction gradients on interior narrowband
            [psix, psiy] = stationaryResults(1).evaluateGradient([X(I)';Y(I)']);
            psix(isnan(psix)) = 0;
            psiy(isnan(psiy)) = 0;
            fullpsix = zeros(n);
            fullpsix(I) = psix;
            fullpsiy = zeros(n);
            fullpsiy(I) = psiy;
            gradpsik = cat(3,fullpsix,fullpsiy);

            [psix, psiy] = stationaryResults(2).evaluateGradient([X(I)';Y(I)']);
            psix(isnan(psix)) = 0;
            psiy(isnan(psiy)) = 0;
            fullpsix = zeros(n);
            fullpsix(I) = psix;
            fullpsiy = zeros(n);
            fullpsiy(I) = psiy;
            gradpsik1 = cat(3,fullpsix,fullpsiy);
            
            % Boundary deformation speed given by
            % C(x) | dΩ = |Ω| ((1 − γ)|∂nψk|^2 + γ|∂nψk+1|^2) - ((1 − γ)λk + γλk+1)
            V = stationaryResults(1).Mesh.area * ((1-gamma) * dot(gradpsik, N, 3).^2 + gamma*dot(gradpsik1, N, 3).^2) ...
                - ((1-gamma)*eigenvalues(1) + gamma * eigenvalues(2));

            % Above calculation only valid for interior narrowband and
            % boundary
            % Eigenfunction gradients defined at same locations
            magpsik = gradpsik(:,:,1).^2 + gradpsik(:,:,2).^2;
            V = (magpsik~=0) .* V;
            
            % Compute closest "boundary" points, grid boundary may not
            % align with domain boundary so those also need extending
            [~, ind] = bwdist(magpsik~=0);           
            
            % Extend deformation field to exterior narrowband and grid
            % boundary
            V = V + ((narrowband & magpsik==0)  .* V(ind));
        end

        function distRegTerm = distRegP2(obj, gradients, s)

            % square grid so arbitrary choice
            gridDelta = obj.coordinateGrid(1,2,1) - obj.coordinateGrid(1,1,1);

            a=(s>=0) & (s<=1);
            b=(s>1);
            ps=a.*sin(2*pi*s)/(2*pi)+b.*(s-1);  
            dps=((ps~=0).*ps+(ps==0))./((s~=0).*s+(s==0));
            
            distRegTerm = LSF_EvaluateGradientVersion.div( ...
                dps.*gradients(:, :, 1) - gradients(:, :, 1), ...
                dps.*gradients(:, :, 2) - gradients(:, :, 2), ...
                gridDelta ...
                ) + 4*del2(obj.values, gridDelta);  
        end

        function obj = reinitialize(obj)
            n = size(obj.values,1);
            cI = rmoutliers(unique(round(obj.contour * n/2), "rows"), "percentile", [20,80]);
            cI(:,1) = cI(:,1) + n/2;
            cI(:,2) =  n/2 - cI(:,2);
            sign = ones(n) - 2 * (obj.values <= 0);
            obj.values = sign .* msfm2d(ones(size(obj.values)),cI',true,true);
        end


        function obj = NeumannBoundCond(obj)
            n = size(obj.values, 1);
            obj.values([1 n],[1 n]) = obj.values([3 n-2],[3 n-2]);  
            obj.values([1 n],2:end-1) = obj.values([3 n-2],2:end-1);          
            obj.values(2:end-1,[1 n]) = obj.values(2:end-1,[3 n-2]);  
        end
    end
    
    methods(Static)

        function gradient = getWENOGradient(diffx, diffy)
            
            % Optimal WENO weights for smooth regions
            w = [.1, .6, .3];
            n = size(diffx,1);

            % example for v1fx
            % dx1 dx2 dx3 dx4 dx5 ->  0 0 dx1 dx2 dx3 
            v1 = cat(3, [zeros(n,2) diffx(:, 1:n-2)], [diffy(3:n, :); zeros(2,n)]);
            v2 = cat(3, [zeros(n,1), diffx(:, 1:n-1)], [diffy(2:n, :); zeros(1,n)]);
            v3 = cat(3, diffx, diffy);
            v4 = cat(3, [diffx(:, 2:n) zeros(n,1)], [zeros(1, n); diffy(1:n-1, :)]);
            v5 = cat(3, [diffx(:, 3:n), zeros(n,2)], [zeros(2,n); diffy(1:n-2, :)]);

            % see page 34 of Level Set Methods Osher Fedkiw
        
            ENO1 = (v1/3 - 7/6 * v2 + 11/6 * v3);
            ENO2 = (-1/6 * v2 + 5/6 * v3 + v4/3);
            ENO3 = (v3/3 + 5/6 * v4 - v5/6);

            gradient = w(1) * ENO1 + w(2) * ENO2 + w(3) * ENO3;

            % Appropriate ENO for edge
            % gx
            gradient(:, 1, 1) = ENO3(:, 1, 1);
            gradient(:, end, 1) = ENO1(:, end, 1);
            % gy
            gradient(1, :, 2) = ENO1(1, :, 2);
            gradient(end, :, 2) = ENO3(end, :, 2);

        end

        function f = div(nx,ny, delta)
            [nxx,~]=gradient(nx, delta);  
            [~,nyy]=gradient(ny, delta);
            f=nxx+nyy;
        end

        % Performs a linesearch with Armijo-Wolfe 
        function alpha = armijoWolfeLinesearch()
            alpha =.1;
        end
    end
end