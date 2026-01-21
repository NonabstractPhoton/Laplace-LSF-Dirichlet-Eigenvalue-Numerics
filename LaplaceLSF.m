%{
    Uses the FEM solution to Δu = 1 as the LSF for a given region
%}
classdef LaplaceLSF < handle
    properties
        % USER PROVIDED
        % PDE model
        model;
        
        % INTERNAL
        % FEM solution to ∆u = 1
        LSFResults;
        extrapolateSolution;
        evaluateGradient;
        
        narrowbandPts {mustBeNumeric}
        narrowbandVals {mustBeNumeric}
        interiorPts {mustBeNumeric}
        interiorNarrowbandPts {mustBeNumeric}

        % Cost function p
        p {mustBeNumeric} = 1000;
        
        % Reduces performance but is necessary for accurate results with
        % large p
        useVariablePrecision {mustBeNumericOrLogical} = true;

        % Multiplier for the gradient descent dynamic coeffecient
        alphaMult {mustBeNumeric} = 1;
        
        % Additional bidirectional narrowband layers 
        nLayers {mustBeNumeric} = 1;
        layerPtCount {mustBeNumeric}

        % NOTE - Vestigial features
        visualizationsLevel {mustBeNumeric} = 0;
    
        useBidirectionalPNorm = false;
        useDirectV = true;

    end
    methods

        function obj = LaplaceLSF(model, nLayers,visLevel)
            obj.model = model;
            obj.nLayers = nLayers;
            obj.visualizationsLevel = visLevel;
        end

        % depends on accurate LSF
        % Let k = obj.layerPtCount, n = obj.nLayers 
        % The first n*k points are interior narrowband, the next k are the boundary, 
        % the final n*k points are the exterior narrowband
        function computeNarrowband(obj)
           
            boundary = abs(obj.LSFResults.NodalSolution) < 10^-3;
            boundaryPts = obj.LSFResults.Mesh.Nodes(:, find(boundary))';
            boundaryGrads = obj.evaluateGradient(boundaryPts);
            
            obj.layerPtCount = size(boundaryPts,1);

            obj.narrowbandPts = [];
            obj.narrowbandVals = [];
 
            % inner
            for i=obj.nLayers:-1:1
                extraPts = boundaryPts + boundaryGrads.*-i.*obj.model.Mesh.MinElementSize;
                obj.narrowbandPts = [obj.narrowbandPts; extraPts];
                obj.narrowbandVals = [obj.narrowbandVals; obj.LSFResults.interpolateSolution(extraPts')];
            end

            obj.narrowbandVals = [obj.narrowbandVals; obj.LSFResults.NodalSolution(boundary)];
            obj.narrowbandPts = [obj.narrowbandPts; boundaryPts];

            obj.interiorPts = setdiff(obj.LSFResults.Mesh.Nodes',obj.narrowbandPts,"rows");

            % outer
            for i=1:obj.nLayers
                extraPts = boundaryPts + boundaryGrads.*i.*obj.model.Mesh.MinElementSize;
                obj.narrowbandPts = [obj.narrowbandPts; extraPts];
                obj.narrowbandVals = [obj.narrowbandVals; obj.extrapolateSolution(extraPts)];
            end

            %{
            if (obj.visualizationsLevel == 2)
                figure();
                hold on;
                scatter(obj.narrowbandPts(:,1),obj.narrowbandPts(:,2));
                scatter(obj.interiorPts(:,1), obj.interiorPts(:,2));
                hold off;
                title('Narrowband (b) Interior (r)');
                pause;
            end
            %}
            
        end

        % depends on accurate mesh
        function computeLSF(obj)

            applyBoundaryCondition(obj.model,"dirichlet", ...
                                         "Edge",1:obj.model.Geometry.NumEdges, ...
                                         "u",0);
            specifyCoefficients(obj.model,"m",0,...
                                      "d",0,...
                                      "c",1,...
                                      "a",0,...
                                      "f",-1);
            obj.LSFResults = solvepde(obj.model);
            obj.extrapolateSolution = scatteredInterpolant(obj.LSFResults.Mesh.Nodes(1,:)', obj.LSFResults.Mesh.Nodes(2,:)',obj.LSFResults.NodalSolution, 'natural','linear');
            obj.evaluateGradient = scatteredInterpolant(obj.LSFResults.Mesh.Nodes(1,:)', obj.LSFResults.Mesh.Nodes(2,:)',[obj.LSFResults.XGradients, obj.LSFResults.YGradients], 'natural','boundary');
            %{
            if (obj.visualizationsLevel == 2)
                pdeplot(obj.model, XYData=obj.LSFResults.NodalSolution);
                title('LSF');
                figure();
                quiver(obj.model.Mesh.Nodes(1,:), obj.model.Mesh.Nodes(2,:),obj.LSFResults.XGradients', obj.LSFResults.YGradients');
                title('LSF Gradients');
                pause;
            end
            %}
        end
        
        % depends on accurate narrowband and LSF
        function updateMesh(obj, generateMeshFunction)             
            % hyperparameters
            chunkCount = 100;
            thicknessMin = (obj.model.Mesh.MinElementSize+obj.model.Mesh.MaxElementSize)/2;
            holeThresh = 3*obj.model.Mesh.MaxElementSize^2;
            regionThresh = 3*obj.model.Mesh.MinElementSize^2;

            pts = cat(1,obj.interiorNarrowbandPts,obj.interiorPts);
            trimmedPts = [];
            shp = alphaShape(pts,"HoleThreshold",holeThresh, "RegionThreshold",regionThresh);
            for i=1:shp.numRegions()
                regionPts = [];
                % Orient along principal axes
                [~,P] = shp.alphaTriangulation(i);
                [~,~,V] = svd(P,"econ");
                rotated = P*V;
                % Select longer axis
                mins = min(rotated);
                [axisLength,dim] = max(max(rotated) - mins);
                
                % [~,I] = sort(rotated(:,dim));
                % ordered = rotated(I,:);

                originLine = rotated(:,dim) - mins(dim);
                
                parfor j=0:(chunkCount-1)
                    chunk = (originLine >= axisLength*j/chunkCount) & (originLine <= axisLength*(j+1)/chunkCount);
                    chunkPts = rotated(chunk,:);
                    otherDim = (dim == 2)  + (dim == 1)*2;
                    if (max(chunkPts(:,otherDim)) - min(chunkPts(:,otherDim)) >= thicknessMin)
                        regionPts = [regionPts; chunkPts];
                    end
                end

                % rotate back
                trimmedPts = gather([trimmedPts; regionPts/V]);
            end
            
            %{
            if (obj.visualizationsLevel == 2)
                figure();
                hold on;
                scatter(pts(:,1), pts(:,2),'b');
                scatter(trimmedPts(:,1), trimmedPts(:,2),'g');
                hold off;
                title("Original (r) vs Trimmed (g)");
                pause;
            end
            %}

            % Rescale to region volume 1
            shp = alphaShape(trimmedPts,"HoleThreshold",holeThresh, "RegionThreshold",regionThresh);
            trimmedPts = ([1/sqrt(shp.area()), 0; 0, 1/sqrt(shp.area())] * trimmedPts')';
            
            shp = alphaShape(trimmedPts,"HoleThreshold",1/shp.area()*holeThresh, "RegionThreshold",1/shp.area()*regionThresh);
            obj.model.Geometry = [];
            geometryFromMesh(obj.model, shp.Points', shp.alphaTriangulation()');
            generateMeshFunction(obj.model);
        end

        function computeInteriorNarrowband(obj)
            
            intband = obj.narrowbandPts(obj.narrowbandVals < 0,:);

            %{
            if (obj.visualizationsLevel == 2)
                figure();
                hold on;
                scatter(intband(:,1), intband(:,2));
                title('Interior Narrowband Before Extension');
                pause;
            end
            %}
            
            pts = reshape(obj.narrowbandPts, obj.layerPtCount, 1+2*obj.nLayers, 2);
            vals = reshape(obj.narrowbandVals, obj.layerPtCount, 1+2*obj.nLayers);

            gaps = diff(sign(vals),1,2);
            gradVals = gradient(vals);
            gradPts = gradient(pts);

            maxDelta = obj.model.Mesh.MinElementSize;

            parfor i=1:size(vals,1)
                gPt = gradPts(i,:,:);
                gVl = gradVals(i,:);
                vl = vals(i,:);
                pt = pts(i,:,:);
                gap = find(gaps(i,:) ~= 0)';
                if (isempty(gap))
                    signgVl = sign(gVl);
                    
                    % Non constant gradient means first order analysis is
                    % not reliable if psi doesn't appear to cross 0
                    if (abs(sum(signgVl)) ~= numel(signgVl))
                        continue;
                    end
                    
                    s = sign(vl(1));
                    sG = signgVl(1);

                    % extrapolates left if heading left makes val -> 0
                    % extrapolates right if heading right makes val -> 0
                    gap = int8(s==sG) + int8(s~=sG)*size(vl,2);

                end
                % V + ∂V/∂x *t = 0
                % xBoundary = x + ∂x*t
                delta = gPt(1,gap,:).*(-vl(gap)./gVl(gap));
                norms = sqrt(delta(:,:,1).^2+delta(:,:,2).^2);
                point = pt(1,gap,:) + (norms > maxDelta) .* maxDelta .*delta./norms ...
                + (norms <= maxDelta).*delta;
                point(isnan(point)) = [];
                intband = [intband; reshape(point,[],2)];
            end

            %{
            if (obj.visualizationsLevel == 2)
                scatter(intband(:,1), intband(:,2),'r');
                title('Interior Narrowband After Extension');
                hold off;
                pause;
            end
            %}

            obj.interiorNarrowbandPts = intband;
            
        end

        function evolve(obj, eigenResults, k, meshGenerator)
            obj.computeLSF();
            obj.computeNarrowband();

            [V, s] = obj.getDeformationField(eigenResults, k);
            
            %{
            if (obj.visualizationsLevel > 0)
                figure();
                scatter(obj.narrowbandPts(:,1),obj.narrowbandPts(:,2),20,V);
                title('Narrowband and Deformation field');
                colorbar;
                pause;
            end
            %}

            alpha = obj.alphaMult*max(abs(obj.narrowbandVals),[],"all")/max(abs(V),[],"all");
            
            obj.narrowbandVals = obj.narrowbandVals - alpha.*s.*V;
            obj.computeInteriorNarrowband();
            obj.updateMesh(meshGenerator);

            %{
            if (obj.visualizationsLevel == 2)
                figure()
                hold on;
                scatter(obj.interiorPts(:,1), obj.interiorPts(:,2),'b');
                scatter(obj.interiorNarrowbandPts(:,1), obj.interiorNarrowbandPts(:,2),'r');
                title('Interior (b) and Interior Narrowband(r)');
                hold off;
                pause;
            end
            %}
            
        end
    
        function [V, s] = getDeformationField(obj, eigenResults, k)
          
            grads = obj.evaluateGradient(obj.narrowbandPts);
            
            % quiver(obj.narrowbandPts(:,1),obj.narrowbandPts(:,2),grads(:,1), grads(:,2));

            intbandRows = 1:obj.layerPtCount*(obj.nLayers+1);

            s = sqrt(sum(grads.^2,2));
            N = grads(intbandRows,:) ./s(intbandRows,:);
                     
            % Evaluate eigenfunction gradients on inner narrowband

            psi = InterpolateEigenGradient(obj.model,eigenResults,1:numel(eigenResults.Eigenvalues),obj.narrowbandPts(intbandRows,:));

            %{
            if (obj.visualizationsLevel > 0)
                figure();
                hold on;
                quiver(obj.narrowbandPts(intbandRows,1),obj.narrowbandPts(intbandRows,2),psik(:,1), psik(:,2),'b');
                quiver(obj.narrowbandPts(intbandRows,1),obj.narrowbandPts(intbandRows,2),psik1(:,1), psik1(:,2),'r');
                hold off;
                title('Eigenfunction Gradients (b=k, r=k+1)');
                pause;
            end
            %}

            % Deformation Field equation provided in image
           
            if (obj.useBidirectionalPNorm)
                V = obj.calculateBidirectionalPNormV(eigenResults,psi,N,k);
            elseif (obj.useDirectV)
                V = obj.calculateDirectV(eigenResults, psi, N, k);
            else
                V = obj.calculatePNormV(eigenResults,psi,N,k);
            end

            exterior = size(V,1):(size(V,1) + obj.layerPtCount * obj.nLayers);
            
            V = cat(1,V, zeros(obj.layerPtCount*obj.nLayers,1));
           
            % O(nlogn) V extension
            m = ceil(log2(obj.layerPtCount));
            m = m + (mod(m,2) == 1);
            distance = zeros(numel(exterior),2);
            distance(:,1) = Inf;
            tempDistance = distance;
            I = wrapIndex([exterior'+(m:-1:1), exterior'-(1:m)], obj.layerPtCount) + obj.layerPtCount*obj.nLayers;
            for i=1:size(I,2)
                tempDistance = [sum((obj.narrowbandPts(I(:,i),:) - obj.narrowbandPts(exterior,:)).^2,2), V(I(:,i))];
                improvements = (tempDistance(:, 1) < distance(:,1)) & (~isnan(tempDistance(:,2)));
                distance(improvements, :) = tempDistance(improvements, :);
            end

            V(exterior,:) = distance(:,2);

        end

        function V = calculateBidirectionalPNormV(obj,eigenResults, psi,N,k)
            if (obj.useVariablePrecision)
                eigenvalues = vpa(eigenResults.Eigenvalues,12);
                lowerSum = vpa(0);
                upperSum = vpa(0);
            else
                eigenvalues = eigenResults.Eigenvalues;
                lowerSum = 0;
                upperSum = 0;
            end

            p = obj.p;
            lambdaSum = (sum((eigenvalues(1:k).^p))^-1 + sum(eigenvalues((k+1):end).^(-p)));

            for i=1:k
                lowerSum = lowerSum + (eigenvalues(i)^(p-1))*(dot(psi(:,:,i),N,2).^2);
            end

            for i=((k+1):numel(eigenvalues))
                upperSum = upperSum + (eigenvalues(i)^(-p-1))*(dot(psi(:,:,i),N,2).^2);                      
            end
            
            V = double(-(lambdaSum^(-1/p)) + eigenResults.Mesh.area * lambdaSum^(-1/p - 1) * ...
                ((sum(eigenvalues(1:k).^p)^-2) * lowerSum + upperSum)); 
        end

        % Note variables with the same name between these two types of objective functions 
        % do not necessarily indicate the same quantity
        function V = calculatePNormV(obj,eigenResults, psi,N,k)
            if (obj.useVariablePrecision)
                eigenvalues = vpa(eigenResults.Eigenvalues,12);
                lowerSum = vpa(0);
            else
                eigenvalues = eigenResults.Eigenvalues;
                lowerSum = 0;
            end

            lambdaSum = sum(eigenvalues(1:k).^obj.p);
            
            for i=1:k
                lowerSum = lowerSum + eigenvalues(i)^(obj.p-1) * (dot(psi(:,:,i),N,2).^2);
            end

            V = double(-lambdaSum^(1/obj.p) + eigenResults.Mesh.area * lambdaSum^(1/obj.p-1) * lowerSum);
        end

        function V = calculateDirectV(obj, eigenResults, psi, N, k)
            V = -eigenResults.Mesh.area * (dot(psi(:,:,k),N,2).^2) + eigenResults.Eigenvalues(k);
        end

        % Note: Uses a modified objective function that helps account for
        % non-simple eigenvalues in a numeric context
        % This is the objective function used in domain evolution calculations
        function J = evaluateModifiedObjectiveFunction(obj, eigenvalues,volume, k)
            
            if (obj.useVariablePrecision)
                eigenvalues = vpa(eigenvalues,12);
            end
            J = double(volume * ((sum((eigenvalues(1:k).^obj.p))^-1 + sum(eigenvalues((k+1):end).^(-obj.p)))^(-1/obj.p)));
        end

        function J = evaluatePNormObjectiveFunction(obj,eigenvalues,volume,k)
            if (obj.useVariablePrecision)
                eigenvalues = vpa(eigenvalues,12);
            end
            J = double(volume * sum(eigenvalues(1:k).^obj.p)^(1/obj.p));
        end

        function J = evaluateStandardObjectiveFunction(~,eigenvalues,volume,k)
            J = volume*eigenvalues(k);
        end

        function setPEffective(obj,p)
            obj.p = p;
            if (p > 100)
                obj.useVariablePrecision = true;
            end
        end
    end
end
