%{ 
    Performs Level Set Evolution to produce Cantidate minimizers for 
    the Dirichlet Eigenvalue Problem.
    Uses the FEM solution to Δu = -1 as the LSF for a given region.

    Inputs:
    - A: LSF representation of initial geometry in a unit meshgrid matrix 
    - k: Target eigenvalue index
    - displayMeshPer: Per how many iterations the Cantidate minimizer
      should be displayed
    - maxIterations: The maximum number of iterations before termination
%}
classdef LaplaceLSFEvolution < handle
    properties
        model;
        lsf;
        k {mustBeNumeric};

        targetRoot;
        
        computedMinimizer;

        currentJ {mustBeNumeric};
        historyJ {mustBeNumeric};

        currentResults;

        meshGenerator = @(model) LaplaceLSFEvolution.defaultGenMesh(model);

        autoCondense = false;

        pMode = false;
    end
    methods
        function obj = LaplaceLSFEvolution(k, targetRoot, p, useBidirectionalPNorm)
            obj.k = k;


            obj.model = createpde();

            obj.lsf = LaplaceLSF(obj.model,7,0);

            argCount = nargin;
            if (nargin == 4)
                obj.lsf.bidirectionalPNormWeights = useBidirectionalPNorm;
                argCount = 3;
            end

            switch (argCount)
                case 3
                    obj.lsf.setPEffective(p);
                    obj.pMode = true;
                    if (isempty(targetRoot))
                        obj.targetRoot = Cantidate.getPRoot(p);
                    else
                        obj.targetRoot = targetRoot;
                    end
                case 2
                    obj.targetRoot = targetRoot;
                case 1
                    obj.targetRoot = Cantidate.defaultRoot;
            end
            
        end

        function setGeometryViaLSFMatrix(obj,A)
            
            tempLsf = LSF(A, 1);
  
            % uses legacy code to generate initial geometry
            tempLsf.updateMesh(obj.model, obj.meshGenerator);

        end

        function setGeometryViaImage(obj,img)

            tempLsf = LSF(img, 1);

            % uses legacy code to generate initial geometry
            tempLsf.updateMeshStrict(obj.model, obj.meshGenerator);
        end

        function setGeometryViaPoints(obj,pts)
            % thresholds subject to change -> make dynamic
            shp = alphaShape(pts,"HoleThreshold",3*.0040, "RegionThreshold",3*.0010);
            obj.model.Geometry = [];
            geometryFromMesh(obj.model, shp.Points', shp.alphaTriangulation()');
            obj.meshGenerator(obj.model);
        end

        function start(obj,displayMeshPer,maxIterations, establishedMinimizer)      
            iterations = -1;
            
            while (iterations < maxIterations)
                %{
                if (mod(iterations, displayMeshPer) == 0)
                    obj.lsf.visualizationsLevel = 1;
                else
                    obj.lsf.visualizationsLevel = 0;
                end
                %}
            
                if (obj.k < 75)
                    obj.currentResults = LDEigSolver(obj.model, [0,20*obj.k], obj.k);
                else
                    obj.currentResults = LDEigSolver(obj.model, [0,20*sqrt(obj.k)], obj.k);
                end

                if (iterations >= 0)
                    obj.lsf.evolve(obj.currentResults, obj.k,obj.meshGenerator);
                end
                iterations = iterations + 1;

                if (obj.pMode)
                    obj.currentJ = obj.lsf.evaluatePNormObjectiveFunction(obj.currentResults.Eigenvalues,obj.currentResults.Mesh.area, obj.k); 
                else
                    obj.currentJ = obj.lsf.evaluateStandardObjectiveFunction(obj.currentResults.Eigenvalues,obj.currentResults.Mesh.area, obj.k);
                end

                if (~isinf(displayMeshPer))
                    disp("Eigenvalues:");
                    disp(obj.currentResults.Eigenvalues);
                    disp("Objective Function: "); 
                    disp(obj.currentJ);
                end

                if (obj.currentJ < establishedMinimizer.J)
                    obj.computedMinimizer = Cantidate(obj.model.Mesh.Nodes, obj.currentJ, obj.currentResults.Eigenvalues,obj.lsf.p);
                    establishedMinimizer = obj.computedMinimizer;
                end
            
                if (~isinf(displayMeshPer) && mod(iterations, displayMeshPer) == 0)
                    figure(1);
                    pdemesh(obj.model);
                    % figure(2);
                    % scatter(obj.lsf.contour(:,1), obj.lsf.contour(:,2))
                    title([num2str(iterations) ' Iterations']);
                    response = input('');
                    if (~(isempty(response)))
                        if (response == -1)
                            dbstop in LaplaceLSFEvolution.m at 143
                            disp('Debugger Triggered');
                        else
                            obj.lsf.visualizationsLevel = response;
                        end
                    end
                end


                checkConditionsAfter = 5;
                %{ 
                    Cancel iteration if, after a sample of given size,
                    a) J is monotonically increasing
                    b) J has a net increase and 
                        i) J is always more than three
                            standard deviations above the established minimizer's J
                        ii) The norm of the second difference across iterations falls below a given threshold
                %}
                obj.historyJ = [obj.historyJ obj.currentJ];
                diffJ = diff(obj.historyJ);
                if (mod(numel(obj.historyJ),checkConditionsAfter) == 0)
                    if (all(diffJ >= 0))
                        disp("Divering J detected");
                        if (isinf(displayMeshPer))
                            return;
                        end
                    elseif ((obj.historyJ(end) - obj.historyJ(1)) > 0)
                        if (all(obj.historyJ - 3*std(obj.historyJ) > establishedMinimizer.J))
                            disp("Insufficently small J detected");
                            if (isinf(displayMeshPer))
                                return;
                            end
                        elseif (norm(diff(diffJ)) < .5)
                            disp("Insufficient change in J detected");
                            if (isinf(displayMeshPer))
                                return;
                            end
                        end
                    end
                    obj.historyJ = [];
                end
            end
        end

        function delete(obj)
            if (~isempty(obj.computedMinimizer))
                
                holeThresh = 3*obj.model.Mesh.MaxElementSize^2;
                regionThresh = 3*obj.model.Mesh.MinElementSize^2;
                shp = alphaShape(obj.computedMinimizer.Points',"HoleThreshold",holeThresh, "RegionThreshold",regionThresh);

                if (shp.numRegions > 1)
                    obj.computedMinimizer.isDisjoint = true;
                end

                if (obj.pMode)
                    obj.computedMinimizer.p = obj.lsf.p;
                end

                obj.computedMinimizer.saveFor(obj.k,[],obj.targetRoot);
                if (obj.autoCondense)
                    Cantidate.condense(obj.k,obj.targetRoot);
                end
            end
        end
    end

    methods (Static)
        function defaultGenMesh(model)
            try
                generateMesh(model,GeometricOrder="quadratic",Hmin=0,Hmax=0.02,Hgrad=2);
            catch E
                disp(E)
                generateMesh(model,GeometricOrder="quadratic",Hgrad=2);
            end
        end
    end
end