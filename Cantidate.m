%{
    The Cantidate Object represents a computed minimizer
    The Cantidate Class provides a variety of functions to help
    find optimal minimizers
%}
classdef Cantidate

    properties (Constant)
        defaultRoot = "CantidateMinimizers/";
    end

    properties
        Points {mustBeNumeric}
        J {mustBeNumeric}
        Eigenvalues {mustBeNumeric}
        p {mustBeNumeric}   
        isDisjoint {mustBeNumericOrLogical} = false;
    end

    methods
        function obj = Cantidate(Points, J, Eigenvalues, p)
            obj.Points = Points;
            obj.J = J;
            obj.Eigenvalues = Eigenvalues;
            if (nargin == 4)
                obj.p = p;
            end
        end

        % Saves the Cantidate minimizer to the appropriate directory
        % Defaults fileName to the next index in the directory
        function saveFor(obj, k, fileName, root)

            if (nargin < 4)
                root = Cantidate.defaultRoot;
            end

            % Invalid or Analytically Known Values
            if (k <= 2 && strcmp(root,Cantidate.defaultRoot))
                disp("Invalid k");
                return;
            end

            name = Cantidate.getDir(k, root);
            directory = dir(name);
            if (nargin < 3 || isempty(fileName))
                name = append(name, num2str(1+sum(~[directory.isdir])),".mat");
            else
                name = append(name,fileName);
            end
            save(name, "obj");
        end
    end

    methods (Static)

        function name = getDir(k, root)
            if (nargin < 2)
                root = Cantidate.defaultRoot;
            end
            name = append(root, num2str(k),'/');
            if (~exist(name, 'dir'))
                mkdir(name);
            end
        end

        % View a plot with relevant information about a Cantidate 
        % minimizer for the kth eigenvalue. Defaults to the primary
        % Cantidate. Identifier is either a numeric for a p value or an
        % explicit file name
        function view(k,identifier,root)
            
            if (nargin < 3)
                root = [];
                if (nargin < 2)
                    identifier = [];
                end
            end

            c = Cantidate.load(k,identifier,root,nargin);

            if (isempty(c))
                error(append("Cantidate k=",num2str(k)," with identifier {",num2str(identifier),"} not found."));
            end

            shp = alphaShape(c.Points',"HoleThreshold",3*.0040,"RegionThreshold",3*.0010);
            figure();
            plot(shp);
            title(append(...
            "k = ", num2str(k), ...
            "    J = ",num2str(c.J),...
            "    Area = ",num2str(shp.area()), ...
            "    p = ",num2str(c.p)...
            ));               
        end

        function c = load(k,identifier,root,argCount)
            
            if (nargin < 4)
                argCount = nargin;
            end

            switch (argCount)
                case 1
                    c = Cantidate.loadInternal(k);
                case 2
                    if (isstring(identifier))
                        c = Cantidate.loadInternal(k,identifier);
                    else
                        c = Cantidate.loadInternal(k,'primary.mat',Cantidate.getPRoot(identifier));
                    end
                case 3
                    c = Cantidate.loadInternal(k,identifier,root);
                otherwise
                    error("Incorrect # of arguments");
            end   
        end

        % Loads the Cantidate object for a computed minimzer 
        % for the kth eigenvalue. Defaults to the primary Cantidate.
        function c = loadInternal(k, fileName, root)
            
            if (nargin < 3)
                root = Cantidate.defaultRoot;
            end

            c = [];

            name = Cantidate.getDir(k,root);
            if (nargin < 2 || isempty(fileName))
                name = append(name,'primary.mat');
            else
                name = append(name, fileName);
            end
            if (exist(name,"file"))
                c = struct2cell(load(name));
                c = c{1};
            end
        end

        % Condense potential computed minimzers to a primary, optimal one
        function condense(kRange,root)
        
            if (nargin < 2)
                root = Cantidate.defaultRoot;
            end

            for k=kRange    
                name = Cantidate.getDir(k,root);
                if (~exist(name,'dir'))
                    disp(append("No directory to condense for k=",num2str(k)));
                    continue;
                end
    
                directory = dir(name);
                files = directory(~[directory.isdir]);
                primary = Cantidate([], Inf,[]);
                primaryConnected = primary;
                for i=1:numel(files)
                    try
                        c = struct2cell(load(append(name,files(i).name)));
                    catch E
                        continue
                    end
                    c = c{1};
                    if (isa(c,"uint32"))
                        continue;
                    end
                    if (c.J <= primary.J)
                        primary = c;
                    end

                    if (~c.isDisjoint && c.J <= primaryConnected.J)
                        primaryConnected = c;
                    end
                end
                if (~isempty(primary.Points))
                    save(append(name,"primary.mat"), "primary");
                end
                if (~isempty(primaryConnected.Points) && primary.J ~= primaryConnected.J)
                    save(append(name,"primaryConnected.mat"), "primaryConnected");
                end
            end

        end

        % Deletes all but the optimal computed minimizer for the kth
        % eigenvalue
        function clean(kRange,root)
            if (nargin < 2)
                root = Cantidate.defaultRoot;
            end
            
            for k=kRange
                name = Cantidate.getDir(k,root);
                if (~exist(name,'dir'))
                    disp(append("No directory to clean for k=",num2str(k)));
                    continue;
                end
    
                directory = dir(name);
                files = directory(~[directory.isdir]);
                
                if (numel(files) == 0)
                    rmdir(name);
                    continue;
                end

                for i=1:numel(files)
                    if (~strcmp(files(i).name,"primary.mat") ...
                            && ~strcmp(files(i).name,"primaryConnected.mat"))
                        delete(append(name,files(i).name));
                    end
                end
            end
        end

        function cleanDirs(pRange,kRange)
            parfor i=numel(pRange)
                p = pRange(i);
                pRoot = Cantidate.getPRoot(p);
                if (dirsize(pRoot) == 0)
                    rmdir(pRoot,"s");
                else
                    Cantidate.clean(kRange,pRoot);
                end
            end
        end

        % Finds the optimal disjoint region for minimizing the kth eigenvalue
        function disjointOptimize(k)
            
            currentMinimizer = Cantidate([], Inf,[]);

            for i=1:int32(k/2)
                j = k-i;
                
                ci = Cantidate.load(i);
                cj = Cantidate.load(j);
            

                if (isempty(ci) || isempty(cj))
                    continue;
                end

                proposedJ = ci.J + cj.J;
                if (proposedJ <= currentMinimizer.J)
                    iPts = sqrt(ci.Eigenvalues(i)/proposedJ)*eye(2)*ci.Points;
                    jPts = sqrt(cj.Eigenvalues(j)/proposedJ)*eye(2)*cj.Points;
                    % Create Disjoint Geometry
                    jPts = jPts - [min(jPts(1,:)); 0] + [max(iPts(1,:)); 0] + [.1;0];
                    
                    unionPts = cat(2,iPts,jPts);
                    unionEigenvalues = sort(cat(1,...
                        double(proposedJ/ci.Eigenvalues(i) * ci.Eigenvalues),...
                        double(proposedJ/cj.Eigenvalues(j) * cj.Eigenvalues) ...
                        ));
                    currentMinimizer = Cantidate(unionPts, proposedJ, unionEigenvalues);
                end
                
            end

            if (~isempty(currentMinimizer.Points))
                currentMinimizer.isDisjoint = true;
                currentMinimizer.saveFor(k, 'disjointOptimized');
            end
        end

        function disjointOptimizePNorm(k,p)
            
            currentMinimizer = Cantidate([], Inf,[],p);   

            for i=1:int32(k/2)
                j = k-i;
                
                ci = Cantidate.load(i,p);
                cj = Cantidate.load(j,p);

                if (isempty(ci) || isempty(cj))
                    continue;
                end

                iNorm = sum(vpa(ci.Eigenvalues(1:i)).^p)^(1/p);
                jNorm = sum(vpa(cj.Eigenvalues(1:j)).^p)^(1/p);

                Vi = ci.J / iNorm;
                Vj = cj.J / jNorm;

                chi = Vi * jNorm / (Vj * iNorm);

                wi = double( 1/(1 + chi^(1/(p+1))) );
                wj = 1 - wi;
                
                proposedJ = double( (wi*Vi + wj*Vj)*(1/(wi^p)*iNorm^p + 1/(wj^p)*jNorm^p)^(1/p) );

                if (proposedJ <= currentMinimizer.J)
                    iPts = sqrt(wi)*eye(2)*ci.Points;
                    jPts = sqrt(wj)*eye(2)*cj.Points;
                    % Create Disjoint Geometry
                    jPts = jPts - [min(jPts(1,:)); 0] + [max(iPts(1,:)); 0] + [.1;0];
                    
                    unionPts = cat(2,iPts,jPts);
                    unionEigenvalues = sort(cat(1,1/wi * ci.Eigenvalues,1/wj * cj.Eigenvalues));
                    
                    currentMinimizer = Cantidate(unionPts, proposedJ, unionEigenvalues,p);
                    
                end
                
            end

            if (~isempty(currentMinimizer.Points))
                currentMinimizer.isDisjoint = true;
                currentMinimizer.saveFor(k, 'disjointOptimized',Cantidate.getPRoot(p));
            end
        end

        % Use to refine the primary Cantidate minimzer for the kth
        % eigenvalue 
        function refine(k,p,displayMeshPer,maxIterations,alphaMult,selectConnected,hmax,a,w,pEffective)
            
            
            if (nargin < 8)
                a = [];
                if (nargin < 7)
                    hmax = [];
                    if (nargin < 6)
                        selectConnected = false;
                        if (nargin <= 2)
                            displayMeshPer = 1;
                            maxIterations = Inf;
                            alphaMult = .5;
                        end
                    end
                end
            end

            if (nargin == 1 || isempty(p))
                if (selectConnected)
                    c = Cantidate.load(k,"primaryConnected.mat");
                else
                    c = Cantidate.load(k);
                end

                if (isempty(c))
                    error("No Cantidate found");
                end
    
                runtime = LaplaceLSFEvolution(k);

                if (~isempty(hmax))
                    runtime.meshGenerator = @(model) generateMesh(model,GeometricOrder="quadratic",Hmax=hmax,Hgrad=2);
                end

                if (isempty(a))
                    runtime.setGeometryViaPoints(c.Points');
                else
                    runtime.setGeometryViaPoints(Cantidate.perturbe(c.Points',a,w))
                end

                runtime.autoCondense = true;
                runtime.lsf.alphaMult = alphaMult;

                if (nargin >= 9)
                    runtime.lsf.setPEffective(pEffective);
                end
                runtime.start(displayMeshPer, maxIterations, c);
            else
                if (selectConnected)
                    c = Cantidate.load(k,"primaryConnected.mat",Cantidate.getPRoot(p));
                else
                    c = Cantidate.load(k,[],Cantidate.getPRoot(p));
                end

                if (isempty(c))
                    error("No Cantidate found");
                end
    
                runtime = LaplaceLSFEvolution(k,[],p);
                
                if (isempty(a))
                    runtime.setGeometryViaPoints(c.Points');
                else
                    runtime.setGeometryViaPoints(Cantidate.perturbe(c.Points',a,w))
                end

                runtime.autoCondense = true;
                if (~isempty(hmax))
                    runtime.meshGenerator = @(model) generateMesh(model,GeometricOrder="quadratic",Hmax=hmax);
                end
                runtime.start(displayMeshPer, maxIterations, c);
            end

        end

        % N x 2 
        function pts = perturbe(pts, a, w)
            
            pts = pts - mean(pts);
            theta = atan2(pts(:,2),pts(:,1));
            pts = pts + a*(pts ./ sqrt(pts(1,:).^2 + pts(2,:).^2) .* cos(w*theta));
        
        end

        function startFrom(k,cantidate,p,displayMeshPer, maxIterations)
            
            if (isempty(cantidate))
                error("Cantidate is empty");
            end
            
            if (nargin < 5)
                displayMeshPer = 1;
                maxIterations = Inf;
            end

            argCount = min(3,nargin);
            if (isempty(p))
                argCount = argCount - 1;
            end

            % This block is included in order to "trust" the previously 
            % computed eigenvalues stored in the the cantidate object
            
            checkPrior = numel(cantidate.Eigenvalues) >= k;
            switch (argCount)
                case 2
                    established = Cantidate.load(k);
                    runtime = LaplaceLSFEvolution(k);
                    runtime.setGeometryViaPoints(cantidate.Points');
                    if (checkPrior)
                        initialJ = runtime.lsf.evaluateStandardObjectiveFunction(cantidate.Eigenvalues,runtime.model.Mesh.area,k);
                    end
                case 3
                    established = Cantidate.load(k,p,[],2);
                    runtime = LaplaceLSFEvolution(k,[],p);
                    runtime.setGeometryViaPoints(cantidate.Points');
                    if (checkPrior)
                        initialJ = runtime.lsf.evaluatePNormObjectiveFunction(cantidate.Eigenvalues,runtime.model.Mesh.area,k);
                    end
            end

            if (isempty(established))
                established = Cantidate([],Inf,[]);
            end

            if (checkPrior && initialJ < established.J)
                cantidate.J = initialJ;
                runtime.computedMinimizer = cantidate;
                established = cantidate;
            end

            runtime.autoCondense = true;
            runtime.start(displayMeshPer,maxIterations,established);
            
        end

        % Creates initial geometry from a grayscale image and begins the
        % evolution process. The final two arguments are optional and
        % default to fully manual iteration.
        %
        % Image Requirements:
        % Black -> Region Interior
        % White -> Region Exterior
        % (Optional) Gray -> Region Boundary
        function startFromImg(k, imgPath, displayMeshPer, maxIterations)
            img = double(imread(imgPath));
            img = img(:, :, 1);
            
            dimDiff = size(img,1) - size(img,2);
            if (dimDiff > 0)
                img = cat(2,img,255+zeros(size(img,1),dimDiff));
            elseif (dimDiff < 0)
                img = cat(1,img,255+zeros(-dimDiff,size(img,2)));
            end
            
            runtime = LaplaceLSFEvolution(k);
            runtime.setGeometryViaImage(img);
            runtime.autoCondense = true;
            c = Cantidate.load(k);
            if (isempty(c))
                c = Cantidate([], Inf,[]);
            end
            if (nargin < 3)
                runtime.start(1,Inf,c);
            else
                runtime.start(displayMeshPer,maxIterations,c);
            end
        end
        
        % Tabulates the weights in the p-Norm based evolution process
        function coeffTable = tabulateCoeffecients(kRange, p, pEffective)
            if (nargin < 3)
                pEffective = p;
            end
            coeffTable = cell(numel(kRange),1);
            for i=1:numel(kRange)
                k = kRange(i);
                if (isempty(p))
                    c = Cantidate.load(k);
                else
                    c = Cantidate.load(k,p);
                end

                if (isempty(c))
                    continue;
                end

                eigenvalues = vpa(c.Eigenvalues);
                lowerSum = sum(eigenvalues(1:k).^pEffective);
                coeffTable{i} = (lowerSum^(1/pEffective - 1)) * (eigenvalues(1:k).^(pEffective-1));
            end
        end

        function plotCoeffecients(kRange, pVals, saveFigures)
            saveFigures = nargin > 2 && saveFigures;
            for k=kRange
                disp(append('Creating ',num2str(k)));
                fig = figure(k);
                data = zeros(k, numel(pVals));
                parfor i=1:numel(pVals)
                    p = pVals(i);
                    c = Cantidate.load(k,p);
                    if (isempty(c))
                        continue;
                    end
                    eigenvalues = vpa(c.Eigenvalues);
                    lowerSum = sum(eigenvalues(1:k).^p);
                    data(:,i) = double((lowerSum^(1/p - 1)) * (eigenvalues(1:k).^(p-1)));
                end
                plot(pVals,data);
                title(append("Coeffecient vs p for k=",num2str(k)));
                legend(string(num2cell(1:k)),'Location','northeastoutside');
                if (saveFigures)
                    saveas(fig,append("Figures/a_vs_p_k_",num2str(k)," pMax=",num2str(pVals(end)),".jpg"));
                end
            end
        end

        function slideshowView(kRange, pRange)
            for i=1:numel(kRange)
                k = kRange(i);
                if (nargin < 2)
                    Cantidate.view(k);
                    in = input("","s");
                    if (~isempty(in))
                        if (strcmp(in,"next"))
                            break;
                        end
                    end
                else
                    for p=pRange
                        Cantidate.view(k,p);
                        in = input("","s");
                        if (~isempty(in))
                            if (strcmp(in,"next"))
                                break;
                            end
                        end
                    end
                end
            end
        end

        function pRoot = getPRoot(p)
            pRoot = append(Cantidate.defaultRoot,"p_",num2str(p),"/");
        end

        
        function displayPercentDifferences(kRange)
        
            for k=[kRange]
                v = Cantidate.load(k).Eigenvalues;
                disp(append("(k-4):k, k=",num2str(k),": "));
                f=v((k-4):k);
                disp(f);
                disp(append("% Difference from k=:",num2str(k)));
                disp((f-f(end))/f(end) * 100);
            end
        end


        function localMins = ballLocalMinTest(kRange)
            
            tol = .6;
            maxIterations = 35;
            n = 256;
            
            % Set Initial Geometry Via Level Set Function
            [X,Y] = meshgrid(1:n);
            X = (X - n/2)*(2/n);
            Y = (n/2 - Y) * (2/n);
            A = X.^2+Y.^2 - 1/pi;
            
            localMins = [];

            runtime = LaplaceLSFEvolution(1);
            runtime.setGeometryViaLSFMatrix(A);
          
            ballEv = LDEigSolver(runtime.model,[0,20*kRange(end)],kRange(end)).Eigenvalues;
            
            parfor k=kRange
                c = Cantidate.load(k,'ball.mat');
                if (isempty(c))
                    c = Cantidate.load(k);
                    if (isempty(c))
                        c = Cantidate([], Inf, []);
                    end
                    runtime = LaplaceLSFEvolution(k);
                    runtime.setGeometryViaLSFMatrix(A);
                    try
                        runtime.start(Inf, maxIterations, c);
                    catch E
                        disp(E)
                        continue;
                    end

        
                    c = Cantidate(runtime.model.Mesh.Nodes,runtime.currentJ,runtime.currentResults.Eigenvalues);
                    c.saveFor(k,'ball.mat');
                end
                    
               
                diff = c.Eigenvalues(max(1,int32(k-log2(k))):k) - ballEv(max(1,int32(k-log2(k))):k);
                if (norm(diff) < tol + .5*10^(log10(c.Eigenvalues(k)) - 2))
                    localMins = [localMins, k];
                end
            end
        end

        function plotZeroSet(k,p,kRange,hmax,tol,saveFigures)
            if (isempty(p))
                runtime = LaplaceLSFEvolution(k);
                runtime.meshGenerator = @(model) generateMesh(model,GeometricOrder="quadratic",Hmin=0,Hmax=hmax,Hgrad=2);
                runtime.setGeometryViaPoints(Cantidate.load(k).Points');
            else
                runtime = LaplaceLSFEvolution(k,[],p);
                runtime.meshGenerator = @(model) generateMesh(model,GeometricOrder="quadratic",Hmin=0,Hmax=hmax,Hgrad=2);
                runtime.setGeometryViaPoints(Cantidate.load(k,p).Points');
            end

            
            results = LDEigSolver(runtime.model, [0,20*kRange(end)], kRange(end));

            for i=1:numel(kRange)
                currentK=kRange(i);
                figure(i);
                scattermat(results.Mesh.Nodes(:,abs(results.Eigenvectors(:,currentK)) < tol),'b.');
                axis equal;
                switch (currentK)
                    case 1
                        suffix = "st";
                    case 2
                        suffix = "nd";
                    case 3
                        suffix = "rd";
                    otherwise
                        suffix= "th";
                end
                title(append("Minimizer for k=",num2str(k),", ",num2str(currentK),suffix," Eigenfunction Zero Set"));
                if (nargin > 5 && saveFigures)
                    name = append("ZeroSetPlots/",num2str(k),"/");
                    if (~exist(name, 'dir'))
                        mkdir(name);
                    end
                    saveas(gcf,append(name,"k",num2str(k),"zs",num2str(currentK),".png"));
                end
            end
            
        end

        % returns difference of counting function and Weyl bounding term
        function testResult = polyaTest(k)
            c = Cantidate.load(k);
            % counting function
            Ndir = sum(c.Eigenvalues <= c.Eigenvalues(k));
            testResult = Ndir - (1/(4*pi) * c.Eigenvalues(k));
        end

        
    end
end