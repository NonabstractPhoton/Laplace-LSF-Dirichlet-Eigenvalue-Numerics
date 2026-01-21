% clutter deletion
close all;

% Grid Size for Initial Geometry
n = 512;

% Set Initial Geometry Via Level Set Function
[X,Y] = meshgrid(1:n);
X = (X - n/2)*(2/n);
Y = (n/2 - Y) * (2/n);
% A = min(cat(3,[sqrt((X-.3).^2 + (Y-.3).^2) - .32],[sqrt((X+.3).^2 + (Y+.3).^2) - .32]),[],3);
% A = X.^2+Y.^2 - 1/pi;

% hyperparameters 
displayMeshPer = Inf;
maxIterations = 150;
% parallelTasks = 10;
useInitMatrix = 1;
oscillationFrequencyRange = [.75,1.156,1.85];
pRange = [1];

% pRange = [1:60,80:20:1000,1050:50:3000,3100:100:10000];
% low error kRange
kRange = [44,47,49,50];
% 14 needs major improvement

% refineRange(kRange);


parfor i=1:(numel(kRange))
    k = kRange(i);
    for j = 1:numel(pRange)
        p = pRange(j);
        
        disp(append("On k=",num2str(k), ", p=",num2str(p)));
        
        
        %{
        attempt(@() Cantidate.startFrom(k,Cantidate.load(k),p,Inf,15));
    
        if ((p <= 40 || p == 60) && k~=17)
            attempt(@() Cantidate.refine(k,p,Inf,5,.3,false,.005));
        end
        continue;
        
        
        c =Cantidate.load(k);
        if (c.Eigenvalues(k) - c.Eigenvalues(k-1) > .001)
            attempt(@() Cantidate.refine(k,[],Inf,40,.5,.0125,.2,5));
            % attempt(@() Cantidate.startFrom(k,Cantidate.load(k-1),[],Inf,70))
            continue
        end
        %}
        k = kRange(i);
        Cantidate.refine(k,[],Inf,10,.5,0,.0025,[],[],2000);
        continue;
        % k = eigenvalueRange(randperm(numel(eigenvalueRange),1));
        w = oscillationFrequencyRange(randperm(numel(oscillationFrequencyRange),1));
        r0 = rand() * (sqrt(1/pi)-1/2) + 1/2;
        a = sqrt(2/pi - 2*r0^2);
    
        A = X.^2 + Y.^2 - (r0+a*sin(w*atan2(Y,X))).^2;
        
        c = Cantidate.load(k);
    
        if (isempty(c) || isa(c,"uint32"))
            c = Cantidate([],Inf,[]);
        end

        runtime = LaplaceLSFEvolution(k);
    
        %{
        if (~useInitMatrix)
            runtime.setGeometryViaPoints(c.Points');
        else
        %}
            runtime.setGeometryViaLSFMatrix(A);
        %end
    
    
        try
            runtime.start(displayMeshPer, maxIterations, c);
        catch E
            disp(E);
        end
    
        runtime.delete();
        
        %{
        try
            % Cantidate.startFrom(k,Cantidate.load(k),p,Inf,5);
            Cantidate.refine(k,p,Inf,10,.5,.03);
        catch E
            disp(E)
        end
         %}
    end
end
%}
for i=1:numel(kRange)
    k = kRange(i);
    Cantidate.condense(i);
    continue
    parfor j = 1:numel(pRange)
        p = pRange(j);
        % Cantidate.disjointOptimizePNorm(k,p);
        Cantidate.condense(k,Cantidate.getPRoot(p));
        Cantidate.clean(k,Cantidate.getPRoot(p));
        Cantidate.disjointOptimize(k);
    end
end


function refineRange(kRange)
    parfor j=1:numel(kRange)
    
        % meshSizes = [.01,.008,.0075,.006,.005,.006,.0075,.008,.009];
        % alphaMults = [.25,.2,.2,.15,.1,.1,.1,.1,.1];
        % iterCount = zeros(1,numel(meshSizes)) + 10;
        meshSizes = [.0018,.002,.0025,.006];
        alphaMults = [.5,.5,.4,.35];
        iterCount = [10,5,5,2];
        k = kRange(j);
        coriginal = Cantidate.load(k);
        for i=1:numel(meshSizes)
            disp(append("On k=",num2str(k),", hmax=",num2str(meshSizes(i))));
            c =Cantidate.load(k);
            iters = iterCount(i);
            attempt(@() Cantidate.refine(k,[],Inf,iters,alphaMults(i),false,meshSizes(i)+.002*(k ~= 9 && k~=13),[],[],2.5*10^3));
            c2 = Cantidate.load(k);
            disp("Eigenvalue delta");
            disp(c2.Eigenvalues(1:k) - c.Eigenvalues(1:k));
            disp("J delta");
            disp(c2.J - c.J);
        end
        disp("Net delta for ");
        disp(num2str(k));
        cend = Cantidate.load(k);
        disp("Eigenvalues delta:");
        disp(cend.Eigenvalues(1:k) - coriginal.Eigenvalues(1:k));
        disp("J delta");
        disp(cend.J - coriginal.J);
    end
end