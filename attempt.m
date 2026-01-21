function attempt(func)
try
    func();
catch E
    disp(E);
end