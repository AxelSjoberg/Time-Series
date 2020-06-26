function model = removeUnsignificant(model, data, threshold, testWhite, plot)

%threshold = norminv(1-signlvl/2);
if nargin < 3
    threshold = 1;
end
if nargin < 4
    testWhite = 0;
end
if nargin < 5
    plot = 0;
end


res = filter(model.A,model.C,data);
models = {model};
FPE = model.Report.Fit.FPE;
white = montiTest(res);

if white == 0
    return
end

Alen = length(model.A);

run = 1;
if length(model.A) > 1
    Afree = model.Structure.A.Free;
end
if length(model.C) > 1
    Cfree = model.Structure.C.Free;
end

i = 1;
while run
    [pvec, pvec_sd] = getpvec(model);
    signvec = abs(pvec./pvec_sd);
    [minvalue, minindex] = min(signvec);
    
    if minvalue < threshold
        A = model.A;
        C = model.C;
        
        if minindex > Alen
            C(minindex - Alen + 2)     = 0;
            Cfree(minindex - Alen + 2) = 0;
        else
            A(minindex+1)     = 0;
            Afree(minindex+1) = 0;
        end    
        
        modelInit = idpoly(A, [], C);
        if length(model.A) > 1
            modelInit.Structure.A.Free = Afree;
        end
        if length(model.C) > 1
            modelInit.Structure.C.Free = Cfree;
        end
        model = pem(data, modelInit);
        
        i = i+1;
        res = filter(model.A,model.C,data);
        models{i} = model;
        FPE(i) = model.Report.Fit.FPE;
        white(i) = int8(montiTest(res));
    else
        run = 0;
    end
end

run = 1;

while run
    [minimum, i] = min(FPE);
    if white == 0
        models{i} = [];
    else
        run = 0;
    end
end

model = models{i};

if testWhite
    res = filter(model.A,model.C,data);
    whitenessTest(res);
end
if plot
    acfpacf(res);
end
end