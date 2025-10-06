# run files for test pilots
using Distributed;
#addprocs(3);
@everywhere include("./src/loadMod.jl");
@everywhere include("./src/decomp_func.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

@everywhere caseList = [13]; 
@everywhere ci = 1;
@everywhere T = 24;

@everywhere fData, bData, dData, sigmaData, gData = readInData(ci, caseList, T);

# generate the random demand portfolio
J = 1:5;
df_Q = CSV.read("./data/$(caseList[ci])/case$(caseList[ci])_Q.csv", DataFrame);
df_r = CSV.read("./data/$(caseList[ci])/case$(caseList[ci])_rdis.csv", DataFrame);

# set up the charger data
NCP = [2,3,5,9,12,13];
c = Dict();
f = Dict();
xbar = Dict();
x0 = Dict();
u0 = Dict();
for i in NCP
    xbar[i] = 15.0;
    x0[i] = 0.0;
    u0[i] = 0;
end
f[2] = 10.0;
f[9] = 12.0;
f[3] = 25.0;
f[5] = 23.0;
f[12] = 15.0;
f[13] = 26.0;
c[2] = 1.0;
c[9] = 1.0;
c[3] = 1.5;
c[5] = 1.5;
c[12] = 1.2;
c[13] = 1.8;

current_charger = [2,9,12];
for i in current_charger
    x0[i] = 1.0;
    u0[i] = 1;
end
ρ = 30;

# calculate r_ij and set Delta
r = Dict();
for j in J
    for i in NCP
        r[i,j] = 100/df_r[i,"x$j"];
    end
end
Delta = 0.1;

# collect the Q data from df_Q
Q = Dict();
for j in J
    for t in 1:T
        Q[j,t] = df_Q[t,"x$j"];
    end
end

# obtain demand parameters
D = Dict();
for i in fData.IDList
    if i in keys(dData.pd)
        for t in 1:T
            D[i,t] = dData.pd[i][t * 4 - 3];
        end
    else
        for t in 1:T
            D[i,t] = 0.0;
        end
    end
end

# obtain generation upper and lower bounds
Pub = Dict();
Plb = Dict();
g = Dict();
for i in fData.IDList
    if i in keys(fData.LocRev)
        Pub[i] = sum(fData.Pmax[g] for g in fData.LocRev[i]);
        Plb[i] = sum(fData.Pmin[g] for g in fData.LocRev[i]);
        g[i] = maximum([fData.cp[g].params[length(fData.cp[g].params) - 1] for g in fData.LocRev[i]]);
    else
        Pub[i] = 0.0;
        Plb[i] = 0.0;
        g[i] = 0.0;
    end
end

# collect all branches
brList = [];
for item in fData.brList
    if !((item[1], item[2]) in brList || (item[2], item[1]) in brList)
        push!(brList, (item[1], item[2]));
    end
end
Wbar = Dict();
θdiff = Dict();
b = Dict();
for item in brList
    Wbar[item] = fData.rateA[(item[1], item[2], 1)];
    θdiff[item] = pi/6;
    b[item] = -fData.b[(item[1], item[2], 1)];
end

# set up the problem data structure
probData = DDUdata(fData.IDList, brList, NCP, J, T, D, Pub, Plb, Q, Wbar, r, θdiff, b, g);

# run the column generation algorithm, linear master problem and cut generation separation
xhat, uhat, UBList, LBList, timeList = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-4, "linear", "piecewise", 1);
xhat, uhat, UBList, LBList, timeList = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-4, "linear", "piecewise");
save("./result/$(caseList[ci])_linear_piecewise.jld", "xhat", xhat, "uhat", uhat, "UBList", UBList, "LBList", LBList, "timeList", timeList);

xhat_biM_lS, uhat_biM_lS, UBList_biM_lS, LBList_biM_lS, timeList_biM_lS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-4, "bilinear", "piecewise", 1);
xhat_biM_lS, uhat_biM_lS, UBList_biM_lS, LBList_biM_lS, timeList_biM_lS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-4, "bilinear", "piecewise");
save("./result/$(caseList[ci])_bilinear_piecewise.jld", "xhat", xhat_biM_lS, "uhat", uhat_biM_lS, "UBList", UBList_biM_lS, "LBList", LBList_biM_lS, "timeList", timeList_biM_lS);

xhat_loS, uhat_loS, UBList_loS, LBList_loS, timeList_loS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-4, "linear", "bilinear_local", 1);
xhat_loS, uhat_loS, UBList_loS, LBList_loS, timeList_loS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-4, "linear", "bilinear_local");
save("./result/$(caseList[ci])_linear_local.jld", "xhat", xhat_loS, "uhat", uhat_loS, "UBList", UBList_loS, "LBList", LBList_loS, "timeList", timeList_loS);

xhat_biM_loS, uhat_biM_loS, UBList_biM_loS, LBList_biM_loS, timeList_biM_loS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-4, "bilinear", "bilinear_local", 1);
xhat_biM_loS, uhat_biM_loS, UBList_biM_loS, LBList_biM_loS, timeList_biM_loS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-4, "bilinear", "bilinear_local");
save("./result/$(caseList[ci])_bilinear_local.jld", "xhat", xhat_biM_loS, "uhat", uhat_biM_loS, "UBList", UBList_biM_loS, "LBList", LBList_biM_loS, "timeList", timeList_biM_loS);

xhat_biS, uhat_biS, UBList_biS, LBList_biS, timeList_biS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-4, "linear", "bilinear_gurobi", 1);
xhat_biS, uhat_biS, UBList_biS, LBList_biS, timeList_biS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-4, "linear", "bilinear_gurobi");
save("./result/$(caseList[ci])_linear_bilinear.jld", "xhat", xhat_biS, "uhat", uhat_biS, "UBList", UBList_biS, "LBList", LBList_biS, "timeList", timeList_biS);

xhat_biM_biS, uhat_biM_biS, UBList_biM_biS, LBList_biM_biS, timeList_biM_biS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-4, "bilinear", "bilinear_gurobi", 1);
xhat_biM_biS, uhat_biM_biS, UBList_biM_biS, LBList_biM_biS, timeList_biM_biS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-4, "bilinear", "bilinear_gurobi");
save("./result/$(caseList[ci])_bilinear_bilinear.jld", "xhat", xhat_biM_biS, "uhat", uhat_biM_biS, "UBList", UBList_biM_biS, "LBList", LBList_biM_biS, "timeList", timeList_biM_biS);
