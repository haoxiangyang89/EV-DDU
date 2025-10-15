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
J = 1:10;
df_Q = CSV.read("./data/$(caseList[ci])/case$(caseList[ci])_Q_new.csv", DataFrame);
df_r = CSV.read("./data/$(caseList[ci])/case$(caseList[ci])_rdis_new.csv", DataFrame);

# set up the charger data
#NCP = [2,3,5,9,11,12,13];
NCP = 1:13;
c = Dict(1 => 1.0, 2 => 1.5, 3 => 1.5, 4 => 1.0, 5 => 1.5, 6 => 1.0, 7 => 1.0, 
         8 => 1.0, 9 => 1.4, 10 => 1.0, 11 => 1.0, 12 => 1.2, 13 => 1.8);
f = Dict(1 => 10.0, 2 => 18.0, 3 => 25.0, 4 => 10.0, 5 => 23.0, 6 => 10.0, 7 => 10.0, 
         8 => 10.0, 9 => 16.0, 10 => 10.0, 11 => 10.0, 12 => 15.0, 13 => 26.0);
xbar = Dict();
x0 = Dict();
u0 = Dict();
for i in NCP
    xbar[i] = 15.0;
    x0[i] = 0.0;
    u0[i] = 0;
end
# f[2] = 10.0;
# f[9] = 12.0;
# f[3] = 25.0;
# f[5] = 23.0;
# f[11] = 10.0;
# f[12] = 15.0;
# f[13] = 26.0;
# c[2] = 1.0;
# c[9] = 1.0;
# c[3] = 1.5;
# c[11] = 1.0;
# c[5] = 1.5;
# c[12] = 1.2;
# c[13] = 1.8;

for i in NCP
    f[i] = f[i] * 30;
    c[i] = c[i] * 100;
end

current_charger = [2,9,12];
for i in current_charger
    x0[i] = 1.0;
    u0[i] = 1;
end
ρ = 50;

# calculate r_ij and set Delta
r = Dict();
cq = Dict();
for j in J
    for i in NCP
        r[i,j] = 100/df_r[i,"x$j"];
        cq[i,j] = 0.1/r[i,j];
        #cq[i,j] = 0.0;
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
probData = DDUdata(fData.IDList, brList, NCP, J, T, D, Pub, Plb, Q, Wbar, r, θdiff, b, g, cq);

# run the column generation algorithm, linear master problem and cut generation separation
xhat, uhat, UBList, LBList, timeList = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "linear", "piecewise", 1);
xhat, uhat, UBList, LBList, timeList = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "linear", "piecewise");
save("./result/$(caseList[ci])_linear_piecewise.jld", "xhat", xhat, "uhat", uhat, "UBList", UBList, "LBList", LBList, "timeList", timeList);

# xhat_biM_lS, uhat_biM_lS, UBList_biM_lS, LBList_biM_lS, timeList_biM_lS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "bilinear", "piecewise", 1);
# xhat_biM_lS, uhat_biM_lS, UBList_biM_lS, LBList_biM_lS, timeList_biM_lS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "bilinear", "piecewise");
# save("./result/$(caseList[ci])_bilinear_piecewise.jld", "xhat", xhat_biM_lS, "uhat", uhat_biM_lS, "UBList", UBList_biM_lS, "LBList", LBList_biM_lS, "timeList", timeList_biM_lS);

xhat_loS, uhat_loS, UBList_loS, LBList_loS, timeList_loS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "linear", "bilinear_local", 1);
xhat_loS, uhat_loS, UBList_loS, LBList_loS, timeList_loS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "linear", "bilinear_local");
save("./result/$(caseList[ci])_linear_local.jld", "xhat", xhat_loS, "uhat", uhat_loS, "UBList", UBList_loS, "LBList", LBList_loS, "timeList", timeList_loS);

# xhat_biM_loS, uhat_biM_loS, UBList_biM_loS, LBList_biM_loS, timeList_biM_loS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "bilinear", "bilinear_local", 1);
# xhat_biM_loS, uhat_biM_loS, UBList_biM_loS, LBList_biM_loS, timeList_biM_loS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "bilinear", "bilinear_local");
# save("./result/$(caseList[ci])_bilinear_local.jld", "xhat", xhat_biM_loS, "uhat", uhat_biM_loS, "UBList", UBList_biM_loS, "LBList", LBList_biM_loS, "timeList", timeList_biM_loS);

xhat_biS, uhat_biS, UBList_biS, LBList_biS, timeList_biS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "linear", "bilinear_gurobi", 1);
xhat_biS, uhat_biS, UBList_biS, LBList_biS, timeList_biS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "linear", "bilinear_gurobi");
save("./result/$(caseList[ci])_linear_bilinear.jld", "xhat", xhat_biS, "uhat", uhat_biS, "UBList", UBList_biS, "LBList", LBList_biS, "timeList", timeList_biS);

# xhat_biM_biS, uhat_biM_biS, UBList_biM_biS, LBList_biM_biS, timeList_biM_biS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "bilinear", "bilinear_gurobi", 1);
# xhat_biM_biS, uhat_biM_biS, UBList_biM_biS, LBList_biM_biS, timeList_biM_biS = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "bilinear", "bilinear_gurobi");
# save("./result/$(caseList[ci])_bilinear_bilinear.jld", "xhat", xhat_biM_biS, "uhat", uhat_biM_biS, "UBList", UBList_biM_biS, "LBList", LBList_biM_biS, "timeList", timeList_biM_biS);

xhat_bid, uhat_bid, UBList_bid, LBList_bid, timeList_bid = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "linear", "bilinear_decomp", 1);
xhat_bid, uhat_bid, UBList_bid, LBList_bid, timeList_bid = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "linear", "bilinear_decomp");
save("./result/$(caseList[ci])_linear_decomp.jld", "xhat", xhat_bid, "uhat", uhat_bid, "UBList", UBList_bid, "LBList", LBList_bid, "timeList", timeList_bid);

xhat_biM_bid, uhat_biM_bid, UBList_biM_bid, LBList_biM_bid, timeList_biM_bid = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "bilinear", "bilinear_decomp", 1);
xhat_biM_bid, uhat_biM_bid, UBList_biM_bid, LBList_biM_bid, timeList_biM_bid = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "bilinear", "bilinear_decomp");
save("./result/$(caseList[ci])_bilinear_decomp.jld", "xhat", xhat_biM_bid, "uhat", uhat_biM_bid, "UBList", UBList_biM_bid, "LBList", LBList_biM_bid, "timeList", timeList_biM_bid);