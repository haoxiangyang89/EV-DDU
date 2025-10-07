# analyze the convergence results and plot them
using Plots;

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
J = 1:6;
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
probData = DDUdata(fData.IDList, brList, NCP, J, T, D, Pub, Plb, Q, Wbar, r, θdiff, b, g, cq);

linear_piecewise_data = load("./result/$(caseList[ci])_linear_piecewise.jld");
bilinear_piecewise_data = load("./result/$(caseList[ci])_bilinear_piecewise.jld");
linear_local_data = load("./result/$(caseList[ci])_linear_local.jld");
bilinear_local_data = load("./result/$(caseList[ci])_bilinear_local.jld");
linear_bilinear_data = load("./result/$(caseList[ci])_linear_bilinear.jld");
bilinear_bilinear_data = load("./result/$(caseList[ci])_bilinear_bilinear.jld");

# plot the convergence progress
linear_piecewise_plot_x = [];
for i in 2:length(linear_piecewise_data["timeList"])
    push!(linear_piecewise_plot_x, sum(linear_piecewise_data["timeList"][1:i]));
end
p = plot(linear_piecewise_plot_x, linear_piecewise_data["LBList"], label="linear_piecewise", 
    linewidth=2, xscale=:log10, xlims=(1, 3600), y_lims = (0,400), xticks=[1,10^1,10^2,10^3,10^3.6], legend=:topleft);

linear_local_plot_x = [];
for i in 2:length(linear_local_data["timeList"])
    push!(linear_local_plot_x, sum(linear_local_data["timeList"][1:i]));
end
plot!(p, linear_local_plot_x, linear_local_data["LBList"], label="linear_local", 
    linewidth=2, xscale=:log10);

linear_bilinear_plot_x = [];
for i in 2:length(linear_bilinear_data["timeList"])
    push!(linear_bilinear_plot_x, sum(linear_bilinear_data["timeList"][1:i]));
end
plot!(p, linear_bilinear_plot_x, linear_bilinear_data["LBList"], label="linear_bilinear", 
    linewidth=2, xscale=:log10);

display(p);

#--------------------------------------------------------------------------------------
# test the solution quality from the linear_piecewise method
xhat = linear_piecewise_data["xhat"];
uhat = linear_piecewise_data["uhat"];
UBhat = linear_piecewise_data["UBList"][length(linear_piecewise_data["UBList"])];
UB_est_sp = build_separation_bilinear(probData, ρ, Delta, xhat, uhat);
set_attribute(UB_est_sp, "TimeLimit", 3600);
optimize!(UB_est_sp);
UB_est = objective_bound(UB_est_sp);
UB_est += sum(c[i] * (xhat[i] - x0[i]) + f[i] * (uhat[i] - u0[i]) for i in NCP);
q_worst = Dict();
for i in probData.NCP
    for j in J
        for t in 1:T
            q_worst[i,j,t] = value(UB_est_sp[:q][i,j,t]);
        end
    end
end
save("./result/$(caseList[ci])_linear_piecewise_solution.jld", "xhat", xhat, "uhat", uhat, "UBhat", UBhat, "UB_est", UB_est, "q_worst", q_worst);
# solve the subproblem with the worst-case q
hsp = solve_h(probData, xhat, q_worst);
optimize!(hsp);
# analyze solution to model (3)
gen_cost = value(hsp[:generation_cost]);
unsat_cost = value(hsp[:unsat_cost]);
cq_cost = sum(sum(probData.cq[i,j] * q_worst[i,j,t] for i in probData.NCP) for j in probData.J for t in 1:probData.T);
investment_cost = sum(c[i] * (xhat[i] - x0[i]) + f[i] * (uhat[i] - u0[i]) for i in probData.NCP);

# generate the demand distribution
dhat = value.(hsp[:d]);
p = plot(1:24, [dhat[2,i] for i in 1:24], label="Bus 2", linewidth=2, xlims=(1, 24), y_lims = (0,2), legend=:topleft);
plot!(p, 1:24, [dhat[9,i] for i in 1:24], label="Bus 9", linewidth=2);
plot!(p, 1:24, [dhat[12,i] for i in 1:24], label="Bus 12", linewidth=2);
display(p);

uzero = Dict();
xzero = Dict();
for i in NCP
    uzero[i] = u0[i];
    xzero[i] = x0[i];
end
hsp_zero = solve_h(probData, xzero, q_worst);
optimize!(hsp_zero);
gen_cost_zero = value(hsp_zero[:generation_cost]);
unsat_cost_zero = value(hsp_zero[:unsat_cost]);
cq_cost_zero = sum(sum(probData.cq[i,j] * q_worst[i,j,t] for i in probData.NCP) for j in probData.J for t in 1:probData.T);