# test the generation + line expansion
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
probData = DDUdata(fData.IDList, brList, NCP, J, T, D, Pub, Plb, Q, Wbar, r, θdiff, b, g, cq);

#------------------------------------------------------------------------------
# test if we further limit the line 7-2, will it result in building a charger at node 7
xhat_rateWbar = Dict();
uhat_rateWbar = Dict();
UBList_rateWbar = Dict();
LBList_rateWbar = Dict();
timeList_rateWbar = Dict();
q_worst_rateWbar = Dict();
gen_cost_rateWbar = Dict();
unsat_cost_rateWbar = Dict();
cq_cost_rateWbar = Dict();
investment_cost_rateWbar = Dict();
UB_est_rateWbar = Dict();
qi_rateWbar = Dict();
P_rateWbar = Dict();
W_rateWbar = Dict();

rateRange = 0.5:-0.1:0.1;
for rateind in eachindex(rateRange)
    rateWbar = rateRange[rateind];
    Wbar[(7,2)] = fData.rateA[(7,2,1)] * rateWbar;
    # set up the problem data structure
    probData = DDUdata(fData.IDList, brList, NCP, J, T, D, Pub, Plb, Q, Wbar, r, θdiff, b, g, cq);
    xhat_rateWbar[rateind], uhat_rateWbar[rateind], UBList_rateWbar[rateind], LBList_rateWbar[rateind], timeList_rateWbar[rateind] = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "bilinear", "bilinear_decomp");
    UB_est_sp, q_worst_rateWbar[rateind] = build_separation_bilinear_decomp(probData, ρ, Delta, xhat_rateWbar[rateind], uhat_rateWbar[rateind], 1e3, 1);
    optimize!(UB_est_sp);
    UB_est = objective_value(UB_est_sp) + sum(c[i] * (xhat_rateWbar[rateind][i] - x0[i]) + f[i] * (uhat_rateWbar[rateind][i] - u0[i]) for i in NCP);
    hsp = solve_h(probData, xhat, q_worst_rateWbar[rateind]);
    optimize!(hsp);
    # analyze solution to model (3)
    gen_cost_rateWbar[rateind] = value(hsp[:generation_cost]);
    unsat_cost_rateWbar[rateind] = value(hsp[:unsat_cost]);
    cq_cost_rateWbar[rateind] = sum(sum(probData.cq[i,j] * q_worst_rateWbar[rateind][i,j,t] for i in probData.NCP) for j in probData.J for t in 1:probData.T);
    investment_cost_rateWbar[rateind] = sum(c[i] * (xhat_rateWbar[rateind][i] - x0[i]) + f[i] * (uhat_rateWbar[rateind][i] - u0[i]) for i in probData.NCP);
    P_rateWbar[rateind] = value.(hsp[:P]);
    W_rateWbar[rateind] = value.(hsp[:W]);
    save("test_rateWbar.jld", "UB_est", UB_est, "xhat_rateWbar", xhat_rateWbar, "uhat_rateWbar", uhat_rateWbar, "UBList_rateWbar", UBList_rateWbar,
             "LBList_rateWbar", LBList_rateWbar, "timeList_rateWbar", timeList_rateWbar,
            "q_worst_rateWbar", q_worst_rateWbar, "gen_cost_rateWbar", gen_cost_rateWbar, "unsat_cost_rateWbar", unsat_cost_rateWbar,
            "cq_cost_rateWbar", cq_cost_rateWbar, "investment_cost_rateWbar", investment_cost_rateWbar,
            "P_rateWbar", P_rateWbar, "W_rateWbar", W_rateWbar);
end

#------------------------------------------------------------------------------
# test if we further decrease the investment cost, will we build a charger anywhere?
Wbar[(7,2)] = fData.rateA[(7,2,1)];
probData = DDUdata(fData.IDList, brList, NCP, J, T, D, Pub, Plb, Q, Wbar, r, θdiff, b, g, cq);
xhat_discount = Dict();
uhat_discount = Dict();
UBList_discount = Dict();
LBList_discount = Dict();
timeList_discount = Dict();
q_worst_discount = Dict();
gen_cost_discount = Dict();
unsat_cost_discount = Dict();
cq_cost_discount = Dict();
investment_cost_discount = Dict();
P_discount = Dict();
W_discount = Dict();

for discount_rate in 1:5
    new_c = Dict();
    new_f = Dict();
    for i in NCP
        new_c[i] = c[i] * (1 - discount_rate * 0.1);
        new_f[i] = f[i] * (1 - discount_rate * 0.1);
    end
    xhat_discount[discount_rate], uhat_discount[discount_rate], UBList_discount[discount_rate], LBList_discount[discount_rate], timeList_discount[discount_rate] = 
        colGen(probData, new_f, new_c, xbar, x0, u0, ρ, Delta, 1e-3, "bilinear", "bilinear_decomp");
    UB_est_sp, q_worst_discount[discount_rate] = build_separation_bilinear_decomp(probData, ρ, Delta, xhat_discount[discount_rate], uhat_discount[discount_rate], 1e3, 1);
    set_attribute(UB_est_sp, "TimeLimit", 3600);
    optimize!(UB_est_sp);
    UB_est = objective_value(UB_est_sp) + sum(c[i] * (xhat_discount[discount_rate][i] - x0[i]) + f[i] * (uhat_discount[discount_rate][i] - u0[i]) for i in NCP);
    hsp = solve_h(probData, xhat_discount[discount_rate], q_worst_discount[discount_rate]);
    optimize!(hsp);
    # analyze solution to model (3)
    gen_cost_discount[discount_rate] = value(hsp[:generation_cost]);
    unsat_cost_discount[discount_rate] = value(hsp[:unsat_cost]);
    cq_cost_discount[discount_rate] = sum(sum(probData.cq[i,j] * q_worst_discount[discount_rate][i,j,t] for i in probData.NCP) for j in probData.J for t in 1:probData.T);
    investment_cost_discount[discount_rate] = sum(c[i] * (xhat_discount[discount_rate][i] - x0[i]) + f[i] * (uhat_discount[discount_rate][i] - u0[i]) for i in probData.NCP);
    P_discount[discount_rate] = value.(hsp[:P]);
    W_discount[discount_rate] = value.(hsp[:W]);
    save("test_discount_rate.jld", "UB_est", UB_est, "xhat_discount", xhat_discount, "uhat_discount", uhat_discount, "UBList_discount", UBList_discount,
             "LBList_discount", LBList_discount, "timeList_discount", timeList_discount,
            "q_worst_discount", q_worst_discount, "gen_cost_discount", gen_cost_discount, "unsat_cost_discount", unsat_cost_discount,
            "cq_cost_discount", cq_cost_discount, "investment_cost_discount", investment_cost_discount,
            "P_discount", P_discount, "W_discount", W_discount);
end

#------------------------------------------------------------------------------
# test if we add a low cost generator somewhere, will it make a difference?
nonGen = [2,3,4,6,8,9,10,11];
for i in nonGen
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

    Pub[i] = 2.0;
    g[i] = 1.0;
    probData = DDUdata(fData.IDList, brList, NCP, J, T, D, Pub, Plb, Q, Wbar, r, θdiff, b, g, cq);
    xhat, uhat, UBList, LBList, timeList = colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "bilinear", "bilinear_decomp");
    UB_est_sp, q_worst = build_separation_bilinear(probData, ρ, Delta, xhat, uhat, 1e3, 1);
    set_attribute(UB_est_sp, "TimeLimit", 3600);
    optimize!(UB_est_sp);
    UB_est = objective_value(UB_est_sp) + sum(c[i] * (xhat[i] - x0[i]) + f[i] * (uhat[i] - u0[i]) for i in NCP);
    hsp = solve_h(probData, xhat, q_worst);
    optimize!(hsp);
    # analyze solution to model (3)
    gen_cost = value(hsp[:generation_cost]);
    unsat_cost = value(hsp[:unsat_cost]);
    cq_cost = sum(sum(probData.cq[i,j] * q_worst[i,j,t] for i in probData.NCP) for j in probData.J for t in 1:probData.T);
    investment_cost = sum(c[i] * (xhat[i] - x0[i]) + f[i] * (uhat[i] - u0[i]) for i in probData.NCP);
end

#------------------------------------------------------------------------------
# test if we increase the load to a certain percentage, will we set up more chargers? 
Q_perc_list = [1.1, 1.2, 1.3, 1.4, 1.5];
expansion_list = [1.5, 2.0, 2.5];
xhat_Q_expansion = Dict();
uhat_Q_expansion = Dict();
UBList_Q_expansion = Dict();
LBList_Q_expansion = Dict();
timeList_Q_expansion = Dict();
q_worst_Q_expansion = Dict();
gen_cost_Q_expansion = Dict();
unsat_cost_Q_expansion = Dict();
cq_cost_Q_expansion = Dict();
investment_cost_Q_expansion = Dict();
P_Q_expansion = Dict();
W_Q_expansion = Dict();

for Q_perc_ind in eachindex(Q_perc_list)
    for j in J
        for t in 1:T
            Q[j,t] = df_Q[t,"x$j"] * Q_perc_ind;
        end
    end
    Pub = probData.Pub;
    # test if we expand generation at bus 7 and the capacity of line 7-2
    for expansion_ind in eachindex(expansion_list)
        Pub[7] = sum(fData.Pmax[g] for g in fData.LocRev[7]) * expansion_list[expansion_ind];
        Wbar[(7,2)] = probData.Wbar[(7,2)] * expansion_list[expansion_ind];
        probData = DDUdata(fData.IDList, brList, NCP, J, T, D, Pub, Plb, Q, Wbar, r, θdiff, b, g, cq);
        xhat_Q_expansion[Q_perc_ind,expansion_ind], uhat_Q_expansion[Q_perc_ind,expansion_ind], UBList_Q_expansion[Q_perc_ind,expansion_ind], 
            LBList_Q_expansion[Q_perc_ind,expansion_ind], timeList_Q_expansion[Q_perc_ind,expansion_ind] = 
            colGen(probData, f, c, xbar, x0, u0, ρ, Delta, 1e-3, "linear", "bilinear_decomp");
        UB_est_sp, q_worst_Q_expansion[Q_perc_ind,expansion_ind] = build_separation_bilinear_decomp(probData, ρ, Delta, xhat_Q_expansion[Q_perc_ind,expansion_ind], 
            uhat_Q_expansion[Q_perc_ind,expansion_ind], 1e3, 1);
            optimize!(UB_est_sp);
        UB_est = objective_value(UB_est_sp) + sum(c[i] * (xhat_Q_expansion[Q_perc_ind,expansion_ind][i] - x0[i]) + f[i] * (uhat_Q_expansion[Q_perc_ind,expansion_ind][i] - u0[i]) for i in NCP);
        hsp = solve_h(probData, xhat_Q_expansion[Q_perc_ind,expansion_ind], q_worst_Q_expansion[Q_perc_ind,expansion_ind]);
        optimize!(hsp);
        # analyze solution to model (3)
        gen_cost_Q_expansion[Q_perc_ind,expansion_ind] = value(hsp[:generation_cost]);
        unsat_cost_Q_expansion[Q_perc_ind,expansion_ind] = value(hsp[:unsat_cost]);
        cq_cost_Q_expansion[Q_perc_ind,expansion_ind] = sum(sum(probData.cq[i,j] * q_worst_Q_expansion[Q_perc_ind,expansion_ind][i,j,t] for i in probData.NCP) for j in probData.J for t in 1:probData.T);
        investment_cost_Q_expansion[Q_perc_ind,expansion_ind] = sum(c[i] * (xhat_Q_expansion[Q_perc_ind,expansion_ind][i] - x0[i]) + f[i] * (uhat_Q_expansion[Q_perc_ind,expansion_ind][i] - u0[i]) for i in probData.NCP);
        P_Q_expansion[Q_perc_ind,expansion_ind] = value.(hsp[:P]);
        W_Q_expansion[Q_perc_ind,expansion_ind] = value.(hsp[:W]);
        save("test_Q_expansion_linear.jld", "UB_est", UB_est, "xhat_Q_expansion", xhat_Q_expansion, "uhat_Q_expansion", uhat_Q_expansion, 
             "UBList_Q_expansion", UBList_Q_expansion, "LBList_Q_expansion", LBList_Q_expansion, "timeList_Q_expansion", timeList_Q_expansion,
            "q_worst_Q_expansion", q_worst_Q_expansion, "gen_cost_Q_expansion", gen_cost_Q_expansion, "unsat_cost_Q_expansion", unsat_cost_Q_expansion,
            "cq_cost_Q_expansion", cq_cost_Q_expansion, "investment_cost_Q_expansion", investment_cost_Q_expansion,
            "P_Q_expansion", P_Q_expansion, "W_Q_expansion", W_Q_expansion);
    end

    # test if we install more RES capacity
end