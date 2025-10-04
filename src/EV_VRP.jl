# script to solve the stochastic EV VRP in Sweden

using JuMP, Gurobi, Combinatorics, LinearAlgebra, JLD, HDF5, Distributions, DelimitedFiles;

const GUROBI_ENV = Gurobi.Env();

function extensive_form(data, scen_data)
    # construct the extensive formulation
    ext_mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)));
    # first stage variables & constraints
    @variable(ext_mp, x[k in data.K] >= 0, Int);
    @expression(ext_mp, deploy_cost, sum(data.c[k] * x[k] for k in data.K));
    @constraint(ext_mp, budget, sum(data.c[k] * x[k] for k in data.K) <= data.B);
    
    # second stage variables & constraints
    @variable(ext_mp, y[a in data.A, k in data.K, s in scen_data.Omega], Bin);
    @variable(ext_mp, z[i in data.I, k in data.K, s in scen_data.Omega], Bin);
end