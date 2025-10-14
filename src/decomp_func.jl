# main script to run column generation algorithm for EV_DDU
function build_masterproblem(probData, f, c, xbar, x0, u0, Delta, K, πList, M = 1e4)
    NCP = probData.NCP;
    J = probData.J;
    T = 1:probData.T;
    λbar = 100 * sum(values(probData.Q));
    # λbar = M;
    # build the relaxed master problem here
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0, "Threads" => 10));
    @variable(mp, x[i in NCP] >= x0[i]); # charger capacity
    @variable(mp, u[i in NCP], Bin); # charger installation indicator
    @variable(mp, V >= 0); # worst-case second-stage cost
    @objective(mp, Min, sum(c[i] * (x[i] - x0[i]) + f[i] * (u[i] - u0[i]) for i in NCP) + V);
    @constraint(mp, capacity_ub[i in NCP], x[i] <= xbar[i] * u[i]);
    @constraint(mp, installation_previous[i in NCP], u[i] >= u0[i]);

    # initialize with columns
    if K != []
        # add z variables and corresponding λ variables
        @variable(mp, z[k in K], container = SparseAxisArray);
        @variable(mp, λQ[k in K, j in J, t in T], container = SparseAxisArray);
        @variable(mp, λu[k in K, i in NCP, j in J, t in T], container = SparseAxisArray, lower_bound = 0);
        @variable(mp, uλu[k in K, i in NCP, j in J, t in T], container = SparseAxisArray, lower_bound = 0, upper_bound = λbar); # bilinear term
        @variable(mp, λqu[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, lower_bound = 0);
        @variable(mp, uλqu1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, lower_bound = 0, upper_bound = λbar); # bilinear term
        @variable(mp, uλqu2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, lower_bound = 0, upper_bound = λbar); # bilinear term
        @variable(mp, λql[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, upper_bound = 0);
        @variable(mp, uλql1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, upper_bound = 0, lower_bound = -λbar); # bilinear term
        @variable(mp, uλql2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, upper_bound = 0, lower_bound = -λbar); # bilinear term
        
        for k in K
            πitem = πList[k];
            # add the V constraint
            @constraint(mp, V_constr[k in K], V >= sum((sum(D[i,t] * πitem["P"][i,t] + probData.Pub[i] * πitem["Pu"][i,t] + probData.Plb[i] * πitem["Pl"][i,t] for i in probData.IDList) + 
                                                    sum(πitem["d"][i,t] * x[i] for i in NCP)) + 
                                                    sum((πitem["Wu"][br[1],br[2],t] - πitem["Wl"][br[1],br[2],t]) * probData.Wbar[br[1],br[2]] +
                                                        (πitem["θu"][br[1],br[2],t] - πitem["θl"][br[1],br[2],t]) * probData.θdiff[br[1],br[2]] for br in probData.brList) for t in T) + 
                                                    z[k], container = SparseAxisArray);
            # add the z constraint
            # @constraint(mp, z_constr[k in K], z[k] >= sum(sum(λQ[k,j,t] * probData.Q[j,t] + sum(λu[k,i,j,t] * probData.Q[j,t] * u[i] for i in NCP) + 
            #                 sum(M * (2 - u[br[1]] - u[br[2]]) * (λqu[k,br[1],br[2],j,t] - λql[k,br[1],br[2],j,t]) for br in probData.brList) for j in J) for t in T),
            #                 container = SparseAxisArray);
            @constraint(mp, z_constr[k in K], z[k] >= sum(sum(λQ[k,j,t] * probData.Q[j,t] + sum(uλu[k,i,j,t] * probData.Q[j,t] for i in NCP) + 
                            sum(M * (2 * (λqu[k,br[1],br[2],j,t] - λql[k,br[1],br[2],j,t]) - uλqu1[k,br[1],br[2],j,t] - uλqu2[k,br[1],br[2],j,t] + uλql1[k,br[1],br[2],j,t] + uλql2[k,br[1],br[2],j,t])
                                for br in probData.brList) for j in J) for t in T), container = SparseAxisArray);
            # add the λ constraints
            @constraint(mp, λ_constr[k in K, i in NCP, j in J, t in T], λQ[k,j,t] + λu[k,i,j,t] + sum(λqu[k,br[1],br[2],j,t] + λql[k,br[1],br[2],j,t] for br in probData.brList if br[1] == i) - 
                            sum((probData.r[br[2],j]/probData.r[br[1],j]) * (1 + Delta) * λqu[k,br[1],br[2],j,t] + 
                                (probData.r[br[2],j]/probData.r[br[1],j]) * (1 - Delta) * λql[k,br[1],br[2],j,t] for br in probData.brList if br[2] == i) >= πitem["q"][i,t] + probData.cq[i,j],
                            container = SparseAxisArray);

            # add the λ-u logic constraints
            @constraint(mp, λqu_logic_constr1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], λqu[k,i,ip,j,t] <= M * u[i], container = SparseAxisArray);
            @constraint(mp, λqu_logic_constr2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], λqu[k,i,ip,j,t] <= M * u[ip], container = SparseAxisArray);
            @constraint(mp, λql_logic_constr1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], λql[k,i,ip,j,t] >= -M * u[i], container = SparseAxisArray);
            @constraint(mp, λql_logic_constr2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], λql[k,i,ip,j,t] >= -M * u[ip], container = SparseAxisArray);
            
            # add McCormick envelopes for bilinear terms
            @constraint(mp, uλu_constr1[k in K, i in NCP, j in J, t in T], uλu[k,i,j,t] <= λbar * u[i], container = SparseAxisArray);
            @constraint(mp, uλu_constr2[k in K, i in NCP, j in J, t in T], uλu[k,i,j,t] <= λu[k,i,j,t], container = SparseAxisArray);
            @constraint(mp, uλu_constr3[k in K, i in NCP, j in J, t in T], uλu[k,i,j,t] >= λu[k,i,j,t] - λbar * (1 - u[i]), container = SparseAxisArray);

            @constraint(mp, uλqu1_constr1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], uλqu1[k,i,ip,j,t] <= λbar * u[i], container = SparseAxisArray);
            @constraint(mp, uλqu1_constr2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], uλqu1[k,i,ip,j,t] <= λqu[k,i,ip,j,t], container = SparseAxisArray);
            @constraint(mp, uλqu1_constr3[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], uλqu1[k,i,ip,j,t] >= λqu[k,i,ip,j,t] - λbar * (1 - u[i]), container = SparseAxisArray);
            @constraint(mp, uλqu2_constr1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], uλqu2[k,i,ip,j,t] <= λbar * u[ip], container = SparseAxisArray);
            @constraint(mp, uλqu2_constr2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], uλqu2[k,i,ip,j,t] <= λqu[k,i,ip,j,t], container = SparseAxisArray);
            @constraint(mp, uλqu2_constr3[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], uλqu2[k,i,ip,j,t] >= λqu[k,i,ip,j,t] - λbar * (1 - u[ip]), container = SparseAxisArray);

            @constraint(mp, uλql1_constr1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], uλql1[k,i,ip,j,t] >= -λbar * u[i], container = SparseAxisArray);
            @constraint(mp, uλql1_constr2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], uλql1[k,i,ip,j,t] >= λql[k,i,ip,j,t], container = SparseAxisArray);
            @constraint(mp, uλql1_constr3[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], uλql1[k,i,ip,j,t] <= λql[k,i,ip,j,t] + λbar * (1 - u[i]), container = SparseAxisArray);
            @constraint(mp, uλql2_constr1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], uλql2[k,i,ip,j,t] >= -λbar * u[ip], container = SparseAxisArray);
            @constraint(mp, uλql2_constr2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], uλql2[k,i,ip,j,t] >= λql[k,i,ip,j,t], container = SparseAxisArray);
            @constraint(mp, uλql2_constr3[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], uλql2[k,i,ip,j,t] <= λql[k,i,ip,j,t] + λbar * (1 - u[ip]), container = SparseAxisArray);
        end
    end
    return mp;
end

function build_masterproblem_bilinear(probData, f, c, xbar, x0, u0, Delta, K, πList, M = 1e4)
    NCP = probData.NCP;
    J = probData.J;
    T = 1:probData.T;
    # build the relaxed master problem here
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0, "Threads" => 10));
    @variable(mp, x[i in NCP] >= x0[i]); # charger capacity
    @variable(mp, u[i in NCP], Bin); # charger installation indicator
    @variable(mp, V >= 0); # worst-case second-stage cost
    @objective(mp, Min, sum(c[i] * (x[i] - x0[i]) + f[i] * (u[i] - u0[i]) for i in NCP) + V);
    @constraint(mp, capacity_ub[i in NCP], x[i] <= xbar[i] * u[i]);
    @constraint(mp, installation_previous[i in NCP], u[i] >= u0[i]);

    # initialize with columns
    if K != []
        # add z variables and corresponding λ variables
        @variable(mp, z[k in K], container = SparseAxisArray);
        @variable(mp, λQ[k in K, j in J, t in T], container = SparseAxisArray);
        @variable(mp, λu[k in K, i in NCP, j in J, t in T], container = SparseAxisArray, lower_bound = 0);
        @variable(mp, λqu[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, lower_bound = 0);
        @variable(mp, λql[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, upper_bound = 0);
        for k in K
            πitem = πList[k];
            # add the V constraint
            @constraint(mp, V_constr[k in K], V >= sum((sum(D[i,t] * πitem["P"][i,t] + probData.Pub[i] * πitem["Pu"][i,t] + probData.Plb[i] * πitem["Pl"][i,t] for i in probData.IDList) + 
                                                    sum(πitem["d"][i,t] * x[i] for i in NCP)) + 
                                                    sum((πitem["Wu"][br[1],br[2],t] - πitem["Wl"][br[1],br[2],t]) * probData.Wbar[br[1],br[2]] +
                                                        (πitem["θu"][br[1],br[2],t] - πitem["θl"][br[1],br[2],t]) * probData.θdiff[br[1],br[2]] for br in probData.brList) for t in T) + 
                                                    z[k], container = SparseAxisArray);
            # add the z constraint
            @constraint(mp, z_constr[k in K], z[k] >= sum(sum(λQ[k,j,t] * probData.Q[j,t] + sum(λu[k,i,j,t] * probData.Q[j,t] * u[i] for i in NCP) + 
                            sum(M * (2 - u[br[1]] - u[br[2]]) * (λqu[k,br[1],br[2],j,t] - λql[k,br[1],br[2],j,t]) for br in probData.brList) for j in J) for t in T),
                            container = SparseAxisArray);
            # add the λ constraints
            @constraint(mp, λ_constr[k in K, i in NCP, j in J, t in T], λQ[k,j,t] + λu[k,i,j,t] + sum(λqu[k,br[1],br[2],j,t] + λql[k,br[1],br[2],j,t] for br in probData.brList if br[1] == i) - 
                            sum((probData.r[br[2],j]/probData.r[br[1],j]) * (1 + Delta) * λqu[k,br[1],br[2],j,t] + 
                                (probData.r[br[2],j]/probData.r[br[1],j]) * (1 - Delta) * λql[k,br[1],br[2],j,t] for br in probData.brList if br[2] == i) >= πitem["q"][i,t] + probData.cq[i,j],
                            container = SparseAxisArray);
        end
    end
    return mp;
end

function addCol(mp, probData, Delta, π_add, K, k, M = 1e4)
    NCP = probData.NCP;
    J = probData.J;
    T = 1:probData.T;
    λbar = 100 * sum(values(probData.Q));
    # λbar = M;

    if k > 1
        # add a new column to the master problem here
        mp[:z][k] = @variable(mp, base_name = "z[$k]");
        for j in J
            for t in T
                mp[:λQ][k,j,t] = @variable(mp, base_name = "λQ[$k,$j,$t]");
                for i in probData.NCP
                    mp[:λu][k,i,j,t] = @variable(mp, lower_bound = 0, base_name = "λu[$k,$i,$j,$t]");
                    mp[:uλu][k,i,j,t] = @variable(mp, lower_bound = 0, upper_bound = λbar, base_name = "uλu[$k,$i,$j,$t]"); # bilinear term
                    for ip in probData.NCP
                        if ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)
                            mp[:λqu][k,i,ip,j,t] = @variable(mp, lower_bound = 0, base_name = "λqu[$k,$i,$ip,$j,$t]");
                            mp[:uλqu1][k,i,ip,j,t] = @variable(mp, lower_bound = 0, upper_bound = λbar, base_name = "uλqu1[$k,$i,$ip,$j,$t]"); # bilinear term
                            mp[:uλqu2][k,i,ip,j,t] = @variable(mp, lower_bound = 0, upper_bound = λbar, base_name = "uλqu2[$k,$i,$ip,$j,$t]"); # bilinear term
                            mp[:λql][k,i,ip,j,t] = @variable(mp, upper_bound = 0, base_name = "λql[$k,$i,$ip,$j,$t]");
                            mp[:uλql1][k,i,ip,j,t] = @variable(mp, upper_bound = 0, lower_bound = -λbar, base_name = "uλql1[$k,$i,$ip,$j,$t]"); # bilinear term
                            mp[:uλql2][k,i,ip,j,t] = @variable(mp, upper_bound = 0, lower_bound = -λbar, base_name = "uλql2[$k,$i,$ip,$j,$t]"); # bilinear term
                        end
                    end
                end
            end
        end
        # add the V constraint
        mp[:V_constr][k] = @constraint(mp, mp[:V] >= sum((sum(D[i,t] * π_add["P"][i,t] + probData.Pub[i] * π_add["Pu"][i,t] + probData.Plb[i] * π_add["Pl"][i,t] for i in probData.IDList) + 
                                                        sum(π_add["d"][i,t] * mp[:x][i] for i in NCP)) + 
                                                        sum((π_add["Wu"][br[1],br[2],t] - π_add["Wl"][br[1],br[2],t]) * probData.Wbar[br[1],br[2]] +
                                                            (π_add["θu"][br[1],br[2],t] - π_add["θl"][br[1],br[2],t]) * probData.θdiff[br[1],br[2]] for br in probData.brList) for t in T) + 
                                                        mp[:z][k]);
        # add the z constraint
        mp[:z_constr][k] = @constraint(mp, mp[:z][k] >= sum(sum(mp[:λQ][k,j,t] * probData.Q[j,t] + sum(mp[:uλu][k,i,j,t] * probData.Q[j,t] for i in NCP) + 
                        sum(M * (2 * (mp[:λqu][k,br[1],br[2],j,t] - mp[:λql][k,br[1],br[2],j,t]) - mp[:uλqu1][k,br[1],br[2],j,t] - mp[:uλqu2][k,br[1],br[2],j,t] + mp[:uλql1][k,br[1],br[2],j,t] + mp[:uλql2][k,br[1],br[2],j,t]) 
                        for br in probData.brList if (br[1] in NCP)&(br[2] in NCP)) for j in J) for t in T));
        # add the λ constraints
        for i in probData.NCP
            for j in J
                for t in T
                    mp[:λ_constr][k,i,j,t] = @constraint(mp, mp[:λQ][k,j,t] + mp[:λu][k,i,j,t] + sum(mp[:λqu][k,br[1],br[2],j,t] + mp[:λql][k,br[1],br[2],j,t] for br in probData.brList if (br[1] in NCP)&(br[2] in NCP)&(br[1] == i)) - 
                                sum((probData.r[br[2],j]/probData.r[br[1],j]) * (1 + Delta) * mp[:λqu][k,br[1],br[2],j,t] + 
                                    (probData.r[br[2],j]/probData.r[br[1],j]) * (1 - Delta) * mp[:λql][k,br[1],br[2],j,t] for br in probData.brList if (br[1] in NCP)&(br[2] in NCP)&(br[2] == i)) >= π_add["q"][i,t] + probData.cq[i,j]);
                    # add McCormick envelopes for bilinear terms
                    mp[:uλu_constr1][k,i,j,t] = @constraint(mp, mp[:uλu][k,i,j,t] <= λbar * mp[:u][i]);
                    mp[:uλu_constr2][k,i,j,t] = @constraint(mp, mp[:uλu][k,i,j,t] <= mp[:λu][k,i,j,t]);
                    mp[:uλu_constr3][k,i,j,t] = @constraint(mp, mp[:uλu][k,i,j,t] >= mp[:λu][k,i,j,t] - λbar * (1 - mp[:u][i]));

                    for ip in probData.NCP
                        # add McCormick envelopes for bilinear terms
                        if ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)
                            # add the λ-u logic constraints
                            mp[:λqu_logic_constr1][k,i,ip,j,t] = @constraint(mp, mp[:λqu][k,i,ip,j,t] <= M * mp[:u][i]);
                            mp[:λqu_logic_constr2][k,i,ip,j,t] = @constraint(mp, mp[:λqu][k,i,ip,j,t] <= M * mp[:u][ip]);
                            mp[:λql_logic_constr1][k,i,ip,j,t] = @constraint(mp, mp[:λql][k,i,ip,j,t] >= -M * mp[:u][i]);
                            mp[:λql_logic_constr2][k,i,ip,j,t] = @constraint(mp, mp[:λql][k,i,ip,j,t] >= -M * mp[:u][ip]);

                            mp[:uλqu1_constr1][k,i,ip,j,t] = @constraint(mp, mp[:uλqu1][k,i,ip,j,t] <= λbar * mp[:u][i]);
                            mp[:uλqu1_constr2][k,i,ip,j,t] = @constraint(mp, mp[:uλqu1][k,i,ip,j,t] <= mp[:λqu][k,i,ip,j,t]);
                            mp[:uλqu1_constr3][k,i,ip,j,t] = @constraint(mp, mp[:uλqu1][k,i,ip,j,t] >= mp[:λqu][k,i,ip,j,t] - λbar * (1 - mp[:u][i]));
                            mp[:uλqu2_constr1][k,i,ip,j,t] = @constraint(mp, mp[:uλqu2][k,i,ip,j,t] <= λbar * mp[:u][ip]);
                            mp[:uλqu2_constr2][k,i,ip,j,t] = @constraint(mp, mp[:uλqu2][k,i,ip,j,t] <= mp[:λqu][k,i,ip,j,t]);
                            mp[:uλqu2_constr3][k,i,ip,j,t] = @constraint(mp, mp[:uλqu2][k,i,ip,j,t] >= mp[:λqu][k,i,ip,j,t] - λbar * (1 - mp[:u][ip]));

                            mp[:uλql1_constr1][k,i,ip,j,t] = @constraint(mp, mp[:uλql1][k,i,ip,j,t] >= -λbar * mp[:u][i]);
                            mp[:uλql1_constr2][k,i,ip,j,t] = @constraint(mp, mp[:uλql1][k,i,ip,j,t] >= mp[:λql][k,i,ip,j,t]);
                            mp[:uλql1_constr3][k,i,ip,j,t] = @constraint(mp, mp[:uλql1][k,i,ip,j,t] <= mp[:λql][k,i,ip,j,t] + λbar * (1 - mp[:u][i]));
                            mp[:uλql2_constr1][k,i,ip,j,t] = @constraint(mp, mp[:uλql2][k,i,ip,j,t] >= -λbar * mp[:u][ip]);
                            mp[:uλql2_constr2][k,i,ip,j,t] = @constraint(mp, mp[:uλql2][k,i,ip,j,t] >= mp[:λql][k,i,ip,j,t]);
                            mp[:uλql2_constr3][k,i,ip,j,t] = @constraint(mp, mp[:uλql2][k,i,ip,j,t] <= mp[:λql][k,i,ip,j,t] + λbar * (1 - mp[:u][ip]));
                        end
                    end
                end
            end
        end
    else
        @variable(mp, z[k in K], container = SparseAxisArray);
        @variable(mp, λQ[k in K, j in J, t in T], container = SparseAxisArray);
        @variable(mp, λu[k in K, i in NCP, j in J, t in T], container = SparseAxisArray, lower_bound = 0);
        @variable(mp, uλu[k in K, i in NCP, j in J, t in T], container = SparseAxisArray, lower_bound = 0, upper_bound = λbar); # bilinear term
        @variable(mp, λqu[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, lower_bound = 0);
        @variable(mp, uλqu1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, lower_bound = 0, upper_bound = λbar); # bilinear term
        @variable(mp, uλqu2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, lower_bound = 0, upper_bound = λbar); # bilinear term
        @variable(mp, λql[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, upper_bound = 0);
        @variable(mp, uλql1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, upper_bound = 0, lower_bound = -λbar); # bilinear term
        @variable(mp, uλql2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, upper_bound = 0, lower_bound = -λbar); # bilinear term

        # add the V constraint
        @constraint(mp, V_constr[k in K], mp[:V] >= sum((sum(D[i,t] * π_add["P"][i,t] + probData.Pub[i] * π_add["Pu"][i,t] + probData.Plb[i] * π_add["Pl"][i,t] for i in probData.IDList) + 
                                                sum(π_add["d"][i,t] * mp[:x][i] for i in NCP)) + 
                                                sum((π_add["Wu"][br[1],br[2],t] - π_add["Wl"][br[1],br[2],t]) * probData.Wbar[br[1],br[2]] +
                                                    (π_add["θu"][br[1],br[2],t] - π_add["θl"][br[1],br[2],t]) * probData.θdiff[br[1],br[2]] for br in probData.brList) for t in T) + 
                                                mp[:z][k], container = SparseAxisArray);
        # add the z constraint
        # @constraint(mp, z_constr[k in K], mp[:z][k] >= sum(sum(mp[:λQ][k,j,t] * probData.Q[j,t] + sum(mp[:λu][k,i,j,t] * probData.Q[j,t] * mp[:u][i] for i in NCP) + 
        #                 sum(M * (2 - mp[:u][br[1]] - mp[:u][br[2]]) * (mp[:λqu][k,br[1],br[2],j,t] - mp[:λql][k,br[1],br[2],j,t]) for br in probData.brList) for j in J) for t in T), container = SparseAxisArray);
        @constraint(mp, z_constr[k in K], mp[:z][k] >= sum(sum(mp[:λQ][k,j,t] * probData.Q[j,t] + sum(mp[:uλu][k,i,j,t] * probData.Q[j,t] for i in NCP) + 
                        sum(M * (2 * (mp[:λqu][k,br[1],br[2],j,t] - mp[:λql][k,br[1],br[2],j,t]) - mp[:uλqu1][k,br[1],br[2],j,t] - mp[:uλqu2][k,br[1],br[2],j,t] + mp[:uλql1][k,br[1],br[2],j,t] + mp[:uλql2][k,br[1],br[2],j,t]) 
                        for br in probData.brList if (br[1] in NCP)&(br[2] in NCP)) for j in J) for t in T), container = SparseAxisArray);
        # add the λ constraints
        @constraint(mp, λ_constr[k in K, i in NCP, j in J, t in T], mp[:λQ][k,j,t] + mp[:λu][k,i,j,t] + sum(mp[:λqu][k,br[1],br[2],j,t] + mp[:λql][k,br[1],br[2],j,t] for br in probData.brList if (br[1] == i)&(br[1] in NCP)&(br[2] in NCP)) - 
                        sum((probData.r[br[2],j]/probData.r[br[1],j]) * (1 + Delta) * mp[:λqu][k,br[1],br[2],j,t] + 
                            (probData.r[br[2],j]/probData.r[br[1],j]) * (1 - Delta) * mp[:λql][k,br[1],br[2],j,t] for br in probData.brList if (br[1] in NCP)&(br[2] in NCP)&(br[2] == i)) >= π_add["q"][i,t] + probData.cq[i,j], container = SparseAxisArray);
        
        # add the λ-u logic constraints
        @constraint(mp, λqu_logic_constr1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:λqu][k,i,ip,j,t] <= M * mp[:u][i], container = SparseAxisArray);
        @constraint(mp, λqu_logic_constr2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:λqu][k,i,ip,j,t] <= M * mp[:u][ip], container = SparseAxisArray);
        @constraint(mp, λql_logic_constr1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:λql][k,i,ip,j,t] >= -M * mp[:u][i], container = SparseAxisArray);
        @constraint(mp, λql_logic_constr2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:λql][k,i,ip,j,t] >= -M * mp[:u][ip], container = SparseAxisArray);

        # add McCormick envelopes for bilinear terms
        @constraint(mp, uλu_constr1[k in K, i in NCP, j in J, t in T], mp[:uλu][k,i,j,t] <= λbar * mp[:u][i], container = SparseAxisArray);
        @constraint(mp, uλu_constr2[k in K, i in NCP, j in J, t in T], mp[:uλu][k,i,j,t] <= mp[:λu][k,i,j,t], container = SparseAxisArray);
        @constraint(mp, uλu_constr3[k in K, i in NCP, j in J, t in T], mp[:uλu][k,i,j,t] >= mp[:λu][k,i,j,t] - λbar * (1 - mp[:u][i]), container = SparseAxisArray);

        @constraint(mp, uλqu1_constr1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:uλqu1][k,i,ip,j,t] <= λbar * mp[:u][i], container = SparseAxisArray);
        @constraint(mp, uλqu1_constr2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:uλqu1][k,i,ip,j,t] <= mp[:λqu][k,i,ip,j,t], container = SparseAxisArray);
        @constraint(mp, uλqu1_constr3[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:uλqu1][k,i,ip,j,t] >= mp[:λqu][k,i,ip,j,t] - λbar * (1 - mp[:u][i]), container = SparseAxisArray);
        @constraint(mp, uλqu2_constr1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:uλqu2][k,i,ip,j,t] <= λbar * mp[:u][ip], container = SparseAxisArray);
        @constraint(mp, uλqu2_constr2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:uλqu2][k,i,ip,j,t] <= mp[:λqu][k,i,ip,j,t], container = SparseAxisArray);
        @constraint(mp, uλqu2_constr3[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:uλqu2][k,i,ip,j,t] >= mp[:λqu][k,i,ip,j,t] - λbar * (1 - mp[:u][ip]), container = SparseAxisArray);

        @constraint(mp, uλql1_constr1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:uλql1][k,i,ip,j,t] >= -λbar * mp[:u][i], container = SparseAxisArray);
        @constraint(mp, uλql1_constr2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:uλql1][k,i,ip,j,t] >= mp[:λql][k,i,ip,j,t], container = SparseAxisArray);
        @constraint(mp, uλql1_constr3[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:uλql1][k,i,ip,j,t] <= mp[:λql][k,i,ip,j,t] + λbar * (1 - mp[:u][i]), container = SparseAxisArray);
        @constraint(mp, uλql2_constr1[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:uλql2][k,i,ip,j,t] >= -λbar * mp[:u][ip], container = SparseAxisArray);
        @constraint(mp, uλql2_constr2[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:uλql2][k,i,ip,j,t] >= mp[:λql][k,i,ip,j,t], container = SparseAxisArray);
        @constraint(mp, uλql2_constr3[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], mp[:uλql2][k,i,ip,j,t] <= mp[:λql][k,i,ip,j,t] + λbar * (1 - mp[:u][ip]), container = SparseAxisArray);
    end
    return mp;
end

function addCol_bilinear(mp, probData, Delta, π_add, K, k, M = 1e4)
    NCP = probData.NCP;
    J = probData.J;
    T = 1:probData.T;

    if k > 1
        # add a new column to the master problem here
        mp[:z][k] = @variable(mp, base_name = "z[$k]");
        for j in J
            for t in T
                mp[:λQ][k,j,t] = @variable(mp, base_name = "λQ[$k,$j,$t]");
                for i in probData.NCP
                    mp[:λu][k,i,j,t] = @variable(mp, lower_bound = 0, base_name = "λu[$k,$i,$j,$t]");
                    for ip in probData.NCP
                        if ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)
                            mp[:λqu][k,i,ip,j,t] = @variable(mp, lower_bound = 0, base_name = "λqu[$k,$i,$ip,$j,$t]");
                            mp[:λql][k,i,ip,j,t] = @variable(mp, upper_bound = 0, base_name = "λql[$k,$i,$ip,$j,$t]");
                        end
                    end
                end
            end
        end
        # add the V constraint
        mp[:V_constr][k] = @constraint(mp, mp[:V] >= sum((sum(D[i,t] * π_add["P"][i,t] + probData.Pub[i] * π_add["Pu"][i,t] + probData.Plb[i] * π_add["Pl"][i,t] for i in probData.IDList) + 
                                                        sum(π_add["d"][i,t] * mp[:x][i] for i in NCP)) + 
                                                        sum((π_add["Wu"][br[1],br[2],t] - π_add["Wl"][br[1],br[2],t]) * probData.Wbar[br[1],br[2]] +
                                                            (π_add["θu"][br[1],br[2],t] - π_add["θl"][br[1],br[2],t]) * probData.θdiff[br[1],br[2]] for br in probData.brList) for t in T) + 
                                                        mp[:z][k]);
        # add the z constraint
        mp[:z_constr][k] = @constraint(mp, mp[:z][k] >= sum(sum(mp[:λQ][k,j,t] * probData.Q[j,t] + sum(mp[:λu][k,i,j,t] * probData.Q[j,t] * mp[:u][i] for i in NCP) + 
                                sum(M * (2 - mp[:u][br[1]] - mp[:u][br[2]]) * (mp[:λqu][k,br[1],br[2],j,t] - mp[:λql][k,br[1],br[2],j,t]) for br in probData.brList if (br[1] in NCP)&(br[2] in NCP)) for j in J) for t in T));
        # add the λ constraints
        for i in probData.NCP
            for j in J
                for t in T
                    mp[:λ_constr][k,i,j,t] = @constraint(mp, mp[:λQ][k,j,t] + mp[:λu][k,i,j,t] + sum(mp[:λqu][k,br[1],br[2],j,t] + mp[:λql][k,br[1],br[2],j,t] for br in probData.brList if (br[1] in NCP)&(br[2] in NCP)&(br[1] == i)) - 
                                sum((probData.r[br[2],j]/probData.r[br[1],j]) * (1 + Delta) * mp[:λqu][k,br[1],br[2],j,t] + 
                                    (probData.r[br[2],j]/probData.r[br[1],j]) * (1 - Delta) * mp[:λql][k,br[1],br[2],j,t] for br in probData.brList if (br[1] in NCP)&(br[2] in NCP)&(br[2] == i)) >= π_add["q"][i,t] + probData.cq[i,j]);
                end
            end
        end
    else
        @variable(mp, z[k in K], container = SparseAxisArray);
        @variable(mp, λQ[k in K, j in J, t in T], container = SparseAxisArray);
        @variable(mp, λu[k in K, i in NCP, j in J, t in T], container = SparseAxisArray, lower_bound = 0);
        @variable(mp, λqu[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, lower_bound = 0);
        @variable(mp, λql[k in K, i in NCP, ip in NCP, j in J, t in T; ((i, ip) in probData.brList)&(i in NCP)&(ip in NCP)], container = SparseAxisArray, upper_bound = 0);
        # add the V constraint
        @constraint(mp, V_constr[k in K], mp[:V] >= sum((sum(D[i,t] * π_add["P"][i,t] + probData.Pub[i] * π_add["Pu"][i,t] + probData.Plb[i] * π_add["Pl"][i,t] for i in probData.IDList) + 
                                                sum(π_add["d"][i,t] * mp[:x][i] for i in NCP)) + 
                                                sum((π_add["Wu"][br[1],br[2],t] - π_add["Wl"][br[1],br[2],t]) * probData.Wbar[br[1],br[2]] +
                                                    (π_add["θu"][br[1],br[2],t] - π_add["θl"][br[1],br[2],t]) * probData.θdiff[br[1],br[2]] for br in probData.brList) for t in T) + 
                                                mp[:z][k], container = SparseAxisArray);
        # add the z constraint
        @constraint(mp, z_constr[k in K], mp[:z][k] >= sum(sum(mp[:λQ][k,j,t] * probData.Q[j,t] + sum(mp[:λu][k,i,j,t] * probData.Q[j,t] * mp[:u][i] for i in NCP) + 
                        sum(M * (2 - mp[:u][br[1]] - mp[:u][br[2]]) * (mp[:λqu][k,br[1],br[2],j,t] - mp[:λql][k,br[1],br[2],j,t]) for br in probData.brList if (br[1] in NCP)&(br[2] in NCP)) for j in J) for t in T), container = SparseAxisArray);
        # add the λ constraints
        @constraint(mp, λ_constr[k in K, i in NCP, j in J, t in T], mp[:λQ][k,j,t] + mp[:λu][k,i,j,t] + sum(mp[:λqu][k,br[1],br[2],j,t] + mp[:λql][k,br[1],br[2],j,t] for br in probData.brList if (br[1] in NCP)&(br[2] in NCP)&(br[1] == i)) - 
                        sum((probData.r[br[2],j]/probData.r[br[1],j]) * (1 + Delta) * mp[:λqu][k,br[1],br[2],j,t] + 
                            (probData.r[br[2],j]/probData.r[br[1],j]) * (1 - Delta) * mp[:λql][k,br[1],br[2],j,t] for br in probData.brList if (br[1] in NCP)&(br[2] in NCP)&(br[2] == i)) >= π_add["q"][i,t] + probData.cq[i,j], container = SparseAxisArray);
    end
    return mp;
end

function build_separation_local(probData, ρ, Delta, xhat, uhat, M = 1e4)
    # build the separation problem here
    NCP = probData.NCP;
    J = probData.J;
    T = 1:probData.T;

    sp = Model(Ipopt.Optimizer);
    set_attribute(sp, "print_level", 0);
    set_attribute(sp, "hsllib", HSL_jll.libhsl_path);
    set_attribute(sp, "linear_solver", "ma27");

    @variable(sp, πP[i in probData.IDList, t in T]); # dual variable for active power balance
    @variable(sp, πPu[i in probData.IDList, t in T] <= 0); # dual variable for active power upper bound
    @variable(sp, πPl[i in probData.IDList, t in T] >= 0); # dual variable for active power lower bound
    @variable(sp, πq[i in NCP, t in T]); # dual variable for EV demand balance
    @variable(sp, πd[i in NCP, t in T] <= 0); # dual variable for EV demand upper bound
    @variable(sp, πW[br in probData.brList, t in T]); # dual variable for branch flow equation
    @variable(sp, πWu[br in probData.brList, t in T] <= 0); # dual variable for branch flow upper bound
    @variable(sp, πWl[br in probData.brList, t in T] >= 0); # dual variable for branch flow lower bound
    @variable(sp, πθu[br in probData.brList, t in T] <= 0); # dual variable for angle difference upper bound
    @variable(sp, πθl[br in probData.brList, t in T] >= 0); # dual variable for angle difference lower bound
    @variable(sp, q[i in NCP, j in J, t in T] >= 0); # EV demand

    # objective function
    @objective(sp, Max, sum((sum(D[i,t] * πP[i,t] + probData.Pub[i] * πPu[i,t] + probData.Plb[i] * πPl[i,t] for i in probData.IDList) + 
                                sum(πd[i,t] * xhat[i] for i in NCP)) + 
                                sum((πWu[br,t] - πWl[br,t]) * probData.Wbar[br[1],br[2]] + (πθu[br,t] - πθl[br,t]) * probData.θdiff[br[1],br[2]] for br in probData.brList) +
                                sum(sum(q[i,j,t] for j in J) * πq[i,t] + sum(probData.cq[i,j] * q[i,j,t] for j in J) for i in NCP) for t in T));

    # add the constraints
    @constraint(sp, d_constr[i in NCP, t in T], πq[i,t] + πd[i,t] - πP[i,t] <= 0);
    @constraint(sp, s_constr[i in NCP, t in T], πq[i,t] <= ρ);
    @constraint(sp, theta_constr[i in probData.IDList, t in T], sum(πθu[br,t] + πθl[br,t] - probData.b[br] * πW[br,t] for br in probData.brList if br[1] == i) - 
                        sum(πθu[br,t] + πθl[br,t] - probData.b[br] * πW[br,t] for br in probData.brList if br[2] == i) == 0);
    @constraint(sp, W_constr[br in probData.brList, t in T], πW[br,t] + πWu[br,t] + πWl[br,t] + πP[br[2],t] - πP[br[1],t] == 0);
    @constraint(sp, P_constr[i in probData.IDList, t in T], πP[i,t] + πPu[i,t] + πPl[i,t] == probData.g[i]);
    @constraint(sp, Q_constr[j in J, t in T], sum(q[i,j,t] for i in NCP) == probData.Q[j,t]);
    @constraint(sp, q_constr[i in NCP, j in J, t in T], q[i,j,t] <= probData.Q[j,t] * uhat[i]);
    @constraint(sp, qu_constr[br in probData.brList, j in J, t in T; (br[1] in NCP)&(br[2] in NCP)], q[br[1],j,t] - (probData.r[br[2],j]/probData.r[br[1],j]) * (1 + Delta) * q[br[2],j,t] <= 
                        M * (2 - uhat[br[1]] - uhat[br[2]]));
    @constraint(sp, ql_constr[br in probData.brList, j in J, t in T; (br[1] in NCP)&(br[2] in NCP)], q[br[1],j,t] - (probData.r[br[2],j]/probData.r[br[1],j]) * (1 - Delta) * q[br[2],j,t] >= 
                        -M * (2 - uhat[br[1]] - uhat[br[2]]));
    
    return sp;
end

function build_separation_bilinear(probData, ρ, Delta, xhat, uhat, M = 1e4)
    # build the separation problem here
    NCP = probData.NCP;
    J = probData.J;
    T = 1:probData.T;

    sp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 1, "TimeLimit" => 600, "Threads" => 10));
    @variable(sp, πP[i in probData.IDList, t in T]); # dual variable for active power balance
    @variable(sp, πPu[i in probData.IDList, t in T] <= 0); # dual variable for active power upper bound
    @variable(sp, πPl[i in probData.IDList, t in T] >= 0); # dual variable for active power lower bound
    @variable(sp, πq[i in NCP, t in T]); # dual variable for EV demand balance
    @variable(sp, πd[i in NCP, t in T] <= 0); # dual variable for EV demand upper bound
    @variable(sp, πW[br in probData.brList, t in T]); # dual variable for branch flow equation
    @variable(sp, πWu[br in probData.brList, t in T] <= 0); # dual variable for branch flow upper bound
    @variable(sp, πWl[br in probData.brList, t in T] >= 0); # dual variable for branch flow lower bound
    @variable(sp, πθu[br in probData.brList, t in T] <= 0); # dual variable for angle difference upper bound
    @variable(sp, πθl[br in probData.brList, t in T] >= 0); # dual variable for angle difference lower bound
    @variable(sp, q[i in NCP, j in J, t in T] >= 0); # EV demand

    # objective function
    @objective(sp, Max, sum((sum(D[i,t] * πP[i,t] + probData.Pub[i] * πPu[i,t] + probData.Plb[i] * πPl[i,t] for i in probData.IDList) + 
                                sum(πd[i,t] * xhat[i] for i in NCP)) + 
                                sum((πWu[br,t] - πWl[br,t]) * probData.Wbar[br[1],br[2]] + (πθu[br,t] - πθl[br,t]) * probData.θdiff[br[1],br[2]] for br in probData.brList) +
                                sum(sum(q[i,j,t] for j in J) * πq[i,t] + sum(probData.cq[i,j] * q[i,j,t] for j in J) for i in NCP) for t in T));

    # add the constraints
    @constraint(sp, d_constr[i in NCP, t in T], πq[i,t] + πd[i,t] - πP[i,t] <= 0);
    @constraint(sp, s_constr[i in NCP, t in T], πq[i,t] <= ρ);
    @constraint(sp, theta_constr[i in probData.IDList, t in T], sum(πθu[br,t] + πθl[br,t] - probData.b[br] * πW[br,t] for br in probData.brList if br[1] == i) - 
                        sum(πθu[br,t] + πθl[br,t] - probData.b[br] * πW[br,t] for br in probData.brList if br[2] == i) == 0);
    @constraint(sp, W_constr[br in probData.brList, t in T], πW[br,t] + πWu[br,t] + πWl[br,t] + πP[br[2],t] - πP[br[1],t] == 0);
    @constraint(sp, P_constr[i in probData.IDList, t in T], πP[i,t] + πPu[i,t] + πPl[i,t] == probData.g[i]);
    @constraint(sp, Q_constr[j in J, t in T], sum(q[i,j,t] for i in NCP) == probData.Q[j,t]);
    @constraint(sp, q_constr[i in NCP, j in J, t in T], q[i,j,t] <= probData.Q[j,t] * uhat[i]);
    @constraint(sp, qu_constr[br in probData.brList, j in J, t in T; (br[1] in NCP)&(br[2] in NCP)], q[br[1],j,t] - (probData.r[br[2],j]/probData.r[br[1],j]) * (1 + Delta) * q[br[2],j,t] <= 
                        M * (2 - uhat[br[1]] - uhat[br[2]]));
    @constraint(sp, ql_constr[br in probData.brList, j in J, t in T; (br[1] in NCP)&(br[2] in NCP)], q[br[1],j,t] - (probData.r[br[2],j]/probData.r[br[1],j]) * (1 - Delta) * q[br[2],j,t] >= 
                        -M * (2 - uhat[br[1]] - uhat[br[2]]));
    
    return sp;
end

function build_separation_piecewise_cuts(probData, ρ, Delta, xhat, uhat, M = 1e4)
    # build the separation problem here
    NCP = probData.NCP;
    J = probData.J;
    T = 1:probData.T;
    qList = [];
    iter_bool = true;
    sp_lblist = [];
    πList = [];

    sp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0, "Threads" => 10));
    @variable(sp, πP[i in probData.IDList, t in T]); # dual variable for active power balance
    @variable(sp, πPu[i in probData.IDList, t in T] <= 0); # dual variable for active power upper bound
    @variable(sp, πPl[i in probData.IDList, t in T] >= 0); # dual variable for active power lower bound
    @variable(sp, πq[i in NCP, t in T]); # dual variable for EV demand balance
    @variable(sp, πd[i in NCP, t in T] <= 0); # dual variable for EV demand upper bound
    @variable(sp, πW[br in probData.brList, t in T]); # dual variable for branch flow equation
    @variable(sp, πWu[br in probData.brList, t in T] <= 0); # dual variable for branch flow upper bound
    @variable(sp, πWl[br in probData.brList, t in T] >= 0); # dual variable for branch flow lower bound
    @variable(sp, πθu[br in probData.brList, t in T] <= 0); # dual variable for angle difference upper bound
    @variable(sp, πθl[br in probData.brList, t in T] >= 0); # dual variable for angle difference lower bound
    @variable(sp, qπ[i in NCP, t in T] >= 0); # bilinear term for π * q

    # objective function
    @objective(sp, Max, sum((sum(D[i,t] * πP[i,t] + probData.Pub[i] * πPu[i,t] + probData.Plb[i] * πPl[i,t] for i in probData.IDList) + 
                                sum(πd[i,t] * xhat[i] for i in NCP)) + 
                                sum((πWu[br,t] - πWl[br,t]) * probData.Wbar[br[1],br[2]] + (πθu[br,t] - πθl[br,t]) * probData.θdiff[br[1],br[2]] for br in probData.brList) for t in T));

    # add the constraints
    @constraint(sp, d_constr[i in NCP, t in T], πq[i,t] + πd[i,t] - πP[i,t] <= 0);
    @constraint(sp, s_constr[i in NCP, t in T], πq[i,t] <= ρ);
    @constraint(sp, theta_constr[i in probData.IDList, t in T], sum(πθu[br,t] + πθl[br,t] - probData.b[br] * πW[br,t] for br in probData.brList if br[1] == i) - 
                        sum(πθu[br,t] + πθl[br,t] - probData.b[br] * πW[br,t] for br in probData.brList if br[2] == i) == 0);
    @constraint(sp, W_constr[br in probData.brList, t in T], πW[br,t] + πWu[br,t] + πWl[br,t] + πP[br[2],t] - πP[br[1],t] == 0);
    @constraint(sp, P_constr[i in probData.IDList, t in T], πP[i,t] + πPu[i,t] + πPl[i,t] == probData.g[i]);

    while iter_bool
        # solve the π problem and obtain a πhat
        optimize!(sp);
        push!(sp_lblist, objective_value(sp));
        πqhat = Dict();
        for i in NCP
            for t in T
                πqhat[i,t] = value(sp[:πq][i,t]);
            end
        end

        if πqhat in πList
            iter_bool = false;
        else
            push!(πList, πqhat);
            # solve q problem to generate a slope for πq
            qp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0));
            @variable(qp, q[i in NCP, j in J, t in T] >= 0); # EV demand
            @constraint(qp, Q_constr[j in J, t in T], sum(q[i,j,t] for i in NCP) == probData.Q[j,t]);
            @constraint(qp, q_constr[i in NCP, j in J, t in T], q[i,j,t] <= probData.Q[j,t] * uhat[i]);
            @constraint(qp, qu_constr[br in probData.brList, j in J, t in T; (br[1] in NCP)&(br[2] in NCP)], q[br[1],j,t] - (probData.r[br[2],j]/probData.r[br[1],j]) * (1 + Delta) * q[br[2],j,t] <= 
                                M * (2 - uhat[br[1]] - uhat[br[2]]));
            @constraint(qp, ql_constr[br in probData.brList, j in J, t in T; (br[1] in NCP)&(br[2] in NCP)], q[br[1],j,t] - (probData.r[br[2],j]/probData.r[br[1],j]) * (1 - Delta) * q[br[2],j,t] >= 
                                -M * (2 - uhat[br[1]] - uhat[br[2]]));
            @objective(qp, Max, sum(sum(sum(q[i,j,t] for j in J) * πqhat[i,t] for i in NCP) for t in T) + sum(sum(sum(probData.cq[i,j] * q[i,j,t] for j in J) for i in NCP) for t in T));
            optimize!(qp);
            qsol = Dict();
            for i in NCP
                for j in J
                    for t in T
                        qsol[i,j,t] = value(qp[:q][i,j,t]);
                    end
                end
            end
            push!(qList, qsol);
            L = length(qList);
            if L == 1
                @variable(sp, po[l in 1:L], Bin, container = SparseAxisArray);
                @variable(sp, 0 <= popi[i in NCP, t in T, l in 1:L] <= ρ, container = SparseAxisArray);
                @constraint(sp, sumOne, sum(po[l] for l in 1:L) == 1);
                @constraint(sp, popi_constr1[i in NCP, t in T, l in 1:L], popi[i,t,l] <= ρ * po[l], container = SparseAxisArray);
                @constraint(sp, popi_constr2[i in NCP, t in T, l in 1:L], popi[i,t,l] <= πq[i,t], container = SparseAxisArray);
                @constraint(sp, popi_constr3[i in NCP, t in T, l in 1:L], popi[i,t,l] >= πq[i,t] + ρ * (po[l] - 1), container = SparseAxisArray);
                @objective(sp, Max, sum((sum(D[i,t] * πP[i,t] + probData.Pub[i] * πPu[i,t] + probData.Plb[i] * πPl[i,t] for i in probData.IDList) + 
                            sum(πd[i,t] * xhat[i] for i in NCP)) + 
                            sum((πWu[br,t] - πWl[br,t]) * probData.Wbar[br[1],br[2]] + (πθu[br,t] - πθl[br,t]) * probData.θdiff[br[1],br[2]] for br in probData.brList) +
                            sum(sum(popi[i,t,l] * sum(qList[l][i,j,t] for j in J) + sum(qList[l][i,j,t] * po[l] * probData.cq[i,j] for j in J) for l in 1:L) for i in NCP) for t in T));
            else
                sp[:po][L] = @variable(sp, binary = true, base_name = "po[$L]");
                for i in NCP
                    for t in T
                        sp[:popi][i,t,L] = @variable(sp, lower_bound = 0, upper_bound = ρ, base_name = "popi[$i,$t,$L]");
                    end
                end
                delete(sp, sp[:sumOne]);
                unregister(sp, :sumOne);
                @constraint(sp, sumOne, sum(sp[:po][l] for l in 1:L) == 1);
                for i in NCP
                    for t in T
                        sp[:popi_constr1][i,t,L] = @constraint(sp, sp[:popi][i,t,L] <= ρ * sp[:po][L], base_name = "popi_constr1[$i,$t,$L]");
                        sp[:popi_constr2][i,t,L] = @constraint(sp, sp[:popi][i,t,L] <= sp[:πq][i,t], base_name = "popi_constr2[$i,$t,$L]");
                        sp[:popi_constr3][i,t,L] = @constraint(sp, sp[:popi][i,t,L] >= sp[:πq][i,t] + ρ * (sp[:po][L] - 1), base_name = "popi_constr3[$i,$t,$L]");
                    end
                end
                @objective(sp, Max, sum((sum(D[i,t] * πP[i,t] + probData.Pub[i] * πPu[i,t] + probData.Plb[i] * πPl[i,t] for i in probData.IDList) + 
                    sum(πd[i,t] * xhat[i] for i in NCP)) + 
                    sum((πWu[br,t] - πWl[br,t]) * probData.Wbar[br[1],br[2]] + (πθu[br,t] - πθl[br,t]) * probData.θdiff[br[1],br[2]] for br in probData.brList) +
                    sum(sum(sp[:popi][i,t,l] * sum(qList[l][i,j,t] for j in J) + sum(qList[l][i,j,t] * sp[:po][l] * probData.cq[i,j] for j in J) for l in 1:L) for i in NCP) for t in T));
            end
        end
    end

    return sp;
end

function build_separation_bilinear_decomp(probData, ρ, Delta, xhat, uhat, M = 1e4, outputOption = 0)
    # solve smaller decomposed bilinear optimization model for the separation problem
    # build the separation problem here
    NCP = probData.NCP;
    J = probData.J;
    T = 1:probData.T;
    qhat = Dict();

    for t in T
        spt = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0, "TimeLimit" => 600, "Threads" => 10));
        @variable(spt, πP[i in probData.IDList]); # dual variable for active power balance
        @variable(spt, πPu[i in probData.IDList] <= 0); # dual variable for active power upper bound
        @variable(spt, πPl[i in probData.IDList] >= 0); # dual variable for active power lower bound
        @variable(spt, πq[i in NCP]); # dual variable for EV demand balance
        @variable(spt, πd[i in NCP] <= 0); # dual variable for EV demand upper bound
        @variable(spt, πW[br in probData.brList]); # dual variable for branch flow equation
        @variable(spt, πWu[br in probData.brList] <= 0); # dual variable for branch flow upper bound
        @variable(spt, πWl[br in probData.brList] >= 0); # dual variable for branch flow lower bound
        @variable(spt, πθu[br in probData.brList] <= 0); # dual variable for angle difference upper bound
        @variable(spt, πθl[br in probData.brList] >= 0); # dual variable for angle difference lower bound
        @variable(spt, q[i in NCP, j in J] >= 0); # EV demand

        # objective function
        @objective(spt, Max, (sum(D[i,t] * πP[i] + probData.Pub[i] * πPu[i] + probData.Plb[i] * πPl[i] for i in probData.IDList) + 
                                sum(πd[i] * xhat[i] for i in NCP)) + 
                                sum((πWu[br] - πWl[br]) * probData.Wbar[br[1],br[2]] + (πθu[br] - πθl[br]) * probData.θdiff[br[1],br[2]] for br in probData.brList) +
                                sum(sum(q[i,j] for j in J) * πq[i] + sum(probData.cq[i,j] * q[i,j] for j in J) for i in NCP));

        # add the constraints
        @constraint(spt, d_constr[i in NCP], πq[i] + πd[i] - πP[i] <= 0);
        @constraint(spt, s_constr[i in NCP], πq[i] <= ρ);
        @constraint(spt, theta_constr[i in probData.IDList], sum(πθu[br] + πθl[br] - probData.b[br] * πW[br] for br in probData.brList if br[1] == i) - 
                            sum(πθu[br] + πθl[br] - probData.b[br] * πW[br] for br in probData.brList if br[2] == i) == 0);
        @constraint(spt, W_constr[br in probData.brList], πW[br] + πWu[br] + πWl[br] + πP[br[2]] - πP[br[1]] == 0);
        @constraint(spt, P_constr[i in probData.IDList], πP[i] + πPu[i] + πPl[i] == probData.g[i]);
        @constraint(spt, Q_constr[j in J], sum(q[i,j] for i in NCP) == probData.Q[j,t]);
        @constraint(spt, q_constr[i in NCP, j in J], q[i,j] <= probData.Q[j,t] * uhat[i]);
        @constraint(spt, qu_constr[br in probData.brList, j in J; (br[1] in NCP)&(br[2] in NCP)], q[br[1],j] - (probData.r[br[2],j]/probData.r[br[1],j]) * (1 + Delta) * q[br[2],j] <= 
                            M * (2 - uhat[br[1]] - uhat[br[2]]));
        @constraint(spt, ql_constr[br in probData.brList, j in J; (br[1] in NCP)&(br[2] in NCP)], q[br[1],j] - (probData.r[br[2],j]/probData.r[br[1],j]) * (1 - Delta) * q[br[2],j] >= 
                            -M * (2 - uhat[br[1]] - uhat[br[2]]));
        optimize!(spt);

        # obtain qhat
        for i in NCP
            for j in J
                qhat[i,j,t] = value(spt[:q][i,j]);
            end
        end
    end

    sp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0, "TimeLimit" => 600, "Threads" => 10));
    @variable(sp, πP[i in probData.IDList, t in T]); # dual variable for active power balance
    @variable(sp, πPu[i in probData.IDList, t in T] <= 0); # dual variable for active power upper bound
    @variable(sp, πPl[i in probData.IDList, t in T] >= 0); # dual variable for active power lower bound
    @variable(sp, πq[i in NCP, t in T]); # dual variable for EV demand balance
    @variable(sp, πd[i in NCP, t in T] <= 0); # dual variable for EV demand upper bound
    @variable(sp, πW[br in probData.brList, t in T]); # dual variable for branch flow equation
    @variable(sp, πWu[br in probData.brList, t in T] <= 0); # dual variable for branch flow upper bound
    @variable(sp, πWl[br in probData.brList, t in T] >= 0); # dual variable for branch flow lower bound
    @variable(sp, πθu[br in probData.brList, t in T] <= 0); # dual variable for angle difference upper bound
    @variable(sp, πθl[br in probData.brList, t in T] >= 0); # dual variable for angle difference lower bound

    # objective function
    @objective(sp, Max, sum((sum(D[i,t] * πP[i,t] + probData.Pub[i] * πPu[i,t] + probData.Plb[i] * πPl[i,t] for i in probData.IDList) + 
                                sum(πd[i,t] * xhat[i] for i in NCP)) + 
                                sum((πWu[br,t] - πWl[br,t]) * probData.Wbar[br[1],br[2]] + (πθu[br,t] - πθl[br,t]) * probData.θdiff[br[1],br[2]] for br in probData.brList) +
                                sum(sum(qhat[i,j,t] for j in J) * πq[i,t] + sum(probData.cq[i,j] * qhat[i,j,t] for j in J) for i in NCP) for t in T));

    # add the constraints
    @constraint(sp, d_constr[i in NCP, t in T], πq[i,t] + πd[i,t] - πP[i,t] <= 0);
    @constraint(sp, s_constr[i in NCP, t in T], πq[i,t] <= ρ);
    @constraint(sp, theta_constr[i in probData.IDList, t in T], sum(πθu[br,t] + πθl[br,t] - probData.b[br] * πW[br,t] for br in probData.brList if br[1] == i) - 
                        sum(πθu[br,t] + πθl[br,t] - probData.b[br] * πW[br,t] for br in probData.brList if br[2] == i) == 0);
    @constraint(sp, W_constr[br in probData.brList, t in T], πW[br,t] + πWu[br,t] + πWl[br,t] + πP[br[2],t] - πP[br[1],t] == 0);
    @constraint(sp, P_constr[i in probData.IDList, t in T], πP[i,t] + πPu[i,t] + πPl[i,t] == probData.g[i]);
    
    if outputOption == 0
        return sp;
    else
        return sp, qhat;
    end
end

function colGen(probData, f, c, xbar, x0, u0, ρ, Delta, ϵ = 1e-4, master_option = "linear", separation_option = "piecewise", time_limit = 3600)
    UB = Inf; 
    LB = 0;
    K = [];
    πList = [];
    iter = 0;
    xhat = Dict();
    uhat = Dict();
    xbest = Dict();
    ubest = Dict();
    UBList = [];
    LBList = [];
    timeList = [0.0];
    NCP = probData.NCP;
    J = probData.J;
    T = 1:probData.T;
    if master_option == "linear"
        mp = build_masterproblem(probData, f, c, xbar, x0, u0, Delta, K, πList);
    else
        mp = build_masterproblem_bilinear(probData, f, c, xbar, x0, u0, Delta, K, πList);
    end
    while ((LB == 0) | ((UB - LB)/LB > ϵ)) & (sum(timeList) < time_limit)
        # solve the relaxed master problem and obtain the incumbent first-stage solutions
        start_iter = time();
        optimize!(mp);
        #optimize!(mp_bilinear);
        LB = objective_value(mp);
        push!(LBList, LB);
        Vhat = value(mp[:V]);
        for i in NCP
            xhat[i] = value(mp[:x][i]);
            uhat[i] = value(mp[:u][i]);
        end

        # build and solve the separation problem
        if separation_option == "bilinear_gurobi"
            sp = build_separation_bilinear(probData, ρ, Delta, xhat, uhat);
        elseif separation_option == "bilinear_local"
            sp = build_separation_local(probData, ρ, Delta, xhat, uhat);
        elseif separation_option == "bilinear_decomp"
            sp = build_separation_bilinear_decomp(probData, ρ, Delta, xhat, uhat);
        else
            sp = build_separation_piecewise_cuts(probData, ρ, Delta, xhat, uhat);
        end
        optimize!(sp);
        Vtilde = objective_value(sp);
        if separation_option == "bilinear_local"
            println("Separation problem objective value: $(Vtilde)");
        end
        πhat = Dict();
        πhat["P"] = Dict();
        πhat["Pu"] = Dict();
        πhat["Pl"] = Dict();
        πhat["q"] = Dict();
        πhat["d"] = Dict();
        πhat["W"] = Dict();
        πhat["Wu"] = Dict();
        πhat["Wl"] = Dict();
        πhat["θu"] = Dict();
        πhat["θl"] = Dict();
        # obtain the πhat values
        for t in T
            for i in probData.IDList
                πhat["P"][i,t] = value(sp[:πP][i,t]);
                πhat["Pu"][i,t] = value(sp[:πPu][i,t]);
                πhat["Pl"][i,t] = value(sp[:πPl][i,t]);
            end
            for i in probData.NCP
                πhat["q"][i,t] = value(sp[:πq][i,t]);
                πhat["d"][i,t] = value(sp[:πd][i,t]);
            end
            for br in probData.brList
                πhat["W"][br[1],br[2],t] = value(sp[:πW][br,t]);
                πhat["Wu"][br[1],br[2],t] = value(sp[:πWu][br,t]);
                πhat["Wl"][br[1],br[2],t] = value(sp[:πWl][br,t]);
                πhat["θu"][br[1],br[2],t] = value(sp[:πθu][br,t]);
                πhat["θl"][br[1],br[2],t] = value(sp[:πθl][br,t]);
            end
        end
        # obtain the πhat values
        push!(πList, πhat);

        # if the worst case V is larger, append the column
        if Vtilde > Vhat
            iter += 1;
            push!(K, iter);
            if master_option == "linear"
                mp = addCol(mp, probData, Delta, πhat, K, iter);
            elseif master_option == "bilinear"
                mp = addCol_bilinear(mp, probData, Delta, πhat, K, iter);
            end
        end

        # update the upper bound
        if UB > Vtilde + sum(c[i] * (xhat[i] - x0[i]) + f[i] * (uhat[i] - u0[i]) for i in NCP);
            UB = Vtilde + sum(c[i] * (xhat[i] - x0[i]) + f[i] * (uhat[i] - u0[i]) for i in NCP);
            for i in NCP
                xbest[i] = xhat[i];
                ubest[i] = uhat[i];
            end
        end
        push!(UBList, UB);
        iter_elapse = time() - start_iter;
        push!(timeList, iter_elapse);

        println("Current LB: $LB, Current UB: $UB, Gap: $((UB - LB)/LB), Iteration Time: $(iter_elapse)");
    end
    return xbest, ubest, UBList, LBList, timeList;
end

function solve_h(probData, xhat, qhat)
    # solve for the h function value in the primal problem form
    NCP = probData.NCP;
    J = probData.J;
    T = 1:probData.T;

    hsp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0));
    # set up the variables
    @variable(hsp, probData.Plb[i] <= P[i in probData.IDList, t in T] <= probData.Pub[i]); # active power generation
    @variable(hsp, -probData.Wbar[br] <= W[br in probData.brList, t in T] <= probData.Wbar[br]); # branch flow
    @variable(hsp, θ[i in probData.IDList, t in T]); # voltage angle
    @variable(hsp, d[i in NCP, t in T] >= 0); # EV satisfied demand
    @variable(hsp, s[i in NCP, t in T] >= 0); # EV unsatisfied demand

    # set up the constraints
    @constraint(hsp, EV_demand[i in NCP, t in T], d[i,t] + s[i,t] == sum(qhat[i,j,t] for j in J));
    @constraint(hsp, EV_demand_limit[i in NCP, t in T], d[i,t] <= xhat[i]);
    @constraint(hsp, DC_flow[br in probData.brList, t in T], W[br,t] == probData.b[br] * (θ[br[1],t] - θ[br[2],t]));
    @constraint(hsp, flow_balance_NEV[i in probData.IDList, t in T; !(i in NCP)], P[i,t] + sum(W[br,t] for br in probData.brList if br[2] == i) - sum(W[br,t] for br in probData.brList if br[1] == i) == probData.D[i,t]);
    @constraint(hsp, flow_balance_EV[i in NCP, t in T], P[i,t] + sum(W[br,t] for br in probData.brList if br[2] == i) - sum(W[br,t] for br in probData.brList if br[1] == i) == probData.D[i,t] + d[i,t]);
    @constraint(hsp, theta_limit[br in probData.brList, t in T], -probData.θdiff[br[1],br[2]] <= θ[br[1],t] - θ[br[2],t] <= probData.θdiff[br[1],br[2]]);
    @constraint(hsp, ref_angle[t in T], θ[7, t] == 0);

    # set up the objective function
    @objective(hsp, Min, sum(sum(probData.g[i] * P[i,t] for i in probData.IDList) + sum(ρ * s[i,t] for i in NCP) for t in T));
    @expression(hsp, generation_cost, sum(sum(probData.g[i] * P[i,t] for i in probData.IDList) for t in T));
    @expression(hsp, unsat_cost, sum(sum(ρ * s[i,t] for i in NCP) for t in T));
    return hsp;
end