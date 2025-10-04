# generate the random demand portfolio
using Distributions, Statistics, Random;

commercial_rate = [19, 7.1, 7.1, 7.1, 7.1, 7.1, 7.1, 5, 5, 9, 13.5, 23.4, 21.1, 23.4, 23.4, 25.7, 31.8, 33.2, 33.2, 33.2, 35.2, 37.5, 29.5, 23.6];
residential_rate = [56.5, 39.2, 31, 24, 18.6, 14.3, 11.1, 8.6, 7, 6, 7, 8.8, 11.1, 16, 22, 28.4, 35.2, 45.1, 55.7, 65, 68.6, 72, 70, 64];

Q0 = [0.8, 1.1, 1.5, 1.2, 0.6];
Qsigma = Q0 ./ 10;
Qj = zeros(24, length(Q0));
for i in eachindex(Q0)
    for j in 1:24
        if i in [1, 5]
            Qj[j,i] = max(0.0, rand(Normal(Q0[i] * commercial_rate[j] / maximum(commercial_rate), Qsigma[i])));
        else
            Qj[j,i] = max(0.0, rand(Normal(Q0[i] * residential_rate[j] / maximum(residential_rate), Qsigma[i])));
        end
    end
end
Qj = round.(Qj,digits = 5)
df = DataFrame(Qj, :auto);
CSV.write("./data/13/case13_Q.csv",df);