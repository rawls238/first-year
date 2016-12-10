N = 1001
rows, cols = N, N
sps_acc = spzeros(N, N)
for i=1:N
    if i <= 11
        if i == 1
            sps_acc[1, 1] = 0.67
        else
            sps_acc[i, i] = 0.33
            sps_acc[i, 1] = 0.34
        end
        sps_acc[i, i+10] = 0.33
    elseif i > 11 && i < 991
        sps_acc[i, 1] = 0.01
        sps_acc[i, i - 10] = 0.33
        sps_acc[i, i] = 0.33
        sps_acc[i, i + 10] = 0.33
    elseif i >= 991 
        sps_acc[i, 1] = 0.01
        sps_acc[i, i - 10] = 0.33
        if i == N
            sps_acc[N, N] = 0.66
        else
            sps_acc[i, i] = 0.33
            sps_acc[i, N] = 0.33
        end
    end
end

dist = sps_acc^10000000000000

sum_of_vec = 0.0
var_of_vec = 0.0
skewness_of_vec_num = 0.0
skewness_of_vec_denom = 0.0
for i=1:N
    sum_of_vec += dist[1, i]
end

mean = sum_of_vec / N

for i=1:N
    var_of_vec += (dist[1, i] - mean)^2
end
var_of_vec = var_of_vec / (N-1)

for i=1:N
    skewness_of_vec_num += (dist[1, i] - mean)^3
    skewness_of_vec_denom += (dist[1, i] - mean)^2
end
skewness_of_vec_num /= N
skewness_of_vec_denom = (skewness_of_vec_denom / (N-1))^(0.5)
skewness = skewness_of_vec_num / skewness_of_vec_denom