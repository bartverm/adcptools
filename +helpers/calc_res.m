function [mse_res, mean_res, std_res, res] = calc_res(rhs, lhs)

res = rhs - lhs;
mse_res = mean(res.^2);
mean_res = mean(res);
std_res = std(res);
end