clear;
clc;

// measurement uncertainties

// u_i -> uExt_i := k_i * u_i, k_i >= 1:
// 1. metrological compatibility: abs(x_i - x_j) / 2 * sqrt(k_i^2 * u_i^2 + k_j^2 * u_j^2) <= 1, 1 <= i, j <= n
// 2. sum_1^n (k_i * u_i)^2 -> min_{k_1, ..., k_n}
// N. A. Burmistrova, Development and study of algorithms for processing inconsistent data in key comparisons of standards, Measuremet Techniques 57 (10) (2015) 1103-1112

// === example ===

// measurement data
x = [1.498, 1.525, 1.554, 1.493, 1.480, 1.500, 1.529, 1.481, 1.535, 1.606];
// measurement uncertainties
u = [0.011, 0.006, 0.012, 0.032, 0.007, 0.011, 0.013, 0.008, 0.008, 0.007];

u2 = u .^ 2;
n = length(x);

function [res] = check(meas, unc)
    tmp = zeros(n, n);
    for i = 1 : n
        for j = (i + 1) : n
            tmp(i, j) = 0.5 * abs(meas(i) - meas(j)) / sqrt(unc(i)^2 + unc(j)^2);
            //tmp(j, i) = tmp(i, j);
        end
    end
    res =  max(tmp);
endfunction


// solving linear optimization problem:

A_cond = [];
b_cond = [];

for i = 1 : n
    for j = (i + 1) : n
        tmp = zeros(1, n);
        tmp(i) = u2(i);
        tmp(j) = u2(j);
        A_cond = [A_cond; tmp];
        b_cond = [b_cond; 0.25 * (x(i) - x(j))^2];
    end
end

ID = diag(ones(1, n));

A_cond = [A_cond; ID];
b_cond = [b_cond; (ones(1, n))'];

nY = length(b_cond);
Y = - diag(ones(1, nY));

A_cond = [A_cond, Y];

c = [u2, zeros(1, nY)]';

x0 = [];
a_opt = karmarkar(A_cond, b_cond, c, x0, 1.e-7);

a_opt = a_opt(1 : n);

// target function
f = u2 * a_opt

// expansion factors
k = sqrt(a_opt);
// rounding...
k = 1.e-4 * round(1.e4 * k)

// new uncertainties
uExt = k .* u'

// > 1
checkBefore = check(x, u)
// 1
checkAfter  = check(x, uExt)
