clear;
clc;

// u_i -> uExt_i := k_i * u_i, k_i >= 1:
// 1. metrological compatibility: abs(x_i - x_j) / 2 * sqrt(k_i^2 * u_i^2 + k_j^2 * u_j^2) <= 1, 1 <= i, j <= n
// 2. sum_1^n (k_i * u_i)^2 -> min_{k_1, ..., k_n}

// measurement data and uncertainties
x = [0.00248, 0.00291, 0.00399, 0.00427, 0.00559, 0.00422, 0.00385, 0.00221, 0.00425, 0.00349, 0.00260, 0.00553, 0.00191, 0.00435, 0.00419, 0.00449];
u = [0.00043, 0.00051, 0.00069, 0.00009, 0.00027, 0.00098, 0.00027, 0.00087, 0.00013, 0.00019, 0.00025, 0.00015, 0.00099, 0.00013, 0.00008, 0.00008];


CHI2CRIT = [3.841, 5.991, 7.815, 9.448, 11.070, 12.592, 14.067, 15.507, 16.919, 18.307, 19.675, 21.026, 22.362, 23.685, 24.996, 26.296, 27.587, 28.869, 30.144, 31.410, 32.671, 33.924, 35.172, 36.415, 37.652, 38.885, 40.113, 41.337, 42.557, 43.773];



// weighted mean
function [xRef] = wMean(x, u)

    w = 1. ./ u.^2;
    kw = 1. / sum(w);
    w = kw * w
    xRef = w * x';

endfunction


// chi-squared test
function [ok, chi2] = chi2Check(x, u)

    xRef = wMean(x, u);
    s2 = (x - xRef) .^ 2;
    u2 = u .^ 2;
    chi2 = sum(s2 ./ u2);

    ok = %F;
    n = length(x);
    if (chi2 <= CHI2CRIT(n - 1)) then
        ok = %T;
    end

endfunction


function [ok] = pairwiseConsistencyCheck(x, u)

    n = length(x);

    ok = %F;
    for i = 1 : n
        for j = (i + 1) : n
            r = 0.5 * abs(x(i) - x(j)) / sqrt(u(i)^2 + u(j)^2);
            if (r > 1.) then
                return;
            end
        end
    end

    ok = %T;

endfunction


function [k, uExt] = pairwiseConsistency(x, u)

    u2 = u .^ 2;
    n = length(x);

    A_cond = [];  b_cond = [];

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

    // objective function
    f = u2 * a_opt

    // expansion factors
    k = sqrt(a_opt);
    // rounding...
    //k = 1.e-4 * round(1.e4 * k)

    // new uncertainties
    uExt = (k .* u')';

endfunction


function [uExt, F, nAccept, chi2] = pairwiseConsistencyAndChi2(x, u0, ku, nTrials, seed)

    rand('seed', seed);

    nAccept = 0;

    F = 1.e10; // objective function
    chi2 = 0.;

    n = length(u0);
    uExt = zeros(1, n);
    kRand = ku * min(u0);

    for i = 1 : nTrials

        u = u0 + kRand * (3. + rand(1, n, "normal"));
        [chi2ok, chi2v] = chi2Check(x, u);

        if (chi2ok) then
            pcOk = pairwiseConsistencyCheck(x, u);
            if (pcOk) then
                nAccept = nAccept + 1;
                s = sum(u .^ 2);
                if (s < F) then
                    F = s;
                    uExt = u;
                    chi2 = chi2v;
                end
            end
        end

        if (modulo(i, 1000) == 0) then
            printf(">> %d trials\n", i);
        end
    end

endfunction




////////////////////////////////////////////////////////////////////////////////

pwc = pairwiseConsistencyCheck(x, u) // false

[k, u0] = pairwiseConsistency(x, u);
[chi2, chi2v] = chi2Check(x, u0) // false

[uExt, F, nAccept, chi2] = pairwiseConsistencyAndChi2(x, u0, 0.04, 1e6, 123456789)

k = uExt ./ u;

[u' k' uExt']

