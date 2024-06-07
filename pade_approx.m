function approx = pade_approx(data, p, q)
    % Padé approximation function
    % Calculate Padé approximation coefficients
    c = zeros(p + q + 1, 1);
    c(1) = 1;
    for k = 1:p
        c(k + 1) = c(k) * (p - k + 1) / (k * (q + k));
    end
    
    % Calculate Padé approximation
    n = length(data);
    approx = zeros(n, 1);
    for i = 1:n
        approx(i) = 0;
        for j = 1:min(i, p + 1)
            approx(i) = approx(i) + c(j) * data(i - j + 1);
        end
    end
end
