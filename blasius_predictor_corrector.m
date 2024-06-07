function [f, f_prime, f_double_prime] = blasius_predictor_corrector(eta, method)
    % Initial conditions
    f0 = 0;
    f_prime0 = 0;
    f_double_prime0 = 0.3320573362151963; % The Blasius value at eta = 0
    
    % Step size
    h = eta(2) - eta(1);
    
    % Initial values
    f = zeros(size(eta));
    f_prime = zeros(size(eta));
    f_double_prime = zeros(size(eta));
    f(1) = f0;
    f_prime(1) = f_prime0;
    f_double_prime(1) = f_double_prime0;
    
    % Predictor-Corrector loop
    for i = 1:length(eta) - 1
        % Predictor step (Euler's method)
        f_pred = f(i) + h * f_prime(i);
        f_prime_pred = f_prime(i) + h * f_double_prime(i);
        f_double_prime_pred = f_double_prime(i) + h * (-0.5 * f(i) * f_double_prime(i));
        
        % Corrector step (Trapezoidal rule)
        f(i + 1) = f(i) + (h / 2) * (f_prime(i) + f_prime_pred);
        f_prime(i + 1) = f_prime(i) + (h / 2) * (f_double_prime(i) + f_double_prime_pred);
        f_double_prime(i + 1) = f_double_prime(i) + (h / 2) * (-0.5 * (f(i) * f_double_prime(i) + f_pred * f_double_prime_pred));
    end
    
    % Apply the specific Padé approximation method
    switch method
        case 'pade4'
            f = pade_approx(f, 2, 2);  % Applying Padé [2/2] approximant
            f_prime = pade_approx(f_prime, 2, 2);
            f_double_prime = pade_approx(f_double_prime, 2, 2);
        case 'pade6'
            f = pade_approx(f, 3, 3);  % Applying Padé [3/3] approximant
            f_prime = pade_approx(f_prime, 3, 3);
            f_double_prime = pade_approx(f_double_prime, 3, 3);
        case 'pade44'
            f = pade_approx(f, 4, 4);  % Applying Padé [4/4] approximant
            f_prime = pade_approx(f_prime, 4, 4);
            f_double_prime = pade_approx(f_double_prime, 4, 4);
        otherwise
            error('Unknown method: %s', method);
    end
end
