function BlasiusPredictorCorrector()
    % Define the domain
    eta_max = 10;
    num_points = 1000;
    eta = linspace(0, eta_max, num_points);
    
    % Solve using predictor-corrector method for different schemes
    [f_4, f_prime_4, f_double_prime_4] = blasius_predictor_corrector(eta, 'pade4');
    [f_6, f_prime_6, f_double_prime_6] = blasius_predictor_corrector(eta, 'pade6');
    [f_44, f_prime_44, f_double_prime_44] = blasius_predictor_corrector(eta, 'pade44');
    
    % Debug: Check values at specific points
    disp('Sample values for f_prime_4:');
    disp(f_prime_4(1:10)); % Display the first 10 values
    disp('Sample values for f_prime_6:');
    disp(f_prime_6(1:10)); % Display the first 10 values
    disp('Sample values for f_prime_44:');
    disp(f_prime_44(1:10)); % Display the first 10 values
    
    % Plot the Blasius boundary layer profiles
    figure;
    
    % Plot f'(η)
    subplot(3, 1, 1);
    plot(eta, f_prime_4, 'b', 'DisplayName', 'Padé 4');
    hold on;
    plot(eta, f_prime_6, 'r', 'DisplayName', 'Padé 6');
    plot(eta, f_prime_44, 'g', 'DisplayName', 'Padé [4/4]');
    title('Velocity Profile: f''(\eta)');
    xlabel('\eta');
    ylabel('f''(\eta)');
    legend show;
    grid on;
    
    % Plot f(η)
    subplot(3, 1, 2);
    plot(eta, f_4, 'b', 'DisplayName', 'Padé 4');
    hold on;
    plot(eta, f_6, 'r', 'DisplayName', 'Padé 6');
    plot(eta, f_44, 'g', 'DisplayName', 'Padé [4/4]');
    title('Stream Function: f(\eta)');
    xlabel('\eta');
    ylabel('f(\eta)');
    legend show;
    grid on;
    
    % Plot f''(η)
    subplot(3, 1, 3);
    plot(eta, f_double_prime_4, 'b', 'DisplayName', 'Padé 4');
    hold on;
    plot(eta, f_double_prime_6, 'r', 'DisplayName', 'Padé 6');
    plot(eta, f_double_prime_44, 'g', 'DisplayName', 'Padé [4/4]');
    title('Shear Stress: f''''(\eta)');
    xlabel('\eta');
    ylabel('f''''(\eta)');
    legend show;
    grid on;
end
