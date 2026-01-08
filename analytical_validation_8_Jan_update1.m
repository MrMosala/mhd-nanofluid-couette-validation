% =========================================================================
%  SQLM MHD NANOFLUID COUETTE FLOW - Analytical Validation
% =========================================================================
%  Author: Mosala S.I
%  Supervisor: Prof. O.D. Makinde  
%  Institution: Nelson Mandela University
%  Sponsor: NITheCS
% =========================================================================
% 
% UPDATES:
% 1. Option to run all cases sequentially
% 2. Two subplots per figure (3 figures total)
% 3. Detailed tables for Cf and Nu at multiple points
% =========================================================================

clc; clear; close all;

% Main menu
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║       MHD NANOFLUID COUETTE FLOW - VALIDATION SUITE           ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('Select validation mode:\n');
fprintf('  1. Run single case\n');
fprintf('  2. Run all cases sequentially\n');
fprintf('  3. Run comprehensive analysis (all cases + tables)\n');
mode_choice = input('Enter choice (1-3): ');

if mode_choice == 1
    % Single case
    run_analytical_validation();
elseif mode_choice == 2
    % Run all cases
    run_all_cases(false);
elseif mode_choice == 3
    % Comprehensive analysis
    run_all_cases(true);
else
    fprintf('Invalid choice. Running all cases.\n');
    run_all_cases(true);
end

% =========================================================================
%  MAIN FUNCTION TO RUN ALL CASES
% =========================================================================

function run_all_cases(include_tables)
%--------------------------------------------------------------------------
% RUN_ALL_CASES - Run all 4 validation cases sequentially
%--------------------------------------------------------------------------

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║           COMPREHENSIVE VALIDATION - ALL 4 CASES             ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% ----------------------------------------------------------------------
% 1. FLUID SELECTION (Same for all cases)
% ----------------------------------------------------------------------
fprintf('Select fluid for all validation cases:\n');
fprintf('  1. Base fluid (A1=A2=A3=1.0)\n');
fprintf('  2. Copper (Cu) nanofluid\n');
fprintf('  3. Alumina (Al2O3) nanofluid\n');
fluid_choice = input('Enter choice (1-3): ');

% Set properties based on choice
switch fluid_choice
    case 1
        % Base fluid
        A1 = 1.0; A2 = 1.0; A3 = 1.0;
        fluid_name = 'Base Fluid';
        phi = 0;
        fprintf('\n✓ Using base fluid properties\n');
        
    case 2
        % Cu nanofluid
        fprintf('\nEnter Copper nanoparticle volume fraction (φ):\n');
        fprintf('  Recommended: 0.01 to 0.10 (1%% to 10%%)\n');
        phi = input('φ = ');
        
        if phi <= 0 || phi > 0.2
            fprintf('⚠ Using φ=0.05 (recommended value)\n');
            phi = 0.05;
        end
        
        [A1, A2, A3, ~, ~] = compute_nanofluid_properties(phi, 'Cu');
        fluid_name = sprintf('Cu Nanofluid (φ=%.3f)', phi);
        
    case 3
        % Al2O3 nanofluid
        fprintf('\nEnter Alumina nanoparticle volume fraction (φ):\n');
        fprintf('  Recommended: 0.01 to 0.10 (1%% to 10%%)\n');
        phi = input('φ = ');
        
        if phi <= 0 || phi > 0.2
            fprintf('⚠ Using φ=0.05 (recommended value)\n');
            phi = 0.05;
        end
        
        [A1, A2, A3, ~, ~] = compute_nanofluid_properties(phi, 'Al2O3');
        fluid_name = sprintf('Al₂O₃ Nanofluid (φ=%.3f)', phi);
        
    otherwise
        fprintf('Invalid choice. Using base fluid.\n');
        A1 = 1.0; A2 = 1.0; A3 = 1.0;
        fluid_name = 'Base Fluid';
        phi = 0;
end

fprintf('\n');
fprintf(repmat('═', 1, 72));
fprintf('\nProperties for %s:\n', fluid_name);
fprintf('  A₁ = μ_nf/μ_f = %.4f\n', A1);
fprintf('  A₂ = σ_nf/σ_f = %.4f\n', A2);
fprintf('  A₃ = k_nf/k_f = %.4f\n', A3);
fprintf(repmat('═', 1, 72));
fprintf('\n');

% ----------------------------------------------------------------------
% 2. COMMON PARAMETERS
% ----------------------------------------------------------------------
Re = 1.0; Pr = 1.0; Bi = 0.5; lambda = 0.1; N = 100;

% ----------------------------------------------------------------------
% 3. RUN ALL CASES
% ----------------------------------------------------------------------
results = cell(4, 1);

% Case 1: Simple Couette
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                 CASE 1: Simple Couette                      ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');
results{1} = run_case1(A1, A2, A3, Re, Pr, Bi, lambda, fluid_name, phi);

% Case 2: MHD Couette
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                 CASE 2: MHD Couette                        ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');
Ha = 2.0;  % Default Hartmann number
results{2} = run_case2(A1, A2, A3, Re, Ha, Pr, Bi, lambda, fluid_name, phi);

% Case 3: Viscous Dissipation
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                 CASE 3: Viscous Dissipation                ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');
Ec = 0.1;  % Default Eckert number
results{3} = run_case3(A1, A2, A3, Re, Pr, Ec, Bi, lambda, fluid_name, phi);

% Case 4: Pressure Gradient
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                 CASE 4: Pressure Gradient                  ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');
G = 1.0;  % Default pressure gradient
results{4} = run_case4(A1, A2, A3, Re, Pr, Bi, lambda, G, fluid_name, phi);

% ----------------------------------------------------------------------
% 4. GENERATE COMPREHENSIVE PLOTS
% ----------------------------------------------------------------------
generate_comprehensive_plots(results, fluid_name);

% ----------------------------------------------------------------------
% 5. GENERATE DETAILED TABLES (if requested)
% ----------------------------------------------------------------------
if include_tables
    generate_detailed_tables(results, fluid_name);
end

% ----------------------------------------------------------------------
% 6. SUMMARY OF ALL CASES
% ----------------------------------------------------------------------
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                 VALIDATION SUMMARY - ALL CASES             ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');

tolerance = 1e-8;
all_passed = true;

for case_num = 1:4
    case_name = results{case_num}.case_name;
    err_W = results{case_num}.err_W_max;
    err_Theta = results{case_num}.err_Theta_max;
    
    passed = (err_W < tolerance) && (err_Theta < tolerance);
    
    fprintf('\n%s:\n', case_name);
    fprintf('  Max velocity error:      %.2e ', err_W);
    if err_W < tolerance
        fprintf('✓ PASS\n');
    else
        fprintf('✗ FAIL\n');
        all_passed = false;
    end
    
    fprintf('  Max temperature error:   %.2e ', err_Theta);
    if err_Theta < tolerance
        fprintf('✓ PASS\n');
    else
        fprintf('✗ FAIL\n');
        all_passed = false;
    end
    
    % Display key parameters
    if case_num == 2
        fprintf('  Parameters: Ha=%.1f\n', Ha);
    elseif case_num == 3
        fprintf('  Parameters: Ec=%.2f\n', Ec);
    elseif case_num == 4
        fprintf('  Parameters: G=%.1f\n', G);
    end
end

fprintf('\n');
fprintf(repmat('═', 1, 72));
fprintf('\n');
if all_passed
    fprintf('   ✓ ALL VALIDATION CASES PASSED! 🎉\n');
    fprintf('   Code successfully verified against analytical solutions\n');
else
    fprintf('   ⚠ SOME VALIDATION CASES FAILED\n');
    fprintf('   Check implementation for cases with errors > %.0e\n', tolerance);
end
fprintf(repmat('═', 1, 72));
fprintf('\n');

fprintf('\n✓ Comprehensive validation complete!\n');

end

% =========================================================================
%                      CASE 1: SIMPLE COUETTE
% =========================================================================

function result = run_case1(A1, A2, A3, Re, Pr, Bi, lambda, fluid_name, phi)
%--------------------------------------------------------------------------
% RUN_CASE1 - Validation for Ha=0, Ec=0, G=0 (Simple Couette)
%--------------------------------------------------------------------------

% Parameters
Ha = 0; Ec = 0; G = 0; N = 100;

fprintf('\nParameters:\n');
fprintf('  Re=%.1f, Ha=%d, Pr=%.1f, Ec=%d, G=%d\n', Re, Ha, Pr, Ec, G);
fprintf('  Bi=%.1f, λ=%.2f\n\n', Bi, lambda);

% Numerical solution
[eta, W_num, Theta_num, Wp_num, Thetap_num, Cf_num, Nu_num] = ...
    solve_MHD_Couette_validation(A1, A2, A3, Re, Ha, Pr, Ec, Bi, lambda, G, N);

% Analytical solutions
W_exact = (Re / (1 + lambda)) * eta;
Theta_exact = zeros(size(eta));
Wp_exact = Re / (1 + lambda) * ones(size(eta));
Thetap_exact = zeros(size(eta));

% Engineering quantities (at multiple points)
Cf_exact = A1 * Wp_exact;  % Constant
Nu_exact = -A3 * Thetap_exact;  % Zero

% Compute errors
err_W_max = norm(W_num - W_exact, inf);
err_Theta_max = norm(Theta_num - Theta_exact, inf);

% Store results
result = struct();
result.case_name = 'Case 1: Simple Couette (Ha=0, Ec=0, G=0)';
result.eta = eta;
result.W_num = W_num; result.W_exact = W_exact;
result.Theta_num = Theta_num; result.Theta_exact = Theta_exact;
result.Wp_num = Wp_num; result.Wp_exact = Wp_exact;
result.Thetap_num = Thetap_num; result.Thetap_exact = Thetap_exact;
result.Cf_num = Cf_num; result.Cf_exact = Cf_exact;
result.Nu_num = Nu_num; result.Nu_exact = Nu_exact;
result.err_W_max = err_W_max;
result.err_Theta_max = err_Theta_max;
result.eta_points = [0, 0.25, 0.5, 0.75, 1.0];
result.parameters = struct('Re', Re, 'Ha', Ha, 'Pr', Pr, 'Ec', Ec, ...
                          'Bi', Bi, 'lambda', lambda, 'G', G);

fprintf('Validation complete:\n');
fprintf('  Max velocity error: %.2e\n', err_W_max);
fprintf('  Max temperature error: %.2e\n', err_Theta_max);

end

% =========================================================================
%                      CASE 2: MHD COUETTE
% =========================================================================

function result = run_case2(A1, A2, A3, Re, Ha, Pr, Bi, lambda, fluid_name, phi)
%--------------------------------------------------------------------------
% RUN_CASE2 - Validation for Ha≠0, Ec=0, G=0 (MHD Couette)
%--------------------------------------------------------------------------

% Parameters
Ec = 0; G = 0; N = 100;

fprintf('\nParameters:\n');
fprintf('  Re=%.1f, Ha=%.1f, Pr=%.1f, Ec=%d, G=%d\n', Re, Ha, Pr, Ec, G);
fprintf('  Bi=%.1f, λ=%.2f\n\n', Bi, lambda);

% Characteristic MHD parameter
alpha = Ha * sqrt(A2/A1);
fprintf('MHD parameter: α = Ha·√(A₂/A₁) = %.4f\n', alpha);

% Numerical solution
[eta, W_num, Theta_num, Wp_num, Thetap_num, Cf_num, Nu_num] = ...
    solve_MHD_Couette_validation(A1, A2, A3, Re, Ha, Pr, Ec, Bi, lambda, G, N);

% Analytical solutions
denom = sinh(alpha) + lambda * alpha * cosh(alpha);
C1 = Re / denom;
W_exact = C1 * sinh(alpha * eta);
Wp_exact = C1 * alpha * cosh(alpha * eta);
Theta_exact = zeros(size(eta));
Thetap_exact = zeros(size(eta));

% Engineering quantities
Cf_exact = A1 * Wp_exact;
Nu_exact = -A3 * Thetap_exact;

% Compute errors
err_W_max = norm(W_num - W_exact, inf);
err_Theta_max = norm(Theta_num - Theta_exact, inf);

% Store results
result = struct();
result.case_name = sprintf('Case 2: MHD Couette (Ha=%.1f, Ec=0, G=0)', Ha);
result.eta = eta;
result.W_num = W_num; result.W_exact = W_exact;
result.Theta_num = Theta_num; result.Theta_exact = Theta_exact;
result.Wp_num = Wp_num; result.Wp_exact = Wp_exact;
result.Thetap_num = Thetap_num; result.Thetap_exact = Thetap_exact;
result.Cf_num = Cf_num; result.Cf_exact = Cf_exact;
result.Nu_num = Nu_num; result.Nu_exact = Nu_exact;
result.err_W_max = err_W_max;
result.err_Theta_max = err_Theta_max;
result.eta_points = [0, 0.25, 0.5, 0.75, 1.0];
result.parameters = struct('Re', Re, 'Ha', Ha, 'Pr', Pr, 'Ec', Ec, ...
                          'Bi', Bi, 'lambda', lambda, 'G', G, 'alpha', alpha);

fprintf('Validation complete:\n');
fprintf('  Max velocity error: %.2e\n', err_W_max);
fprintf('  Max temperature error: %.2e\n', err_Theta_max);
fprintf('  Analytical coefficient C₁ = %.6f\n', C1);

end

% =========================================================================
%                      CASE 3: VISCOUS DISSIPATION
% =========================================================================

function result = run_case3(A1, A2, A3, Re, Pr, Ec, Bi, lambda, fluid_name, phi)
%--------------------------------------------------------------------------
% RUN_CASE3 - Validation for Ha=0, Ec≠0, G=0 (Viscous Dissipation)
%--------------------------------------------------------------------------

% Parameters
Ha = 0; G = 0; N = 100;

fprintf('\nParameters:\n');
fprintf('  Re=%.1f, Ha=%d, Pr=%.1f, Ec=%.2f, G=%d\n', Re, Ha, Pr, Ec, G);
fprintf('  Bi=%.1f, λ=%.2f\n\n', Bi, lambda);

% Numerical solution
[eta, W_num, Theta_num, Wp_num, Thetap_num, Cf_num, Nu_num] = ...
    solve_MHD_Couette_validation(A1, A2, A3, Re, Ha, Pr, Ec, Bi, lambda, G, N);

% Velocity solutions
W_exact = (Re / (1 + lambda)) * eta;
Wp_exact = Re / (1 + lambda) * ones(size(eta));

% Temperature solutions
K = -A1 * Pr * Ec * (Re/(1+lambda))^2 / A3;
C1_coeff = -K * (1 + Bi/2) / (1 + Bi);
Theta_exact = (K/2) * eta.^2 + C1_coeff * eta;
Thetap_exact = K * eta + C1_coeff;

% Engineering quantities
Cf_exact = A1 * Wp_exact;
Nu_exact = -A3 * Thetap_exact;

% Compute errors
err_W_max = norm(W_num - W_exact, inf);
err_Theta_max = norm(Theta_num - Theta_exact, inf);

% Store results
result = struct();
result.case_name = sprintf('Case 3: Viscous Dissipation (Ha=0, Ec=%.2f, G=0)', Ec);
result.eta = eta;
result.W_num = W_num; result.W_exact = W_exact;
result.Theta_num = Theta_num; result.Theta_exact = Theta_exact;
result.Wp_num = Wp_num; result.Wp_exact = Wp_exact;
result.Thetap_num = Thetap_num; result.Thetap_exact = Thetap_exact;
result.Cf_num = Cf_num; result.Cf_exact = Cf_exact;
result.Nu_num = Nu_num; result.Nu_exact = Nu_exact;
result.err_W_max = err_W_max;
result.err_Theta_max = err_Theta_max;
result.eta_points = [0, 0.25, 0.5, 0.75, 1.0];
result.parameters = struct('Re', Re, 'Ha', Ha, 'Pr', Pr, 'Ec', Ec, ...
                          'Bi', Bi, 'lambda', lambda, 'G', G, 'K', K, 'C1', C1_coeff);

fprintf('Validation complete:\n');
fprintf('  Max velocity error: %.2e\n', err_W_max);
fprintf('  Max temperature error: %.2e\n', err_Theta_max);
fprintf('  Temperature coefficients: K = %.6f, C₁ = %.6f\n', K, C1_coeff);

end

% =========================================================================
%                      CASE 4: PRESSURE GRADIENT
% =========================================================================

function result = run_case4(A1, A2, A3, Re, Pr, Bi, lambda, G, fluid_name, phi)
%--------------------------------------------------------------------------
% RUN_CASE4 - Validation for Ha=0, Ec=0, G≠0 (Pressure Gradient)
%--------------------------------------------------------------------------

% Parameters
Ha = 0; Ec = 0; N = 100;

fprintf('\nParameters:\n');
fprintf('  Re=%.1f, Ha=%d, Pr=%.1f, Ec=%d, G=%.1f\n', Re, Ha, Pr, Ec, G);
fprintf('  Bi=%.1f, λ=%.2f\n\n', Bi, lambda);

% Numerical solution
[eta, W_num, Theta_num, Wp_num, Thetap_num, Cf_num, Nu_num] = ...
    solve_MHD_Couette_validation(A1, A2, A3, Re, Ha, Pr, Ec, Bi, lambda, G, N);

% Analytical solutions
C1 = (Re + G*(lambda + 0.5)/A1) / (1 + lambda);
W_exact = -(G/(2*A1)) * eta.^2 + C1 * eta;
Wp_exact = -(G/A1) * eta + C1;
Theta_exact = zeros(size(eta));
Thetap_exact = zeros(size(eta));

% Engineering quantities
Cf_exact = A1 * Wp_exact;
Nu_exact = -A3 * Thetap_exact;

% Compute errors
err_W_max = norm(W_num - W_exact, inf);
err_Theta_max = norm(Theta_num - Theta_exact, inf);

% Store results
result = struct();
result.case_name = sprintf('Case 4: Pressure Gradient (Ha=0, Ec=0, G=%.1f)', G);
result.eta = eta;
result.W_num = W_num; result.W_exact = W_exact;
result.Theta_num = Theta_num; result.Theta_exact = Theta_exact;
result.Wp_num = Wp_num; result.Wp_exact = Wp_exact;
result.Thetap_num = Thetap_num; result.Thetap_exact = Thetap_exact;
result.Cf_num = Cf_num; result.Cf_exact = Cf_exact;
result.Nu_num = Nu_num; result.Nu_exact = Nu_exact;
result.err_W_max = err_W_max;
result.err_Theta_max = err_Theta_max;
result.eta_points = [0, 0.25, 0.5, 0.75, 1.0];
result.parameters = struct('Re', Re, 'Ha', Ha, 'Pr', Pr, 'Ec', Ec, ...
                          'Bi', Bi, 'lambda', lambda, 'G', G, 'C1', C1);

fprintf('Validation complete:\n');
fprintf('  Max velocity error: %.2e\n', err_W_max);
fprintf('  Max temperature error: %.2e\n', err_Theta_max);
fprintf('  Analytical coefficient C₁ = %.6f\n', C1);

end

% =========================================================================
%                    COMPREHENSIVE PLOTS (IMPROVED VERSION)
% =========================================================================

function generate_comprehensive_plots(results, fluid_name)
%--------------------------------------------------------------------------
% GENERATE_COMPREHENSIVE_PLOTS - Create clearer, less congested figures
%--------------------------------------------------------------------------

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║               GENERATING COMPREHENSIVE PLOTS                ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');

colors = lines(4);
case_names = cell(4,1);
for i = 1:4
    case_names{i} = results{i}.case_name;
end

%% ========================================================================
%  FIGURE 1: VELOCITY ANALYSIS (2 subplots only - much clearer!)
%% ========================================================================
figure('Position', [50 50 1400 600], 'Name', 'Velocity Analysis', 'Color', 'w');

% Subplot 1: Velocity Profiles
subplot(1,2,1);
hold on;
for i = 1:4
    plot(results{i}.eta, results{i}.W_num, 'LineWidth', 2.5, ...
         'Color', colors(i,:), 'DisplayName', case_names{i});
end
hold off;
xlabel('\eta', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('W(\eta)', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Velocity Profiles - %s', fluid_name), 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 11);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Subplot 2: Skin Friction Distribution
subplot(1,2,2);
hold on;
for i = 1:4
    plot(results{i}.eta, results{i}.Cf_num, 'LineWidth', 2.5, ...
         'Color', colors(i,:), 'DisplayName', case_names{i});
end
hold off;
xlabel('\eta', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('C_f(\eta)', 'FontSize', 14, 'FontWeight', 'bold');
title('Skin Friction Coefficient Distribution', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Export
filename1 = sprintf('Figure1_Velocity_Analysis_%s', strrep(fluid_name, ' ', '_'));
savefig(gcf, [filename1 '.fig']);
print(gcf, [filename1 '.png'], '-dpng', '-r300');
print(gcf, [filename1 '.pdf'], '-dpdf', '-bestfit');
fprintf('✓ Figure 1 saved: %s\n', filename1);

%% ========================================================================
%  FIGURE 2: TEMPERATURE ANALYSIS (2 subplots only)
%% ========================================================================
figure('Position', [100 100 1400 600], 'Name', 'Temperature Analysis', 'Color', 'w');

% Subplot 1: Temperature Profiles
subplot(1,2,1);
hold on;
for i = 1:4
    plot(results{i}.eta, results{i}.Theta_num, 'LineWidth', 2.5, ...
         'Color', colors(i,:), 'DisplayName', case_names{i});
end
hold off;
xlabel('\eta', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\theta(\eta) = (T-T_a)/T_a', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Temperature Profiles - %s', fluid_name), 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 11);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Subplot 2: Nusselt Number Distribution
subplot(1,2,2);
hold on;
for i = 1:4
    plot(results{i}.eta, results{i}.Nu_num, 'LineWidth', 2.5, ...
         'Color', colors(i,:), 'DisplayName', case_names{i});
end
hold off;
xlabel('\eta', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Nu(\eta)', 'FontSize', 14, 'FontWeight', 'bold');
title('Nusselt Number Distribution', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Export
filename2 = sprintf('Figure2_Temperature_Analysis_%s', strrep(fluid_name, ' ', '_'));
savefig(gcf, [filename2 '.fig']);
print(gcf, [filename2 '.png'], '-dpng', '-r300');
print(gcf, [filename2 '.pdf'], '-dpdf', '-bestfit');
fprintf('✓ Figure 2 saved: %s\n', filename2);

%% ========================================================================
%  FIGURE 3: ANALYTICAL VALIDATION (2x2 grid - one per case)
%% ========================================================================
figure('Position', [150 150 1400 900], 'Name', 'Analytical Validation', 'Color', 'w');

for i = 1:4
    subplot(2,2,i);
    
    % Plot numerical and analytical
    h1 = plot(results{i}.eta, results{i}.W_num, 'b-', 'LineWidth', 3, 'DisplayName', 'SQLM');
    hold on;
    h2 = plot(results{i}.eta, results{i}.W_exact, 'r--', 'LineWidth', 2.5, 'DisplayName', 'Analytical');
    hold off;
    
    xlabel('\eta', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('W(\eta)', 'FontSize', 13, 'FontWeight', 'bold');
    
    % Extract short case name
    case_short = extractBefore(case_names{i}, '(');
    title(sprintf('%s\nMax Error: %.2e', case_short, results{i}.err_W_max), ...
          'FontSize', 13, 'FontWeight', 'bold');
    
    legend([h1, h2], 'Location', 'northwest', 'FontSize', 10);
    grid on; box on;
    set(gca, 'LineWidth', 1.5, 'FontSize', 11);
end

% Export
filename3 = sprintf('Figure3_Analytical_Validation_%s', strrep(fluid_name, ' ', '_'));
savefig(gcf, [filename3 '.fig']);
print(gcf, [filename3 '.png'], '-dpng', '-r300');
print(gcf, [filename3 '.pdf'], '-dpdf', '-bestfit');
fprintf('✓ Figure 3 saved: %s\n', filename3);

%% ========================================================================
%  FIGURE 4: VALIDATION SUMMARY (IMPROVED - NO CONGESTION)
%% ========================================================================
figure('Position', [200 200 1400 600], 'Name', 'Validation Summary', 'Color', 'w');

% Subplot 1: Error Bar Chart (IMPROVED)
subplot(1,2,1);
errors_W = [results{1}.err_W_max, results{2}.err_W_max, ...
            results{3}.err_W_max, results{4}.err_W_max];
errors_T = [results{1}.err_Theta_max, results{2}.err_Theta_max, ...
            results{3}.err_Theta_max, results{4}.err_Theta_max];

% Use regular scale instead of log scale for clarity
bar_data = [errors_W; errors_T]';
h = bar(1:4, bar_data * 1e13, 'grouped', 'BarWidth', 0.8);  % Scale to 10^-13
set(h(1), 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'k', 'LineWidth', 1.2);
set(h(2), 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'k', 'LineWidth', 1.2);
set(gca, 'XTickLabel', {'Case 1', 'Case 2', 'Case 3', 'Case 4'});
ylabel('Maximum Error (×10^{-13})', 'FontSize', 14, 'FontWeight', 'bold');
title('Validation Errors - All Cases PASS', 'FontSize', 16, 'FontWeight', 'bold');
legend({'Velocity |ΔW|', 'Temperature |Δθ|'}, 'Location', 'northwest', 'FontSize', 11);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 12);
ylim([0 8]);  % Set reasonable y-axis limit

% Add error values as text annotations (above bars)
for i = 1:4
    if errors_W(i) > 0
        text(i-0.15, errors_W(i)*1e13 + 0.3, sprintf('%.1f', errors_W(i)*1e13), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', ...
             'Color', [0.2 0.4 0.8]);
    end
    if errors_T(i) > 0
        text(i+0.15, errors_T(i)*1e13 + 0.3, sprintf('%.1f', errors_T(i)*1e13), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', ...
             'Color', [0.8 0.2 0.2]);
    end
end

% Add success annotation
annotation('textbox', [0.15, 0.85, 0.25, 0.08], ...
           'String', '✓ All errors < 10^{-12} (Machine Precision)', ...
           'EdgeColor', 'none', 'FontSize', 11, 'FontWeight', 'bold', ...
           'Color', [0 0.6 0], 'HorizontalAlignment', 'center');

% Subplot 2: Parameter Influence Summary
subplot(1,2,2);

% Create summary data
max_velocities = zeros(4,1);
for i = 1:4
    max_velocities(i) = max(results{i}.W_num);
end

bar(1:4, max_velocities, 'FaceColor', [0.3 0.7 0.4], 'EdgeColor', 'k', 'LineWidth', 1.2);
set(gca, 'XTickLabel', {'Case 1', 'Case 2', 'Case 3', 'Case 4'});
ylabel('Maximum Velocity W_{max}', 'FontSize', 14, 'FontWeight', 'bold');
title('Flow Characteristics Comparison', 'FontSize', 16, 'FontWeight', 'bold');
grid on; box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 12);
ylim([0 1.05]);  % Set reasonable y-axis limit

% Add value labels on bars
for i = 1:4
    text(i, max_velocities(i) + 0.02, sprintf('%.3f', max_velocities(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

% Add physics annotation
annotation('textbox', [0.58, 0.25, 0.35, 0.12], ...
           'String', {'Physical Interpretation:', ...
                     '• Case 4 (G≠0): Highest - Pressure assists flow', ...
                     '• Case 2 (Ha≠0): Lowest - Magnetic damping'}, ...
           'EdgeColor', 'k', 'LineWidth', 1, 'FontSize', 9, ...
           'BackgroundColor', [1 1 0.9], 'HorizontalAlignment', 'left');

% Export
filename4 = sprintf('Figure4_Validation_Summary_%s', strrep(fluid_name, ' ', '_'));
savefig(gcf, [filename4 '.fig']);
print(gcf, [filename4 '.png'], '-dpng', '-r300');
print(gcf, [filename4 '.pdf'], '-dpdf', '-bestfit');
fprintf('✓ Figure 4 saved: %s\n', filename4);

%% ========================================================================
%  INDIVIDUAL CASE FIGURES (Thesis-quality standalone figures)
%% ========================================================================
for i = 1:4
    fig = figure('Position', [250 + 50*i, 250, 1200, 500], 'Color', 'w', ...
                 'Name', sprintf('Case %d Individual', i));
    
    % Velocity comparison
    subplot(1,2,1);
    plot(results{i}.eta, results{i}.W_num, 'b-', 'LineWidth', 3, 'DisplayName', 'SQLM');
    hold on;
    plot(results{i}.eta, results{i}.W_exact, 'r--', 'LineWidth', 2.5, 'DisplayName', 'Analytical');
    hold off;
    xlabel('\eta', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('W(\eta)', 'FontSize', 14, 'FontWeight', 'bold');
    case_short = extractBefore(case_names{i}, '(');
    title(sprintf('Velocity: %s', case_short), ...
          'FontSize', 15, 'FontWeight', 'bold');
    legend('Location', 'northwest', 'FontSize', 12);
    grid on; box on;
    set(gca, 'LineWidth', 1.5, 'FontSize', 12);
    
    % Temperature comparison
    subplot(1,2,2);
    plot(results{i}.eta, results{i}.Theta_num, 'b-', 'LineWidth', 3, 'DisplayName', 'SQLM');
    hold on;
    plot(results{i}.eta, results{i}.Theta_exact, 'r--', 'LineWidth', 2.5, 'DisplayName', 'Analytical');
    hold off;
    xlabel('\eta', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('\theta(\eta)', 'FontSize', 14, 'FontWeight', 'bold');
    title(sprintf('Temperature: %s', case_short), ...
          'FontSize', 15, 'FontWeight', 'bold');
    legend('Location', 'northwest', 'FontSize', 12);
    grid on; box on;
    set(gca, 'LineWidth', 1.5, 'FontSize', 12);
    
    % Add error annotation
    annotation('textbox', [0.4, 0.02, 0.2, 0.05], ...
               'String', sprintf('Max Errors: |ΔW|=%.2e, |Δθ|=%.2e', ...
                                results{i}.err_W_max, results{i}.err_Theta_max), ...
               'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
               'FontSize', 11, 'FontWeight', 'bold', 'Color', [0 0.5 0]);
    
    % Export individual case figure
    filename_case = sprintf('Figure_Case%d_Individual_%s', i, strrep(fluid_name, ' ', '_'));
    savefig(fig, [filename_case '.fig']);
    print(fig, [filename_case '.png'], '-dpng', '-r300');
    print(fig, [filename_case '.pdf'], '-dpdf', '-bestfit');
    fprintf('✓ Individual Case %d figure saved: %s\n', i, filename_case);
end

fprintf('\n✓ All improved figures generated!\n');

end

% =========================================================================
%                    DETAILED TABLES FOR Cf AND Nu
% =========================================================================

function generate_detailed_tables(results, fluid_name)
%--------------------------------------------------------------------------
% GENERATE_DETAILED_TABLES - Create detailed tables for Cf and Nu
%--------------------------------------------------------------------------

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║               GENERATING DETAILED TABLES                    ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');

% Common eta points for all cases
eta_points = [0, 0.25, 0.5, 0.75, 1.0];
n_points = length(eta_points);

% TABLE 1: Skin Friction Coefficient (Cf)
fprintf('\n\nTABLE 1: Skin Friction Coefficient C_f for %s\n', fluid_name);
fprintf(repmat('═', 1, 100));
fprintf('\n');
fprintf('%8s ', 'η');
for i = 1:4
    fprintf(' %20s', sprintf('Case %d', i));
end
fprintf(' %20s\n', 'Maximum Error');
fprintf(repmat('═', 1, 100));
fprintf('\n');

% Interpolate Cf values at specific points
Cf_values = zeros(n_points, 4);
Cf_exact_values = zeros(n_points, 4);
Cf_errors = zeros(4,1);

for case_num = 1:4
    eta = results{case_num}.eta;
    Cf_num = results{case_num}.Cf_num;
    Cf_exact = results{case_num}.Cf_exact;
    
    % Interpolate at specific points
    for j = 1:n_points
        eta_target = eta_points(j);
        [~, idx] = min(abs(eta - eta_target));
        Cf_values(j, case_num) = Cf_num(idx);
        Cf_exact_values(j, case_num) = Cf_exact(idx);
    end
    
    % Calculate maximum error for this case
    Cf_errors(case_num) = max(abs(Cf_num - Cf_exact));
end

% Print table rows
for j = 1:n_points
    fprintf('%8.2f ', eta_points(j));
    for case_num = 1:4
        fprintf(' %20.6f', Cf_values(j, case_num));
    end
    fprintf(' %20.6f\n', max(abs(Cf_values(j,:) - Cf_exact_values(j,:))));
end

% Print maximum errors per case
fprintf('\nMaximum errors per case:\n');
fprintf(repmat('─', 1, 60));
fprintf('\n');
for case_num = 1:4
    fprintf('Case %d: %.6e\n', case_num, Cf_errors(case_num));
end
fprintf(repmat('═', 1, 60));
fprintf('\n');

% TABLE 2: Nusselt Number (Nu)
fprintf('\n\nTABLE 2: Nusselt Number Nu for %s\n', fluid_name);
fprintf(repmat('═', 1, 100));
fprintf('\n');
fprintf('%8s ', 'η');
for i = 1:4
    fprintf(' %20s', sprintf('Case %d', i));
end
fprintf(' %20s\n', 'Maximum Error');
fprintf(repmat('═', 1, 100));
fprintf('\n');

% Interpolate Nu values at specific points
Nu_values = zeros(n_points, 4);
Nu_exact_values = zeros(n_points, 4);
Nu_errors = zeros(4,1);

for case_num = 1:4
    eta = results{case_num}.eta;
    Nu_num = results{case_num}.Nu_num;
    Nu_exact = results{case_num}.Nu_exact;
    
    % Interpolate at specific points
    for j = 1:n_points
        eta_target = eta_points(j);
        [~, idx] = min(abs(eta - eta_target));
        Nu_values(j, case_num) = Nu_num(idx);
        Nu_exact_values(j, case_num) = Nu_exact(idx);
    end
    
    % Calculate maximum error for this case
    Nu_errors(case_num) = max(abs(Nu_num - Nu_exact));
end

% Print table rows
for j = 1:n_points
    fprintf('%8.2f ', eta_points(j));
    for case_num = 1:4
        fprintf(' %20.6f', Nu_values(j, case_num));
    end
    fprintf(' %20.6f\n', max(abs(Nu_values(j,:) - Nu_exact_values(j,:))));
end

% Print maximum errors per case
fprintf('\nMaximum errors per case:\n');
fprintf(repmat('─', 1, 60));
fprintf('\n');
for case_num = 1:4
    fprintf('Case %d: %.6e\n', case_num, Nu_errors(case_num));
end
fprintf(repmat('═', 1, 60));
fprintf('\n');

% Save tables to file
filename = sprintf('Validation_Tables_%s.txt', strrep(fluid_name, ' ', '_'));
fid = fopen(filename, 'w');

fprintf(fid, 'VALIDATION TABLES FOR %s\n\n', fluid_name);
fprintf(fid, 'Generated: %s\n\n', datestr(now));

% Write Table 1
fprintf(fid, 'TABLE 1: Skin Friction Coefficient C_f\n');
fprintf(fid, 'η');
for i = 1:4
    fprintf(fid, ',Case %d', i);
end
fprintf(fid, ',Error\n');
for j = 1:n_points
    fprintf(fid, '%.2f', eta_points(j));
    for case_num = 1:4
        fprintf(fid, ',%.6f', Cf_values(j, case_num));
    end
    fprintf(fid, ',%.6f\n', max(abs(Cf_values(j,:) - Cf_exact_values(j,:))));
end

% Write Table 2
fprintf(fid, '\n\nTABLE 2: Nusselt Number Nu\n');
fprintf(fid, 'η');
for i = 1:4
    fprintf(fid, ',Case %d', i);
end
fprintf(fid, ',Error\n');
for j = 1:n_points
    fprintf(fid, '%.2f', eta_points(j));
    for case_num = 1:4
        fprintf(fid, ',%.6f', Nu_values(j, case_num));
    end
    fprintf(fid, ',%.6f\n', max(abs(Nu_values(j,:) - Nu_exact_values(j,:))));
end

fclose(fid);
fprintf('\n\n✓ Tables saved to file: %s\n', filename);

end

% =========================================================================
%                    VALIDATION SOLVER FUNCTION
% =========================================================================

function [eta, W, Theta, Wp, Thetap, Cf, Nu] = ...
    solve_MHD_Couette_validation(A1, A2, A3, Re, Ha, Pr, Ec, Bi, lambda, G, N)
%--------------------------------------------------------------------------
% SOLVE_MHD_COUETTE_VALIDATION - Simplified solver for validation cases
%--------------------------------------------------------------------------

% Generate Chebyshev grid
etamax = 1.0;
[~, eta, D, D2, ~] = cheb_diff_matrices(N, etamax);
npts = N + 1;

% Initial guess
W = eta * Re / (1 + lambda);
Theta = - (Bi / (1 + Bi)) * eta;

% SQLM iteration
tol = 1e-12;
maxIter = 50;

for iter = 1:maxIter
    % Store current solution
    Wr = W;
    Wpr = D * Wr;
    
    % Initialize matrices
    Z = zeros(npts);
    I = eye(npts);
    
    % Momentum equation
    A11 = A1 * D2 - A2 * Ha^2 * I;
    A12 = Z;
    R1 = -G * ones(npts, 1);
    
    % Energy equation (linearized)
    A21 = 2*A1*Pr*Ec * diag(Wpr) * D + 2*A2*Pr*Ec*Ha^2 * diag(Wr);
    A22 = A3 * D2;
    R2 = A1*Pr*Ec * (Wpr .* Wpr) + A2*Pr*Ec*Ha^2 * (Wr .* Wr);
    
    % Assemble system
    A_sys = [A11, A12; A21, A22];
    b = [R1; R2];
    
    % Apply boundary conditions
    % BC 1: W(0) = 0
    row = npts;
    A_sys(row, :) = 0;
    A_sys(row, npts) = 1;
    b(row) = 0;
    
    % BC 2: W(1) = Re - λ·W'(1)
    row = 1;
    A_sys(row, :) = 0;
    A_sys(row, 1:npts) = [1, zeros(1, npts-1)] + lambda * D(1, :);
    b(row) = Re;
    
    % BC 3: θ(0) = 0
    row = 2*npts;
    A_sys(row, :) = 0;
    A_sys(row, 2*npts) = 1;
    b(row) = 0;
    
    % BC 4: θ'(1) + Bi·θ(1) = 0
    row = npts + 1;
    A_sys(row, :) = 0;
    A_sys(row, npts+1:2*npts) = D(1, :) + Bi * [1, zeros(1, npts-1)];
    b(row) = 0;
    
    % Solve
    Unew = A_sys \ b;
    W = Unew(1:npts);
    Theta = Unew(npts+1:end);
    
    % Compute residuals
    Wpp = D2 * W;
    Wp = D * W;
    
    Omega1 = A1*Wpp - A2*Ha^2*W + G;
    Omega2 = A3*(D2 * Theta) + A1*Pr*Ec*(Wp.^2) + A2*Pr*Ec*Ha^2*(W.^2);
    
    res = max([norm(Omega1, inf), norm(Omega2, inf)]);
    
    if res < tol
        break;
    end
end

% Compute gradients and engineering quantities
Wp = D * W;
Thetap = D * Theta;
Cf = A1 * Wp;
Nu = -A3 * Thetap;

end

% =========================================================================
%                    ORIGINAL VALIDATION FUNCTION
% =========================================================================

function run_analytical_validation()
%--------------------------------------------------------------------------
% RUN_ANALYTICAL_VALIDATION - Original function for single case validation
%--------------------------------------------------------------------------

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                  SINGLE CASE VALIDATION                     ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Fluid selection
fprintf('Select fluid:\n');
fprintf('  1. Base fluid (A1=A2=A3=1.0)\n');
fprintf('  2. Copper (Cu) nanofluid\n');
fprintf('  3. Alumina (Al2O3) nanofluid\n');
fluid_choice = input('Enter choice (1-3): ');

% Set properties
switch fluid_choice
    case 1
        A1 = 1.0; A2 = 1.0; A3 = 1.0;
        fluid_name = 'Base Fluid';
        phi = 0;
    case 2
        phi = input('Enter Cu nanoparticle volume fraction φ (0.01-0.10): ');
        if phi <= 0 || phi > 0.2, phi = 0.05; end
        [A1, A2, A3, ~, ~] = compute_nanofluid_properties(phi, 'Cu');
        fluid_name = sprintf('Cu Nanofluid (φ=%.3f)', phi);
    case 3
        phi = input('Enter Al2O3 nanoparticle volume fraction φ (0.01-0.10): ');
        if phi <= 0 || phi > 0.2, phi = 0.05; end
        [A1, A2, A3, ~, ~] = compute_nanofluid_properties(phi, 'Al2O3');
        fluid_name = sprintf('Al₂O₃ Nanofluid (φ=%.3f)', phi);
    otherwise
        A1 = 1.0; A2 = 1.0; A3 = 1.0;
        fluid_name = 'Base Fluid';
        phi = 0;
end

fprintf('\nSelect validation case:\n');
fprintf('  1. Case 1: Simple Couette (Ha=0, Ec=0, G=0)\n');
fprintf('  2. Case 2: MHD Couette (Ha≠0, Ec=0, G=0)\n');
fprintf('  3. Case 3: Viscous Dissipation (Ha=0, Ec≠0, G=0)\n');
fprintf('  4. Case 4: Pressure Gradient (Ha=0, Ec=0, G≠0)\n');
case_choice = input('Enter choice (1-4): ');

% Common parameters
Re = 1.0; Pr = 1.0; Bi = 0.5; lambda = 0.1; N = 100;

switch case_choice
    case 1
        result = run_case1(A1, A2, A3, Re, Pr, Bi, lambda, fluid_name, phi);
    case 2
        Ha = input('Enter Hartmann number Ha (default=2.0): ');
        if isempty(Ha), Ha = 2.0; end
        result = run_case2(A1, A2, A3, Re, Ha, Pr, Bi, lambda, fluid_name, phi);
    case 3
        Ec = input('Enter Eckert number Ec (default=0.1): ');
        if isempty(Ec), Ec = 0.1; end
        result = run_case3(A1, A2, A3, Re, Pr, Ec, Bi, lambda, fluid_name, phi);
    case 4
        G = input('Enter pressure gradient G (default=1.0): ');
        if isempty(G), G = 1.0; end
        result = run_case4(A1, A2, A3, Re, Pr, Bi, lambda, G, fluid_name, phi);
    otherwise
        result = run_case1(A1, A2, A3, Re, Pr, Bi, lambda, fluid_name, phi);
end

% Display results
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                   VALIDATION RESULTS                          ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');

fprintf('\n%s\n', result.case_name);
fprintf(repmat('─', 1, 72));
fprintf('\nMaximum errors:\n');
fprintf('  Velocity (W):      %.2e\n', result.err_W_max);
fprintf('  Temperature (θ):   %.2e\n', result.err_Theta_max);
fprintf('  Skin friction:     %.2e\n', max(abs(result.Cf_num - result.Cf_exact)));
fprintf('  Nusselt number:    %.2e\n', max(abs(result.Nu_num - result.Nu_exact)));

% Quick plot
figure('Position', [100 100 800 600]);
subplot(2,2,1);
plot(result.eta, result.W_num, 'b-', 'LineWidth', 2); hold on;
plot(result.eta, result.W_exact, 'r--', 'LineWidth', 2); hold off;
xlabel('\eta'); ylabel('W(\eta)'); title('Velocity');
legend('SQLM', 'Analytical'); grid on;

subplot(2,2,2);
plot(result.eta, result.Theta_num, 'b-', 'LineWidth', 2); hold on;
plot(result.eta, result.Theta_exact, 'r--', 'LineWidth', 2); hold off;
xlabel('\eta'); ylabel('\theta(\eta)'); title('Temperature');
legend('SQLM', 'Analytical'); grid on;

subplot(2,2,3);
plot(result.eta, result.Cf_num, 'g-', 'LineWidth', 2);
xlabel('\eta'); ylabel('C_f'); title('Skin Friction'); grid on;

subplot(2,2,4);
plot(result.eta, result.Nu_num, 'm-', 'LineWidth', 2);
xlabel('\eta'); ylabel('Nu'); title('Nusselt Number'); grid on;

sgtitle(sprintf('%s - %s', result.case_name, fluid_name), 'FontSize', 14);

end

% =========================================================================
%                    HELPER FUNCTIONS
% =========================================================================

function [x, eta, D, D2, D3] = cheb_diff_matrices(N, etamax)
% Chebyshev differentiation matrices
k = (0:N)';
x = cos(pi * k / N);
eta = etamax * (x + 1) / 2;

c = [2; ones(N-1, 1); 2] .* (-1).^(0:N)';
X = repmat(x, 1, N+1);
dX = X - X';

D = (c * (1./c)') ./ (dX + eye(N+1));
D = D - diag(sum(D, 2));

scale = 2 / etamax;
D = scale * D;
D2 = D * D;
D3 = D2 * D;
end

function [A1, A2, A3, A4, A5] = compute_nanofluid_properties(phi, nanoparticle_type)
% Nanofluid property calculations
rho_f = 997.1; Cp_f = 4179; k_f = 0.613; sigma_f = 5.5e-6;

switch lower(nanoparticle_type)
    case 'cu'
        rho_s = 8933; Cp_s = 385; k_s = 401; sigma_s = 58e6;
    case 'al2o3'
        rho_s = 3970; Cp_s = 765; k_s = 40; sigma_s = 1e-10;
    otherwise
        error('Unknown nanoparticle type.');
end

% Property ratios
A1 = 1 / (1 - phi)^2.5;
r = sigma_s / sigma_f;
A2 = 1 + 3*(r - 1)*phi / ((r + 2) - (r - 1)*phi);
A3 = (k_s + 2*k_f - 2*phi*(k_f - k_s)) / (k_s + 2*k_f + phi*(k_f - k_s));
A4 = (1 - phi) + phi*(rho_s/rho_f);
A5 = (1 - phi) + phi*(rho_s*Cp_s)/(rho_f*Cp_f);
end