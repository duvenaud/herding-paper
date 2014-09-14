% Quickly produces examples of most of the figures in the paper.
%
% David Duvenaud
% March 2012

function demo

close all;

% Demo parameters:
sigma=0.5;          % Gaussian kernel width.
num_samples = 40;   % Number of samples in total.
num_queries = 1000; % Number of points to consider at each step.
                    % In the paper, we used num_queries = 10000.

% Fix the seed of the random generators.
seed=0;
randn('state',seed);
rand('state',seed);

% Load the mixture of Gaussians used in the super-samples by kernel herding
% paper (many thanks to Yutian!).
load yutian_mixture.mat obj X Y f

% Convert mixture of Gaussians to a different structure.
mix.weights = obj.PComponents';
mix.means = obj.mu;
mix.covs = obj.Sigma;

% This demo searches for possible next locations by drawing from a uniform prior
% with the following range.  This is a bad idea in high dimensions,
% and is only done here so that all the code will be really simple.
range = [ -6, 6; -5 3];

% Perform kernel herding and then BQ ( aka BMC) to choose the sample locations.
kernel.height = 1;
kernel.covariance = [sigma 0; 0 sigma].^2;

fprintf('\nComputing %d herding samples...\n', num_samples );
[herding_samples, herding_errors] = ...
    herding_sequential_sample( mix, kernel, num_samples, range, num_queries);

fprintf('\nComputing %d BQ samples...\n', num_samples );
[bmc_samples, bmc_variances] = ...
    bmc_sequential_sample( mix, kernel, num_samples, range, num_queries);



% Plot the input distribution.
% ==============================================
N_1d = 600;
xrange = linspace( range(1,1), range(1,2), N_1d);   % Choose a set of x locations.
yrange = linspace( range(1,1), range(1,2), N_1d);   % Choose a set of x locations.
[xvals, yvals] = meshgrid( xrange, yrange);
gridvals = [xvals(:) yvals(:)];
figure(1); clf;
dist_vals = mix_gaussians_pdf(gridvals, mix );
dh = contour( xvals, yvals, reshape(dist_vals, N_1d, N_1d ), 30, ...
    'LineWidth', .5); hold on;
drawnow;

% Plot herding samples.
linewidth = 1;
hold on; hsh = plot(herding_samples(1:20,1), herding_samples(1:20,2), 'rs',...
                    'LineWidth', linewidth);

% Plot the first 8 BMC samples.
% First, compute the weight of each sample.
[ev, weights] = bmc_expected_variance( mix, kernel, bmc_samples(1:8, :) );
for i = 1:8
    % Change size of marker depending on weight.
    bsh = plot(bmc_samples(i,1), bmc_samples(i,2), ...
        'kx', 'LineWidth', linewidth, 'Markersize', (weights(i).*230));
end

legend( [ hsh, bsh], {'Herding Samples', 'SBQ Samples'},... );
        'Fontsize', 10, 'Interpreter','latex')
legend boxoff

% Make the figure nice-looking.
set( gca, 'XTick', [] );
set( gca, 'yTick', [] );
set( gca, 'XTickLabel', '' );
set( gca, 'yTickLabel', '' );
set(gcf, 'color', 'white');
set(gca, 'YGrid', 'off');



% Plot weight histogram
% ===================================
fontsize = 10;
[ev, weights] = bmc_expected_variance( mix, kernel, bmc_samples );
figure(2); clf; hist(weights, 20); hold on;
    
y_limits = ylim;
onn_h = plot( [1/i, 1/i], [y_limits(1), y_limits(2)], 'r--' ); hold on;    
z_h = plot( [0, 0], [y_limits(1), y_limits(2)], 'k--' ); hold on;    
title(sprintf('weights at %d samples', num_samples));
xlabel('SBQ weight');
ylabel('count');
set(get(gca,'XLabel'),'Rotation',0,'Interpreter','latex', 'Fontsize', fontsize);
set(get(gca,'YLabel'),'Rotation',90,'Interpreter','latex', 'Fontsize', fontsize);
set(gca,'Fontsize', fontsize - 2);

legend( [ onn_h, z_h], {'1/n weight', 'zero weight'}, ...
    'Fontsize', 10, 'Interpreter','latex', 'Location', 'Best')
legend boxoff 




% Plot sum of weights
% ==================================

fprintf('\n computing weight sums... \n');
for i = 1:num_samples
    [ev, cur_weights] = ...
        bmc_expected_variance( mix, kernel, herding_samples(1:i, :) );
    weight_sum(i) = sum(cur_weights);
    fprintf('.');
end

figure(3); clf;
bqw_h = plot(weight_sum, 'b-'); hold on;
x_limits = xlim;
onn_h = plot( [x_limits(1), x_limits(2)], [1 1], 'k--' ); hold on;    
ylim( [ 0 1.1] );
xlabel( 'number of samples' );
ylabel( 'sum of weights' );
set(get(gca,'XLabel'),'Rotation',0,'Interpreter','latex', 'Fontsize', fontsize);
set(get(gca,'YLabel'),'Rotation',90,'Interpreter','latex', 'Fontsize', fontsize);
set(gca,'Fontsize', fontsize - 2 );

legend( [ onn_h, bqw_h], {'sum of herding weights', 'sum of BQ weights'}, ...
    'Fontsize', 10, 'Interpreter','latex', 'Location', 'SouthEast')
legend boxoff 



% Plot expected error cuves
% ==================================

% Compute the MMD of using BQ weights on herding samples.
fprintf('\nComputing expected error');
for i = 1:num_samples
    expected_variances_herding_points(i) = ...
        bmc_expected_variance( mix, kernel, herding_samples(1:i, :) );
end

figure(4); clf;
heh = loglog( herding_errors, 'b-' ); hold on;
hsv = loglog( expected_variances_herding_points, 'r-' ); hold on;
bvh = loglog( bmc_variances, 'g-' ); hold on;
legend( [ heh, hsv, bvh], {'Herding with 1/N weights', 'Herding with BQ weights', 'SBQ with BQ weights' }, ...
        'Fontsize', 10, 'Interpreter','latex')
legend boxoff  

% Make the plot pretty.
xlabel( 'number of samples' );
ylabel( 'MMD or $\epsilon^{2}_{BQ}$' );
set(get(gca,'XLabel'),'Rotation',0,'Interpreter','latex', 'Fontsize', fontsize);
set(get(gca,'YLabel'),'Rotation',90,'Interpreter','latex', 'Fontsize', fontsize);
set(gcf, 'color', 'white');
set(gca, 'YGrid', 'off');




% Plot in-model error
% ==================================

fprintf('\nComputing in-model empirical error\n');
num_funcs = 50;

for f = 1:num_funcs
    % Generate a function with RHKS norm of 1.
    rkhs_func.weights = randn(10,1);
    rkhs_func.means = mvnrnd([0,0], [3 0; 0 3], 10);
    rkhs_func.covs = repmat(kernel.covariance, [1, 1, 10]);
    
    % renormalize
    c = 0;
    for i = 1:length(rkhs_func.weights)
        for j = 1:length(rkhs_func.weights)
           c = c + rkhs_func.weights(i) * rkhs_func.weights(j) * mvnpdf(rkhs_func.means(i,:), rkhs_func.means(j,:), kernel.covariance);
        end
    end
    rkhs_func.weights = rkhs_func.weights ./ sqrt(c);

    cur_f = @(x) mix_gaussians_pdf(x, rkhs_func );
    cur_truth = 0;
    for i = 1:length(mix.weights)
        for j = 1:length(rkhs_func.weights)
            cur_truth = cur_truth + mix.weights(i) * rkhs_func.weights(j) ...
                * mvnpdf( mix.means(i, :), rkhs_func.means(j, :), mix.covs(:,:,i) + rkhs_func.covs(:,:,j));
        end
    end

    funcs{f} = cur_f;
    f_truths(f) = cur_truth;
end

random_draws = mix_gaussians_draw( mix, num_samples );

% Compute errors.
sample_sets = cell(0);
sample_sets{1} = random_draws;
sample_sets{2} = herding_samples;
sample_sets{3} = herding_samples;
sample_sets{4} = bmc_samples;
num_methods = length(sample_sets);

f_vals = NaN(num_samples, num_methods, num_funcs);
prediction = NaN(num_samples, num_methods, num_funcs);
truth = NaN(num_samples, num_methods, num_funcs);
errors = NaN(num_samples, num_methods, num_funcs);

for i = 1:num_samples
    for m = 1:num_methods
        cur_samples = sample_sets{m};

        % Evaluate weights.
        if m < 3
            % Herding and random weights.
            cur_weights = repmat(1/i, i, 1);
        else
            % BMC weights
            [expected_variances, cur_weights] = ...
                bmc_expected_variance( mix, kernel, cur_samples(1:i, :) );
        end

        for f = 1:num_funcs
            cur_func = funcs{f};
            
            % Evaluate samples.
            f_vals(i, m, f) = cur_func( cur_samples(i, :) );
                        
            % Make prediction (linear combination).
            prediction( i, m, f ) = cur_weights(:)' * f_vals(1:i, m, f);
            truth(i, m, f) = f_truths(f);
            errors( i, m, f ) = abs(prediction( i, m, f ) - truth(i, m, f));
        end
    end
    fprintf('.');
end

figure(5); clf;
plotrange = 1:num_samples;
hre = loglog( mean(errors(plotrange, 1, :), 3), 'k-' ); hold on;
heh = loglog( mean(errors(plotrange, 2, :), 3), 'b-' ); hold on;
hsv = loglog( mean(errors(plotrange, 3, :), 3), 'r-' ); hold on;
bvh = loglog( mean(errors(plotrange, 4, :), 3), 'g-' ); hold on;

xlabel( 'number of samples' );
ylabel( 'mean absolute error' );
set(get(gca,'XLabel'),'Rotation',0,'Interpreter','latex', 'Fontsize', fontsize);
set(get(gca,'YLabel'),'Rotation',90,'Interpreter','latex', 'Fontsize', fontsize);
xlim([1 num_samples])
set(gcf, 'color', 'white');
set(gca, 'YGrid', 'off');
legend( [ hre, heh, hsv, bvh], { 'i.i.d. sampling', 'Herding with $1/N$ weights', 'Herding with BQ weights', 'SBQ with BQ weights' }, ...
        'Fontsize', 10, 'Interpreter','latex', 'Location', 'SouthWest')
legend boxoff  
