function [samples, expected_errs] = ...
    herding_sequential_sample( mix, kernel, num_samples, range, num_queries)
%
% Sequentially selects sample locations with respect to a base distribution in
% order to minimize the expected variance of an integral
% using Kernel Herding.
%
% mix specifies a mixture of Gaussians.

[K,D] = size(mix.means);

samples = NaN(0, D);
expected_errs = NaN(0, D);

for i = 1:num_samples
    [samples(i, :), expected_errs(i)] = ...
        herding_next_sample( mix, kernel, samples, range, num_queries );
    fprintf('h');
end
