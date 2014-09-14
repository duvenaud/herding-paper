function [sample, expected_variance] = ...
    bmc_next_sample( mix, kernel, existing_samples, range, num_queries )
%
% Chooses the next sample location, given a set of existing samples,
% and a kernel.
%
% range:  the range over which to search.
%    dx: the granularity of the search.
%

[K,D] = size(mix.means);

% Generate random sample locations in the box specified by range.
sample_locations = NaN( num_queries, D );
for d = 1:D
    sample_locations(:, d) = unifrnd( range(d, 1), range(d, 2), num_queries, 1);
end

% Evaluate the expected variance at each point.
expected_variances = bmc_expected_variance_next( mix, kernel, existing_samples, sample_locations );

% Choose the best one.
[expected_variance, best_ix] = min( expected_variances );
sample = sample_locations( best_ix, : );
