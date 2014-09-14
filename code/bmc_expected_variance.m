function [expected_variance, weights] = ...
    bmc_expected_variance( mix, kernel, existing_samples )
%
% Compute the expected variance of the estimate of Z given a set of existing
% samples.
%
% David Duvenaud
% March 2012

[K,D] = size(mix.means);
num_existing = size(existing_samples, 1);

% Compute pre-existing z's.
existing_zs = zeros( num_existing, 1);
for k = 1:K
    existing_zs = existing_zs + mix.weights(k) ...
        * mvnpdf( existing_samples, mix.means(k, :), ...
                  mix.covs(:, :, k) + kernel.covariance );
end
existing_zs = existing_zs .* kernel.height;

% Compute Gram matrix at existing sample locations.
existing_K = NaN(num_existing, num_existing);
for i = 1:num_existing
    existing_K(i, :) = kernel.height ...
                     .*  mvnpdf( existing_samples, existing_samples(i,:), ...
                                 kernel.covariance );
end

% Improve conditioning.
existing_K = existing_K + diag(ones(num_existing, 1) .* 1e-6 .* max(max(existing_K)));

zKz = existing_zs' * (existing_K \ existing_zs);


prior_variance = prior_variance_mix( mix, kernel );
expected_variance = prior_variance - zKz;

assert( all( expected_variance >= 0 ));

weights = existing_zs' / existing_K;
