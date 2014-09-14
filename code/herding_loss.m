function expected_loss = herding_loss( mix, kernel, existing_samples )
%
% Compute the expected variance of the estimate of Z given a set of existing
% samples, and a proposed new one.
%
%
% David Duvenaud
% March 2012

prior_variance = prior_variance_mix( mix, kernel );

[K, D] = size(mix.means);
num_existing = size(existing_samples, 1);

% Compute pre-existing z's.
existing_zs = zeros( num_existing, 1);
for k = 1:K
    existing_zs = existing_zs + mix.weights(k) ...
        * mvnpdf( existing_samples, mix.means(k, :), ...
                  mix.covs(:, :, k) + kernel.covariance );
end
existing_zs = existing_zs .* kernel.height;

  % All weights the same for now.
old_weights = ones(num_existing, 1) ./ (num_existing);

% Compute Gram matrix at existing sample locations.
existing_K = NaN(num_existing, num_existing);
for i = 1:num_existing
    existing_K(i, :) = kernel.height ...
                     .*  mvnpdf( existing_samples, existing_samples(i,:), ...
                                 kernel.covariance );
end

expected_loss = prior_variance - 2*old_weights'*existing_zs ...
           + old_weights' * existing_K * old_weights;

assert( all( expected_loss >= 0 ));
