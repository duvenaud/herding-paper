function expected_loss = ...
    herding_loss_next( mix, kernel, existing_samples, sample_locations )
%
% Compute the expected variance of the estimate of Z given a set of existing
% samples, and a proposed new one.
%
%
% David Duvenaud
% March 2012

prior_variance = prior_variance_mix( mix, kernel );

[K,D] = size(mix.means);
num_existing = size(existing_samples, 1);
n_sample_locs = size(sample_locations, 1);

% Compute pre-existing z's.
existing_zs = zeros( num_existing, 1);
for k = 1:K
    existing_zs = existing_zs + mix.weights(k) ...
        * mvnpdf( existing_samples, mix.means(k, :), ...
                  mix.covs(:, :, k) + kernel.covariance );
end
existing_zs = existing_zs .* kernel.height;

% Compute new z's for points we might add.
new_zs = zeros( n_sample_locs, 1);
for k = 1:K
    new_zs = new_zs + mix.weights(k) ...
                    * mvnpdf( sample_locations, mix.means(k, :), ...
                              mix.covs(:, :, k) + kernel.covariance );
end
new_zs = new_zs .* kernel.height;

  % All weights the same for now.
old_weights = ones(num_existing, 1) ./ (num_existing + 1);
new_weight = 1 / (num_existing + 1);

% Compute Gram matrix at existing sample locations.
existing_K = NaN(num_existing, num_existing);
for i = 1:num_existing
    existing_K(i, :) = kernel.height ...
                     .*  mvnpdf( existing_samples, existing_samples(i,:), ...
                                 kernel.covariance );
end

old_term = prior_variance - 2*old_weights'*existing_zs ...
           + old_weights' * existing_K * old_weights;


% Now, for each new sample location, see what the expected variance would be if
% we added that sample.
new_term = NaN(n_sample_locs, 1);
for i = 1:n_sample_locs
    new_gram_col = kernel.height ...
                     .*  mvnpdf( existing_samples, sample_locations(i,:), ...
                                 kernel.covariance );
    new_gram_diag = kernel.height ...
                     .*  mvnpdf( sample_locations(i,:), sample_locations(i,:), ...
                                 kernel.covariance );
      
    new_term(i) = -2 * new_weight * new_zs(i) ...
                  + new_weight * 2 * old_weights' * new_gram_col ...
                  + new_weight^2 * new_gram_diag;
end

expected_loss = old_term + new_term;

assert( all( expected_loss >= 0 ));
