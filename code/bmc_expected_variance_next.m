function expected_variances = ...
    bmc_expected_variance_next( mix, kernel, existing_samples, sample_locations )
%
% Compute the expected variance of the estimate of Z given a set of existing
% samples, and a proposed new one.
%
%
% David Duvenaud
% March 2012


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


% Compute Gram matrix at existing sample locations.
existing_K = NaN(num_existing, num_existing);
for i = 1:num_existing
    existing_K(i, :) = kernel.height ...
                     .*  mvnpdf( existing_samples, existing_samples(i,:), ...
                                 kernel.covariance );
end


% Now, for each new sample location, see what the expected variance would be if
% we added that sample.
zKz = NaN(n_sample_locs, 1);
for i = 1:n_sample_locs
    new_gram_row = kernel.height ...
                     .*  mvnpdf( existing_samples, sample_locations(i,:), ...
                                 kernel.covariance );
    new_gram_diag = kernel.height ...
                     .*  mvnpdf( sample_locations(i,:), sample_locations(i,:), ...
                                 kernel.covariance );
    
    % Slow way.  Todo:  use fast updates.
    new_K = [ [existing_K, new_gram_row]; [new_gram_row', new_gram_diag] ];
    %new_chol = chol( new_K );
    
    % Improve conditioning.
    new_K = new_K + diag(ones(num_existing + 1, 1) .* 1e-6 .* max(max(new_K)));
    
    new_z_vec = [existing_zs; new_zs(i)];
    
    zKz(i) = new_z_vec' * (new_K \ new_z_vec);
end

prior_variance = prior_variance_mix( mix, kernel );
expected_variances = prior_variance - zKz;

assert( all( expected_variances >= 0 ));
