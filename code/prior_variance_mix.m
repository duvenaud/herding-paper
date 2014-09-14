function prior_variance = bmc_prior_variance_mix( mix, kernel )
%
% Compute the prior variance of BMC when the input distribution is a mixture of
% Gaussians, and the kernel is a scaled Gaussian.
%
% David Duvenaud
% March 2012

K = size(mix.means, 1);

prior_variance = 0;
for k = 1:K
    for j = 1:K
        cur_covariance = kernel.covariance + mix.covs(:, :, k) + mix.covs(:, :, j);
        prior_variance = prior_variance + mix.weights(j) * mix.weights(k) ...
                       * mvnpdf( mix.means(k,:), mix.means(j,:), cur_covariance);
    end
end
prior_variance = prior_variance * kernel.height;
