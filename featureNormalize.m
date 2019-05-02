function [X_norm, mu, sigma] = featureNormalize(X)
%FEATURENORMALIZE Normalizes the features in X 
%   FEATURENORMALIZE(X) returns a normalized version of X where
%   the mean value of each feature is 0 and the standard deviation
%   is 1. This is often a good preprocessing step to do when
%   working with learning algorithms.

% You need to set these values correctly
X_norm = X;
mu = zeros(1, size(X, 2));
sigma = zeros(1, size(X, 2));

% ====================== YOUR CODE HERE ======================
% Instructions: First, for each feature dimension, compute the mean
%               of the feature and subtract it from the dataset,
%               storing the mean value in mu. Next, compute the 
%               standard deviation of each feature and divide
%               each feature by it's standard deviation, storing
%               the standard deviation in sigma. 
%
%               Note that X is a matrix where each column is a 
%               feature and each row is an example. You need 
%               to perform the normalization separately for 
%               each feature. 
%
% Hint: You might find the 'mean' and 'std' functions useful.
%       

% for i = 1:2   %for i in 1 to the number of columns (features) of X
%     feat_mean = mean(X(:,i));           %mean of feature i
%     mu(i) = sum(X(:,i)) - feat_mean;    %mu is the feature i - the mean of feature i
%     feat_std = std(X(:,i));             % standard deviation of feature i 
%     sigma(i) = sum(X(:,i)) / feat_std;  %sigma is the feature i / std of i
% end

mu = mean(X); %row vector (1xn)
sigma = std(X); %row vector (1xn)
m = size(X,1);  %number of samples (rows)
mu_matrix = ones(m,1) * mu; %mxn matrix with the mean in every cell
sigma_matrix = ones(m,1) * sigma;%mxn matrix with sig in every cell
%the above is because mu is 1xn, & mx1 * 1xn = mxn matrix
X_norm = (X - mu_matrix) ./ sigma_matrix;   
% normalised X matrix is X minus mu elementwise, divided element wise by
% sigma matrix
% ============================================================

end
