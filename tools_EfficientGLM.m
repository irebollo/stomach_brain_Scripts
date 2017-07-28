function res = tools_EfficientGLM(nn_data, nn_desmat)
%EFFICIENTGLM computes a multiple linear regression in parallel along
%several dimensions.
%   - "data": matrix of data (rows: replication, columns: data to be 
%       processed in parallel).
%   - "desmat": design matrix (rows: replication, column: regressors).
%
% Script history:
% * Florent Meyniel: initial version.
% * Maxime Maheu: removing NaN checking to improve efficiency.

% Use the rank-revealing QR to remove dependent columns of X.
ndata        = size(nn_data, 2);
[n, nreg]    = size(nn_desmat);
[Q, R, perm] = qr(nn_desmat, 0);
p = sum(abs(diag(R)) > max(n, nreg)*eps(R(1)));
if p < nreg
    R = R(1:p, 1:p);
    Q = Q(:, 1:p);
    perm = perm(1:p);
end

% Compute the LS coefficients, filling in zeros in elements corresponding
% to rows of X that were thrown out.
b = zeros(nreg, ndata);
b(perm, :) = R \ (Q'*nn_data);

res = b;

end