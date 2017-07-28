%% Bayesian t-test 
%
%  xdat corresponds to the data to be tested against the null hypothesis.
%  If you want to test if condition A differs from condition B, compute the
%  following:
%  >> xdat = <mean condition A for each subject> - <mean condition B for
%  each subject>;
%
%  xref corresponds to the reference effect size (a constant). If it comes
%  from a paired difference vector, compute the following:
%  >> xref = mean(x)/std(x)
%  xref may correspond to the smallest effect size you could find for a
%  p-value of 0.05. In that case of nobs = 16, for instance, enter:
%  >> xref = 0.533;
%
%  This function returns the Bayes factor in log-10 units.
%  See Kass & Raftery 1995 for the interpretation of the value of the Bayes
%  Factor.
%
%  Works with Matlab2012 (not 2009, 2010, 2011).
%
%  Valentin Wyart <valentin.wyart@ens.fr>


function bf = my_ttest_bayes(xdat, xref)

% standardize data
xdat = xdat/std(xdat);

nobs = length(xdat);
sref = 1/sqrt(nobs); % reference effect size s.d. (do not change!)

% define prior function
prior_type = 'two-sided';
switch prior_type
    case 'one-sided'
        fprior = @(m)normpdf(m,xref,sref);
    case 'two-sided'
        fprior = @(m)0.5*(normpdf(m,+xref,sref)+normpdf(m,-xref,sref));
    otherwise
        error('undefined prior type!');
end

%%%%%%%%%%%%%%%%%%%% test data
% compute log(model evidence) for hypothesis H0 (no effect)
logev_h0 = sum(log(normpdf(xdat,0)));

% compute log(model evidence) for hypothesis H1 (true effect)
logev_h1 = log(integral(@(m)arrayfun(@(mi)exp(sum(log(normpdf(xdat,mi)))).*fprior(mi),m),-inf,+inf));

% compute Bayes factor (in log-10 units)
bf = (logev_h1-logev_h0)/log(10);

