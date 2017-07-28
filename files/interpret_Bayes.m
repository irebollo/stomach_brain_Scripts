%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% interpret_Bayes
%%%
%%% This function gives the interpretation of a Bayes Factor, in terms of
%%% how much it supports H1 against H0.
%%% The classification applied here corresponds to the one developped by
%%% Kass & Raftery, Journal of the American Statistical Association, 1995.
%%%
%%% Caution:
%%% The interpretation of the Bayes Factor depends on what is considered as
%%% H0 and what is considered as H1. In the current implementation of the
%%% Bayes Factor computation (in log10), a BF<0 shows evidence in favor of
%%% H1 (there is a true effect between conditions), and BF>0 shows evidence
%%% in favor of H0 (the null hypothesis is true). 
%%%
%%% Input:
%%% - bf_log10: Bayes Factor value in log10
%%%
%%% Output:
%%% - res_Bayes: the result of the Bayes Factor test, i.e. what theory does
%%% our data support and how strongly
%%% - bf: Bayes Factor (not in log10)
%%%
%%% Mariana Babo-Rebelo, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [res_Bayes, bf] = interpret_Bayes(bf_log10)

abs_bf_log10 = abs(bf_log10);
bf = 10^(abs_bf_log10);

if abs_bf_log10 < 0.5
    res_Bayes = 'Anecdotal evidence for ';
elseif abs_bf_log10 >= 0.5 && abs_bf_log10 < 1
    res_Bayes = 'Substantial evidence for ';
elseif abs_bf_log10 >= 1 && abs_bf_log10 < 2
    res_Bayes = 'Strong evidence for ';
elseif abs_bf_log10 >= 2
    res_Bayes = 'Decisive evidence for ';
else
    error('Something wrong with Bayes Factor computation');
end

if bf_log10 < 0
    res_Bayes = [res_Bayes 'H1'];
else
    res_Bayes = [res_Bayes 'H0'];
end

end