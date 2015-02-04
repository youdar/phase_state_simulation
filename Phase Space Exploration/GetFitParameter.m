function FitVal = GetFitParameter(gof)

% Evaluate gof and make sure we don't get NaN as a value

if isnan(gof.adjrsquare)
    FitVal = 1;
else
    FitVal = gof.adjrsquare;
end