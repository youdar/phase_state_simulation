function [negativeFlag] = testForNegVal(inputMatrix,basePnum,adjPnum)

    % checking if any variables are negative
    % Parameters numbers in inputMatrix
    % 1: klpf       elongation (monomer attachment) for free aggregates
    % 2: klmf       monomer dettachment, for free aggregates 
    % 3: klpm       elongation (monomer attachment) for membrane bounnd aggregates
    % 4: klmm       monomer dettachment, for membrane bounnd aggregates
    % 5: kap        attachemnt rate, of aggregates to te membrane
    % 6: kam        dettachment rate, of aggregates to te membrane
    % 7: kn         free aggregates nucleation rate
    % 8: cnmRatio   cnmRatio = (cnm equivalent)/cn = cnm*p0/cn
    % 9: kbf        free aggregates brakage rate
    % 10: kbm       mambrane aggregates brakage rate
    % 11: y00       initial A-Beta concentration
    % 12: p00       initial Prion concentration
    % 13: rs        Sequestration rate
    % 14: rfm       Out-flow rate
    % 15: n         size of maximum aggregates length
    % 16: m         number of iterations

    negativeFlag = 0;

    for i=1:4
        if inputMatrix(1,basePnum)<0 || inputMatrix(1,adjPnum)<0
            negativeFlag = 1;
        end
    end


end

