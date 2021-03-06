function [TestVal,AccomulationFlag] = GetTestValues(TestVal,TestFlag,yLast,gLast,alarmFlag)

% Get test values for ExplorPhaseSpace program
% TestVal is an array (1,4). it idicates wether a test passed or not. pass test 0, fail test 1
% AccomulationFlag check if we have accomulation at the end of the y or x arrays
% If we have aggregates at the and we need to increase the maximum size of y and x arrays 

AccomulationFlag = 0;

for j=1:4
    TestVal(j) = TestFlag{j};             % test fail if we have aggregates for at leat FailCounter steps
    AccomulationFlag = AccomulationFlag + (yLast{j}>0 || gLast{j}>0);
    disp(['Lab number: ' num2str(j)])
    disp('==========================')
    disp(['Alarm Value (Particle number factor): ' num2str(alarmFlag{j})])
    disp(['Test Value: ' num2str(TestFlag{j})])
    disp(['yLast: ' num2str(yLast{j})])
    disp(['gLast: ' num2str(gLast{j})])
    disp('--------------------------')
end
disp(['AccomulationFlag: ' num2str(max(AccomulationFlag))])
disp('==========================')