% script that calls several regression tests for comparing
% multiscale and finescale solution by looking at relative error norms. 
val = 1;
clear;
close all;

load prevRes
load prevCases
compare = true;
testCases = what(fullfile(SAMSIMROOTDIR, 'regression', 'tests'));
testCases = testCases.m;
htmldir   = fullfile(SAMSIMROOTDIR, 'regression', 'html');

numCases = length(testCases);

threshold = 0.1;
%% run test-cases
disp('**************************************************');
disp('  Running analysis on multiscale simulation code ')
disp('**************************************************');
%%{


msg_error   = 'SYNTAX ERROR';
msg_ok      = 'OK';
msg_warning = 'CHECK THIS';

T = cell(numCases, 4);

for i = 1:numCases
    %%  Test number i_
    test_name = testCases{i}(1:end-2);
    T{i,1} = sprintf('<a href=%s.html> %s </a>\n', ...
                     test_name, test_name);
    str = [' Running ',testCases{i}, '...'];
    str = [str, repmat(' ', [1, 35-numel(str)])];
    fprintf(str);
    try
        fun    = str2func(test_name);
        T{i,2} = fun(false);
        publish(testCases{i}(1:end-2), struct('showCode', false, 'outputDir', htmldir));
        if T{i,2} < threshold 
            T{i,4} = msg_ok;
            fprintf('%s\n',msg_ok);
        else
            T{i,4} = msg_warning;
            fprintf('%s\n',msg_warning);
        end
    catch ME
        T{i,4} = msg_error;
        fprintf('%s\n',msg_error);
        T{i,1} = sprintf('<a> %s </a>\n', testCases{i});
    end
end
%}
%% Plot and display results
disp('* Results *');


head  = {'Testcase','Current result', 'Stored result', 'Status'};
color = cell(numCases, 1);

[color{strcmp(msg_ok,      {T{:,4}})}] = deal('#66FF00');
[color{strcmp(msg_warning, {T{:,4}})}] = deal('#FF9900');
[color{strcmp(msg_error,   {T{:,4}})}] = deal('#FF0000');





%
%  Save checksums instead, to detect when a test has been changed.
%  This should trigger a warning message.
%
if ~isequal(prevCases, testCases)
   disp(['Tests not compatible for comparison ', ...
         'with previous test-run.'])
   compare = false;   
   [T{:,3}] = deal('NA');
end
writeSummary(head, T, color, [htmldir, filesep, 'summary.html'])

normRes = norm(res);
%{
if compare %compare result with result from a previous run
   diff = prevRes'-res;
   figure;
   subplot(2,1,1)
      plot(1:length(res), res); hold on;
      plot(1:length(res), prevRes, 'r*');
      title('Plot of result against previous result');
      ylabel('Relative error norm')
      xlabel('Test case')
      set(gca, 'xtick', 1:length(res));
      legend('result', 'previous result', 'Location', 'Best');
   subplot(2,1,2)
      plot(1:length(res), diff); 
      title('Plot of difference');
      set(gca, 'xtick', 1:length(res))
      ylabel('Difference in relative error')
      xlabel('Test case')
   
   disp(['Norm of result   : ', num2str(normRes)])
   disp(['norm(prevRes-res): ', num2str(norm(diff))]);
   
   if norm(diff)~=0
      disp('Check your results!!');
      disp(['Consider saving ''res'' as ''prevRes''', ...
            ' if results are improved:']);
      disp('(prevRes = res; save regression/prevRes;)'); 
   end

else %when results can not be compared with those of previos run
   figure;
   plot(1:length(res), res);
   set(gca, 'xtick', 1:length(res))
   title('Norm of testing');
   ylabel('Relative error norm')
   xlabel('Test case')
   
   disp(['Norm of result: ', num2str(normRes)])
   disp(['NB: Consider saving ''res'' and ''testCases''', ...
         ' if test cases have been changed:']);
   disp('(prevRes   = res; save regression/prevRes;)');
   disp('(prevCases = testCases; save regression/prevCases;)');
end


%}


