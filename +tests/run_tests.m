import matlab.unittest.TestRunner;
import matlab.unittest.plugins.CodeCoveragePlugin;
import matlab.unittest.plugins.TestReportPlugin;
import matlab.unittest.plugins.XMLPlugin;
import matlab.unittest.plugins.codecoverage.CoverageReport;

suite = testsuite('+tests', 'IncludeSubfolders', false);

[~,~] = mkdir('code-coverage');
[~,~] = mkdir('test-results');

runner = TestRunner.withTextOutput();
runner.addPlugin(TestReportPlugin.producingPDF('test-results/results.pdf'));
runner.addPlugin(CodeCoveragePlugin.forFolder({'.'}, 'IncludingSubfolders', true,...
  "Producing",CoverageReport('code-coverage')));

results = runner.run(suite);
display(results);

assertSuccess(results);
