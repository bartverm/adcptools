classdef test_documentation < matlab.unittest.TestCase
% Tests the generation of the documentation
    methods(Test)
        function test_generate_documentation(testCase)
            import matlab.unittest.fixtures.CurrentFolderFixture
            testCase.applyFixture(CurrentFolderFixture("doc"))
            publish_all_doc
        end
    end
end