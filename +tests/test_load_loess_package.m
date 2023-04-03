classdef test_load_loess_package < matlab.unittest.TestCase
    methods(Test)
        function test_loess(testCase)
            % tests adding loess package to path
            testCase.assertTrue(helpers.load_loess_package);
        end
        function test_compile(testCase)
            % tests compiling loess
            testCase.applyFixture(tests.RemoveFileFixture(...
                fullfile('loess_submodule',['loess.',mexext])));
            testCase.assertTrue(helpers.load_loess_package);
        end
        function test_preload(testCase)
            % tests adding nmea package to path while package is already
            % loaded
            helpers.load_loess_package;
            testCase.assertTrue(helpers.load_loess_package);
        end
    end
end