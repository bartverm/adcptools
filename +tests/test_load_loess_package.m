classdef test_load_loess_package < matlab.unittest.TestCase
    methods(Test)
        function test_loess(testCase)
            % tests adding loess package to path
            testCase.assertTrue(helpers.load_loess_package);
        end
        function test_preload(testCase)
            % tests adding nmea package to path while package is already
            % loaded
            helpers.load_loess_package;
            loess
            testCase.assertTrue(helpers.load_loess_package);
        end
    end
end