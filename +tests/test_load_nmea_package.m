classdef test_load_nmea_package < matlab.unittest.TestCase
    methods(Test)
        function test_load(testCase)
            % tests adding nmea package to path
            testCase.assertTrue(helpers.load_nmea_package);
            msg = nmea.Message;
            testCase.assertClass(msg,'nmea.Message');
        end
        function test_preload(testCase)
            % tests adding nmea package to path while package is already
            % loaded
            helpers.load_nmea_package;
            testCase.assertTrue(helpers.load_nmea_package);
            msg = nmea.Message;
            testCase.assertClass(msg,'nmea.Message');
        end
        
    end
end