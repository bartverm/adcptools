classdef test_readDeployment < matlab.unittest.TestCase
    properties
        data_path = '~/src/adcptools/doc/demo/raw_data/'
    end
    methods(Test)
        function read_testdata(testCase)
            dat = readDeployment('trans', testCase.data_path);
            testCase.assertClass(dat,'struct');
        end
        function wrong_depname(testCase)
            testCase.verifyError(...
                @()readDeployment('wrong', testCase.data_path), ...
                'readDeployment:NoFileFound')
        end
        function wrong_path(testCase)
            testCase.verifyError(...
                @()readDeployment('trans', 'whatever'), ...
                'readDeployment:InexistentPath')
        end
    end
end