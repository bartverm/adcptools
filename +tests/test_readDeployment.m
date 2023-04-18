classdef test_readDeployment < matlab.unittest.TestCase
% Test for readDeployment function
    properties
        data_path = 'doc/demo/raw_data/'
    end
    methods(Test)
        function read_testdata(testCase)
            dat = rdi.readDeployment('trans', testCase.data_path);
            testCase.assertClass(dat,'struct');
        end
        function wrong_depname(testCase)
            testCase.verifyError(...
                @()rdi.readDeployment('wrong', testCase.data_path), ...
                'readDeployment:NoFileFound')
        end
        function wrong_path(testCase)
            testCase.verifyError(...
                @()rdi.readDeployment('trans', 'whatever'), ...
                'readDeployment:InexistentPath')
        end
    end
end