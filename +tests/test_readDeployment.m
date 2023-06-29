classdef test_readDeployment < matlab.unittest.TestCase
% Test for readDeployment function
    properties
        data_path = ''
        dep_name = 'trans'
    end
    methods
        function obj = test_readDeployment(varargin)
            obj.data_path = fullfile(helpers.adcptools_root,...
                'doc','sample_data',...
                'rdi_muara_muntai_bend');
        end
    end
    methods(Test)
        function read_testdata(testCase)
            dat = rdi.readDeployment(testCase.dep_name,...
                testCase.data_path);
            testCase.assertClass(dat,'struct');
        end
        function wrong_depname(testCase)
            testCase.verifyError(...
                @()rdi.readDeployment('wrong',...
                testCase.data_path), ...
                'readDeployment:NoFileFound')
        end
        function wrong_path(testCase)
            testCase.verifyError(...
                @()rdi.readDeployment('trans', 'whatever'), ...
                'readDeployment:InexistentPath')
        end
    end
end