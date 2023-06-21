classdef RemoveFileFixture < matlab.unittest.fixtures.Fixture
    properties(Access = private)
        tmp_filename(1,1) string
        filename(1,1) string
    end
    methods
        function obj=RemoveFileFixture(fname)
            assert(exist(fname,"file"), 'File does not exist.')
            obj.filename = fname; 
        end
        function setup(fixture)
            fixture.tmp_filename = tempname;
            movefile(fixture.filename,fixture.tmp_filename);
            fixture.SetupDescription = "Moved file " +...
                fixture.filename + " to " + fixture.tmp_filename;
        end
        function teardown(fixture)
            movefile(fixture.tmp_filename,fixture.filename);
            fixture.TeardownDescription = "Moved file " +...
                fixture.tmp_filename + " back to " + fixture.filename;
        end
    end
    methods (Access=protected)
        function tf = isCompatible(fixture1,fixture2)
            tf = strcmp(fixture1.filename,fixture2.filename);
        end
    end
end