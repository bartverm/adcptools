function publish_all_doc(varargin)
% publishes all adcptools documentation
%
%   publish_all_doc() publishes all ADCPtools documentation in the html
%   folder. Only publishes files that were modified after html creation.
%   Files that have a corresponding html that is newer, are run, but not
%   published
%
%   publish_all_doc(["name1 name2"]) Do the above, but only for the given
%   files
%
%   publish_all_doc(..., "force_publish",...) force publication of all
%   files, regardless of what was generated before.

    close all
    force_publish = false;
    names = [ ...
        "main",...
        "reading_data",...
        "initial_inspection",...
        "vertical_positioning",...
        "horizontal_positioning",...
        "repeat_transect_processing",...
        "transect_data_selection",...
        "transect_bathymetry",...
        "transect_cross_section",...
        "transect_mesh_construction",...
        "transect_post_processing",...
        "transect_solving_data"...
    ];

    for carg = 1:nargin
        curarg = varargin{carg};
        if isstring(curarg)
            if strcmp(curarg,"force_publish")
                force_publish = true;
            else
                assert(all(ismember(curarg, names),"all"),...
                    "All given names should match a filename")
                names = curarg;
            end
        else
            error("Can only handle string input")
        end
    end


    % note that the order of publishing the files matters, since the code
    % in the help files builds on code from previous files
    is_modified = publish_or_run(force_publish, names);
    close all
    
    if is_modified
        hdir = fullfile(helpers.adcptools_root, 'html');
        builddocsearchdb(hdir)
    end
    
    open_adcptools_documentation
end

function mod = publish_or_run(fpublish, name)
    mdir = fullfile(helpers.adcptools_root, 'doc');
    hdir = fullfile(helpers.adcptools_root, 'html');
    mod = false;
    n_files = numel(name);
    for cf = 1:n_files
        generate_doc = false;
        % check if html exist
        html_name = fullfile(hdir, name(cf) + ".html");
        m_name = fullfile(mdir, name(cf) + ".m");
        assert(exist(m_name,'file'), ...
            m_name + " does not exist")
        if exist(html_name,'file')
            d = dir(html_name);
            html_time = d.datenum;
            d = dir(m_name);
            m_time = d.datenum;
            if html_time < m_time
                generate_doc = true;
            end
        else
            generate_doc = true;
        end
        if fpublish || generate_doc
            publish(m_name, outputDir = hdir);
            mod = true;
        else
            evalin('base', "run('" + m_name + "')")
        end
    end
end

function recursive_publish(path, outpath)
    files = dir(fullfile(path, '*.m'));
    for cf = 1:numel(files)
        publish(files(cf).name, outputDir=outpath,...
            evalCode=false,showCode=false);
    end
    subpkg = dir(fullfile(path, '+*'));
    subpkg(~[subpkg.isdir])=[];
    for cs = 1:numel(subpkg)
        html_dir = fullfile(outpath, subpkg(cs).name(2:end));
        if ~exist(html_dir,"dir")
            mkdir(html_dir)
        end
        recursive_publish(subpkg(cs).name, html_dir)
    end
end

