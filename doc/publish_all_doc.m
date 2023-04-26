function publish_all_doc()
%     recursive_publish('./', '../html')
    hdir = fullfile(helpers.adcptools_root, 'html');
    publish('main.m', outputDir=hdir);
    publish('reading_data.m', outputDir=hdir);
    publish('initial_inspection.m', outputDir=hdir);
    publish('repeat_transect_processing.m', outputDir=hdir);
    publish('transect_bathymetry.m', outputDir=hdir);
    publish('transect_cross_section.m', outputDir=hdir);
    publish('transect_data_selection.m', outputDir=hdir);
    publish('transect_mesh_construction.m', outputDir=hdir);
    publish('transect_post_processing.m', outputDir=hdir);
    publish('transect_solving_data.m', outputDir=hdir);

    builddocsearchdb(hdir)

    web(fullfile(hdir, 'main.html'))

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

end