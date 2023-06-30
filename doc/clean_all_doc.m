function clean_all_doc()
    rmdir(fullfile('../html', 'helpsearch*'))
    delete(fullfile('../html','*.html'))
end