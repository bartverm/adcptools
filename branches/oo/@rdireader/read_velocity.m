function read_velocity(obj)
for cntens=1:nens                                                  % Loop for all ensembles
    CurrentNbins=dataout.nbins(dataout.FileNumber(cntens));        % Read the number of bins in current file
    Ndatablock=find(DataHeader{cntens}(:,1)==256, 1);              % Search for velocity data in current ensemble
    if isempty(Ndatablock)                                         % If data is not found
        continue                                                   % Leave matrix empty (actually full of zeros, change to -32768???)
    else                                                           % Otherwise
        fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);      % Calculate position of velocity data in file
        dataout.VEL(1:CurrentNbins,cntens,:)=...
            reshape(readVEL(fileid(cntens),...
            fpos,CurrentNbins),[],1,4);                            % Read velocity data of current ensemble
    end                                                            % Note: for explanation on velocity data see under the reading function
end
