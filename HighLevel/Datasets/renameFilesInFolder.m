function [ ] = renameFilesInFolder( path, ext )
    
    filespec = sprintf('%s%s%s%s', './', path, '/*.', ext);
    files = dir(filespec);
    
    i = 1;
    filebasepath = sprintf('%s%s%s', './', path, '/');
    for file = files'
        filepathold = sprintf('%s%s', filebasepath, file.name);
        filepathnew = sprintf('%s%03d%s%s', filebasepath, i, '.', ext);
        movefile(filepathold, filepathnew);
        i = i + 1;
    end
end