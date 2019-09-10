# Performa  file compression with gzip command
SC.Gzip<-function(file.name, remove=T) {
    # file.name; the name of the file to be compressed
    # remove; whether keep the original file
    if (remove) system(paste('gzip -f ', file.name))
    else system(paste('gzip -c -f ', file.name, ' > ', file.name, '.gz', sep='' ));
}