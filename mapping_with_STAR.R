fqsplit <- readRDS('~/projects/OB_RNAseq/resources/fqsplit_list.rds')
mydir <- '/earth/OB_RNAseq/'
system(paste('/opt/STAR-2.7.10b/bin/Linux_x86_64/STAR --genomeLoad LoadAndExit --genomeDir /earth/public_data/STAR_index/grcm39'))
mapply(function(x, y) {
    system(paste0(
        '/opt/STAR-2.7.8a/bin/Linux_x86_64/STAR --runThreadN 70',
        ' --genomeDir /earth/public_data/STAR_index/grcm39_APP-NLGF',
        ' --genomeLoad LoadAndKeep --readFilesCommand zcat',
        ' --readFilesIn ', x[1], ' ', x[2],
        ' --outFileNamePrefix ', mydir, 'aligned/', y, '_',
        ' --quantMode GeneCounts',
        ' --outSAMtype BAM SortedByCoordinate',
        ' --outSAMattributes All',
        ' --outReadsUnmapped Fastx',
        ' --outBAMsortingThreadN 8',
        ' --chimSegmentMin 30',
        ' --chimOutType WithinBAM',
        ' --peOverlapNbasesMin 25',
        ' --peOverlapMMp 0.01',
        ' --bamRemoveDuplicatesType UniqueIdentical',
        ' --limitBAMsortRAM 100000000000'
    ))
}, x = fqsplit, y = names(fqsplit))
system(paste('/opt/STAR-2.7.10b/bin/Linux_x86_64/STAR --genomeLoad Remove --genomeDir /earth/public_data/STAR_index/grcm39'))