test_that("train expression works", {

  dir = system.file('extdata', package='MOSTWAS')
  
  a = trainExpression(geneInt = 'test',
                      snps = data.table::fread(file.path(dir,'snps.txt')),
                      snpLocs = data.table::fread(file.path(dir,'snpLocs.txt')),
                      mediator = data.table::fread(file.path(dir,'mediators.txt')),
                      medLocs = read.table(file.path(dir,'medLocs.txt'),header=T),
                      covariates = NULL,
                      qtlFull = rbind(data.table::fread(file.path(dir,'testCis.txt')),
                                      data.table::fread(file.path(dir,'testTra.txt'))),
                      numMed = 5,
                      seed = 1,
                      k = 5,
                      cisDist = 1e6,
                      parallel = F,
                      prune = F,
                      cores = 1)

  truth = data.table::fread(file.path(dir,'truth.txt'))

})

