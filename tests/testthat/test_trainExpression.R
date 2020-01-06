test_that("train expression works", {

  a = trainExpression(geneInt = 'test',
                      snps = data.table::fread('sample_data/snps.txt'),
                      snpLocs = data.table::fread('sample_data/snpLocs.txt'),
                      mediator = data.table::fread('sample_data/mediators.txt'),
                      medLocs = read.table('sample_data/medLocs.txt',header=T),
                      covariates = NULL,
                      qtlFull = rbind(data.table::fread('sample_data/testCis.txt'),
                                      data.table::fread('sample_data/testTra.txt')),
                      numMed = 5,
                      seed = 1,
                      k = 5,
                      fileName = 'test',
                      cisDist = 1e6,
                      parallel = F,
                      prune = F,
                      cores = 1)

  truth = data.table::fread('sample_data/truth.txt')

})

