test_that("train expression works", {

  a = trainExpression(geneInt = 'test',
                      snps = fread('data/snps.txt'),
                      snpLocs = fread('data/snpLocs.txt'),
                      mediator = fread('data/mediators.txt'),
                      medLocs = read.table('data/medLocs.txt',header=T),
                      covariates = NULL,
                      qtlFull = rbind(fread('data/testCis.txt'),
                                      fread('data/testTra.txt')),
                      numMed = 5,
                      seed = 1,
                      k = 5,
                      fileName = 'test',
                      cisDist = 1e6,
                      parallel = F,
                      prune = F,
                      cores = 1)

  truth = fread('data/truth.txt')

})

