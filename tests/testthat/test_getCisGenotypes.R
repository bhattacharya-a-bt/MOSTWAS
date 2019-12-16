test_that("get cis-genotypes works", {

  set.seed(1)
  biomInt = 'test'
  locs = data.frame(medID = c('test','test2','test3'),
                    chr = sample(1:22,3),
                    pos = sample(1:1e10,3))

  nS.in = 100
  nS.out = 150
  snpLocs = data.frame(snpid = paste0('snp',1:(nS.in + nS.out)),
                       chr = rep(locs$chr[1],(nS.in + nS.out)),
                       pos = c(runif(nS.in,locs$pos[1]-1e6,locs$pos[1]+1e6),
                               runif(nS.out,locs$pos[1]-2e6,locs$pos[1]-1e6)))

  p = 200
  snpsShell = matrix(rnorm(200 * (nS.out + nS.in)),ncol = 200)
  snps = as.data.frame(cbind(sample(paste0('snp',1:(nS.in + nS.out))),
                             snpsShell))
  colnames(snps) = c('SNP',paste0('sample',1:p))
  snps$SNP = as.character(snps$SNP)

  cisGeno = getCisGenotypes(biomInt = 'test',
                            locs = locs,
                            snps = snps,
                            snpLocs = snpLocs,
                            cisDist = 1e6)

  expect_equal(nrow(cisGeno$snpCur),nS.in)
  expect_equal(length(cisGeno$snpList),nS.in)
  expect_equal(nrow(cisGeno$thisSNP),nS.in)


})
