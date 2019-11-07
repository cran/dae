#devtools::test("dae")
context("design")

cat("#### Test for designRandomize\n")
test_that("randomize", {
 skip_on_cran()
 library(dae)
 
  #Generate 5 x 5 Latin square
  Treatments <- factor(rep(1:5, times=5))
  RCBD.lay <- designRandomize(allocated = Treatments, 
                              recipient = list(Rows=5, Columns=5), 
                              nested.recipients = list(Columns = "Rows"), 
                              seed = 521814)
  testthat::expect_equal(nrow(RCBD.lay), 25)
  testthat::expect_equal(ncol(RCBD.lay), 3)
  testthat::expect_equal(as.numfac(RCBD.lay$Treatments), 
                         c(1,4,2,3,5,4,2,5,1,3,5,1,3,2,4,4,1,2,5,3,3,4,2,5,1))
  
  Treatments <- factor(designLatinSqrSys(5))
  LSD.lay <- designRandomize(allocated = Treatments,  
                             recipient = list(Rows=5, Columns=5), 
                             seed = 154381)
  testthat::expect_equal(nrow(RCBD.lay), 25)
  testthat::expect_equal(ncol(RCBD.lay), 3)
  testthat::expect_equal(as.numfac(LSD.lay$Treatments),  
                         c(5,4,2,3,1,2,1,4,5,3,3,2,5,1,4,4,3,1,2,5,1,5,3,4,2))
 
 
  #'## Generate a layout for a standard athlete training experiment
  eg1.lay <- designRandomize(allocated = fac.gen(list(Intensities = 3, Surfaces = 3), 
                                                 times = 4),
                             recipient = list(Months = 4, Athletes = 3, Tests = 3), 
                             nested.recipients = list(Athletes = "Months", 
                                                      Tests = c("Months", "Athletes")),
                             seed = 2598)
  testthat::expect_equal(nrow(eg1.lay), 36)
  testthat::expect_equal(ncol(eg1.lay), 5)
  #'## Generate a layout for a simple two-phase athlete training experiment
  #'## Phase 1 - the split-plot design that has already been generated. 
  #'## Phase 2 - randomize tests (and training conditions) to locations, 
  #'##           but Months assigned systematicaly to Batches 
  #'##           so except Batches from the randomization
  eg2.lay <- designRandomize(allocated = eg1.lay,
                             recipient = list(Batches = 4, Locations = 9), 
                             nested.recipients = list(Locations = "Batches"),
                             except = "Batches", 
                             seed = 71230)
  testthat::expect_equal(nrow(eg2.lay), 36)
  testthat::expect_equal(ncol(eg2.lay), 7)
  testthat::expect_equal(as.numfac(eg2.lay$Months), rep(1:4, each=9))
  testthat::expect_equal(eg2.lay$Months,  eg2.lay$Batches)

  
  #A small example to test randomization
  Exp.unit <- list(Squares=2, Rows=3, Columns=3, Halfplots=2, Reps=2)
  Exp.nest <-  list(Columns="Squares", Halfplots=c("Squares","Rows","Columns"),
                    Reps=c("Squares","Rows","Columns","Halfplots"))
  Exp.unit.dat <- fac.gen(Exp.unit)
  Exp.alloc.dat <- data.frame(Trellis = factor(rep(c(1,2,3, 2,3,1, 3,1,2), each=4, times=2)),
                              Method = factor(rep(1:2, each=2, times=18)))

  #randomize supplying list
  Exp.rand.dat <- designRandomize(recipient=Exp.unit, 
                                  nested.recipient=Exp.nest,
                                  allocated=Exp.alloc.dat, 
                                  seed = 646154)
  testthat::expect_equal(nrow(Exp.rand.dat), 72)
  testthat::expect_equal(ncol(Exp.rand.dat), 7)
  Exp.rand.canon <- designAnatomy(list(unit = ~ ((Squares/Columns)*Rows)/Halfplots/Reps,
                                       trt = ~ Trellis*Method),
                                  data = Exp.rand.dat)
  summ.rand <- summary(Exp.rand.canon)
  testthat::expect_equivalent(na.omit(summ.rand$decomp$aefficiency), c(1,1,1))

  #recipient factors in standard order in data frame instead of list
  Exp.std.dat <- designRandomize(recipient=Exp.unit.dat, 
                                 nested.recipients=Exp.nest,
                                 allocated=Exp.alloc.dat, 
                                 seed = 78125, unit.permutation = TRUE)
  testthat::expect_equal(nrow(Exp.std.dat), 72)
  testthat::expect_equal(ncol(Exp.std.dat), 9)
  #Test that all unit factors are the same before and after randomization
  testthat::expect_true(all(unlist(lapply(names(Exp.unit.dat), 
                                          function(facname, lay1, lay2)
                                          {
                                            fac <- lay1[facname]
                                            other.fac <- lay2[facname][1]
                                            all(fac == other.fac)
                                          }, lay1 = Exp.unit.dat, lay2 = Exp.std.dat))))
  Exp.std.canon <- designAnatomy(list(unit = ~ ((Squares/Columns)*Rows)/Halfplots/Reps,
                                      trt = ~ Trellis*Method),
                                 data = Exp.std.dat)
  summ.std <- summary(Exp.std.canon)
  testthat::expect_equivalent(na.omit(summ.std$decomp$aefficiency), c(1,1,1))
  
  #recipient factors in permuted order in a data frame
  Exp.unit.perm.dat <- Exp.unit.dat[Exp.std.dat$.Permutation,]
  Exp.alloc.perm.dat <- Exp.alloc.dat[Exp.std.dat$.Permutation,]
  Exp.perm.dat <- designRandomize(recipient=Exp.unit.perm.dat, 
                                  nested.recipients=Exp.nest,
                                  allocated=Exp.alloc.perm.dat, 
                                  seed = 64614, unit.permutation = TRUE)
  testthat::expect_equal(nrow(Exp.perm.dat), 72)
  testthat::expect_equal(ncol(Exp.perm.dat), 9)
  #Test that all unit factors are the same before and after randomization
  testthat::expect_true(all(unlist(lapply(names(Exp.unit.perm.dat), 
                                          function(facname, lay1, lay2)
                                          {
                                            fac <- lay1[facname]
                                            other.fac <- lay2[facname][1]
                                            all(fac == other.fac)
                                          }, lay1 = Exp.unit.perm.dat, lay2 = Exp.perm.dat))))
  #derandomize the allocated factors and check that have the same as before randomization
  Exp.derand.dat <- Exp.perm.dat[Exp.perm.dat$.Permutation, ]
  testthat::expect_true(all(unlist(lapply(names(Exp.alloc.perm.dat), 
                                          function(facname, lay1, lay2)
                                          {
                                            fac <- lay1[facname]
                                            other.fac <- lay2[facname][1]
                                            all(fac == other.fac)
                                          }, lay1 = Exp.alloc.perm.dat, lay2 = Exp.derand.dat))))
  testthat::expect_equal(as.numfac(Exp.perm.dat$Trellis)[1:12], rep(c(1,3,2), each=4))
  Exp.perm.canon <- designAnatomy(list(unit = ~ ((Squares/Columns)*Rows)/Halfplots/Reps,
                                      trt = ~ Trellis*Method),
                                 data = Exp.perm.dat)
  summ.perm <- summary(Exp.perm.canon)
  testthat::expect_equivalent(na.omit(summ.perm$decomp$aefficiency), c(1,1,1))

  
  #Test factors when not in order of columns does not match nesting 
  RCBD.sys <- cbind(fac.gen(list(rep = 2, plot=1:3, block = c("I","II"))),
                     tr = factor(rep(1:3, each=2, times=2)))
  ## obtain randomized layout, speciying 
  RCBD.lay <- designRandomize(allocated = RCBD.sys["tr"], 
                              recipient = RCBD.sys[c("rep", "block", "plot")], 
                              nested.recipients = list(plot = c("block","rep"), 
                                                       block="rep"), 
                              seed = 9719532, 
                              unit.permutation = TRUE)
  RCBD.canon <- designAnatomy(list(unit = ~ rep/block/plot, trt = ~ tr),
                              data = RCBD.lay)
  summ.RCBD <- summary(RCBD.canon)
  testthat::expect_equal(summ.RCBD$decomp$aefficiency[3], 1)
  #Test that the order of unit factors are the same in RCBD.sys and RCBD.lay
  testthat::expect_true(all(unlist(lapply(c("rep", "block", "plot"), 
                                          function(facname, lay1, lay2)
                                          {
                                            fac <- lay1[facname]
                                            other.fac <- lay2[facname][1]
                                            all(fac == other.fac)
                                          }, lay1 = RCBD.sys, lay2 = RCBD.lay))))
  #Test derandomized treatments the same as in RCBD.sys
  testthat::expect_true(all(RCBD.lay[RCBD.lay$.Permutation, "tr"] == RCBD.sys$tr))
  
  # Test with recipient columns listed in the same order as nesting, tr a factor 
  tr = factor(rep(1:3, each=2, times=2))
  RCBD.unit <- list(rep = 2, plot=c(0,2,4), block = c("I","II"))
  RCBD.unit <- fac.gen(RCBD.unit)
  RCBD.nest <- list(plot = c("block","rep"), block="rep")
  RCBD.lay <- designRandomize(recipient=RCBD.unit, nested.recipients=RCBD.nest, 
                             allocated=tr, seed=7197132)
  RCBD.canon <- designAnatomy(list(unit = ~ rep/block/plot, trt = ~ tr),
                              data = RCBD.lay)
  summ.RCBD <- summary(RCBD.canon)
  testthat::expect_equal(summ.RCBD$decomp$aefficiency[3], 1)
  #Test that the order of unit factors are the same in RCBD.sys and RCBD.lay
  testthat::expect_true(all(unlist(lapply(c("rep", "block", "plot"), 
                                          function(facname, lay1, lay2)
                                          {
                                            fac <- lay1[facname]
                                            other.fac <- lay2[facname][1]
                                            all(fac == other.fac)
                                          }, lay1 = RCBD.unit, lay2 = RCBD.lay))))
  #Test derandomized treatments the same as in RCBD.sys
  testthat::expect_true(all(RCBD.lay[RCBD.lay$.Permutation, "tr"] == tr))
  
  #Test except
  LS.std.unit <- list(row = c("I","II","III","IV"), col = 4)
  treat <- factor(designLatinSqrSys(4))
  LS.std.lay <- designRandomize(recipient=LS.std.unit, 
                                allocated=treat, 
                                seed=7197132, unit.permutation = TRUE) 
  testthat::expect_equal(LS.std.lay$.Permutation, 
                         c(1,3,4,2,13,15,16,14,9,11,12,10,5,7,8,6))
  
  #check except for non-nested factors
  LS.noran <- designRandomize(recipient=LS.std.unit, 
                              allocated=treat, 
                              except=c("row","col"), 
                              seed=7197132, unit.permutation = TRUE) 
  testthat::expect_equal(LS.noran$.Permutation, 1:16)
  
  #Complicated nesting
  recip <- list(S = 2, r = 2, c = 3, R = 2, C = 2)
  nest <- list(r = c("R","C","S"), c = c("S","C"), R = "S", C = "S")
  alloc <- data.frame(tr = factor(rep(1:4, times=12)))
  r <- designRandomize(recipient=recip, 
                       nested.recipients=nest, 
                       allocated=alloc, seed=7197132)
  r.canon <- designAnatomy(list(unit = ~ S/(R*(C/c)) + S:R:C:r, trt = ~ tr),
                           data = r)
  summ.r <- summary(r.canon)
  testthat::expect_equal(summ.r$decomp$aefficiency[c(2,4,7)], c(1,1,1))
  
  #Some examples to test except
  l <- designRandomize(recipient=recip, 
                       nested.recipients=nest, 
                       except=c("R","C"), 
                       allocated=alloc, seed=7197132)
  testthat::expect_equal(as.numfac(l$tr), rep(1:4, times = 12))

  #RCBD example
  RCBD.unit <- list(rep = 2, plot=c(0,2,4), block = c("I","II"))
  RCBD.unit <- fac.gen(RCBD.unit)
  unrand.rows <-as.numeric(rownames(with(RCBD.unit, RCBD.unit[order(rep,block,plot),])))
  RCBD.nest <- list(plot = c("block","rep"), block="rep")
  RCBD.alloc <- data.frame(tr = factor(rep(1:3, each=2, times=2)))
  RCBD.lay <- designRandomize(recipient=RCBD.unit, 
                              nested.recipients=RCBD.nest,
                              allocated=RCBD.alloc, 
                              seed=7197132, unit.permutation = TRUE)
  RCBD.lay <- with(RCBD.lay, RCBD.lay[order(rep,block,plot),])
  testthat::expect_equal(RCBD.lay$.Permutation, c(6,4,2,1,3,5,7,11,9,10,8,12))
  
  RCBD.noblk.ran <- designRandomize(recipient=RCBD.unit, nested.recipients=RCBD.nest, 
                                    except="block", 
                                    allocated=RCBD.alloc, 
                                    seed=7197132, unit.permutation = TRUE)
  RCBD.noblk.ran <- with(RCBD.noblk.ran, RCBD.noblk.ran[order(rep,block,plot),])
  #Block I should have 1,3,5,7,9,11; II should have 2,4,6,8,10,12
  testthat::expect_true(all(c(1,3,5,7,9,11) %in% 
                              RCBD.noblk.ran$.Permutation[RCBD.noblk.ran$block == "I"]))
  testthat::expect_true(all(c(2,4,6,8,10,12) %in% 
                              RCBD.noblk.ran$.Permutation[RCBD.noblk.ran$block == "II"]))

  RCBD.noplt.ran <- designRandomize(recipient=RCBD.unit, nested.recipients=RCBD.nest, 
                                    except="plot", 
                                    allocated=RCBD.alloc, 
                                    seed=7197132, unit.permutation = TRUE)
  RCBD.noplt.ran <- with(RCBD.noplt.ran, RCBD.noplt.ran[order(rep,block,plot),])
  testthat::expect_equal(as.numfac(RCBD.noplt.ran$tr), rep(1:3, times = 4))
  
  RCBD.plt.ran <- designRandomize(recipient=RCBD.unit, nested.recipients=RCBD.nest, 
                                  except=c("rep","block"), 
                                  allocated=RCBD.alloc, 
                                  seed=7197132, unit.permutation = TRUE)
  RCBD.plt.ran <- with(RCBD.plt.ran, RCBD.plt.ran[order(rep,block,plot),])
  #triples should contain the same nos as in the rownames
  testthat::expect_equal(as.numfac(RCBD.plt.ran$.Permutation), 
                         c(3,1,5,2,6,4,9,7,11,8,10,12))
  
})
