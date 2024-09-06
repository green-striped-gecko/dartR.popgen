library(dartR.popgen)
library(dartR.base)
possums.gl <- dartR.data::possums.gl
platypus.gl <- dartR.data::platypus.gl
testset.gl <- dartR.data::testset.gl
testset.gs <- dartR.data::testset.gs
LBP <- dartR.data::LBP
EYR <- dartR.data::EYR
gl.run.popcluster(x = testset.gs, popcluster.path = "/home/lau/popcluster_bin", 
                  output.path = "/home/lau/chingching_dev_Lnx_dartR.popgen/testing_res",
                  filename = "test_lnx", minK = 1, maxK = 4, rep = 1, PopData = 1, location = 0, loc_admixture = 1) 

Q <- gl.plot.popcluster(testset.gs, filename="test_lnx", 
                   input.dir="/home/lau/chingching_dev_Lnx_dartR.popgen/testing_res", plot.K=4,
                   color_clusters=NULL)


Sys.getenv("MALLOC_CHECK_")

print(Sys.setenv(MALLOC_CHECK_ = 2))

Sys.unsetenv("MALLOC_CHECK_")

old.path <- "~/chingching_dev_Lnx_dartR.popgen/"