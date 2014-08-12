all: roxygen build

roxygen: CoxHD
	R -e 'library(roxygen2); roxygenize("CoxHD")'
 
check: CoxHD
	R CMD check CoxHD

build: CoxHD
	R CMD build --no-build-vignettes CoxHD

install: CoxHD
	R CMD INSTALL CoxHD
