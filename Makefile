all: roxygen build

roxygen: CoxHD
	R -e 'library(roxygen2); roxygenize("CoxHD")'
 
check: CoxHD
	R CMD check CoxHD

bump: CoxHD
	git describe --abbrev=0 | xargs -I%  git rev-list %..HEAD --count | awk '{print $$0 +1}' | xargs -I% sed -i'~' '/^Version/ s/\([0-9]*\.[0-9]*\)\(.[0-9]*\)/\1.%/g' CoxHD/DESCRIPTION
	date | xargs -I% sed -i'~' 's/Date:.*/Date: %/g' CoxHD/DESCRIPTION

build: CoxHD clean
	cp -rf CoxHD CoxHD.build
	R CMD build --no-build-vignettes CoxHD.build

install: CoxHD build
	R CMD INSTALL CoxHD.build
	
clean: 
	rm -rf CoxHD.build CoxHD.Rcheck
