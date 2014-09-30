all: roxygen build

roxygen: CoxHD
	R -e 'library(roxygen2); roxygenize("CoxHD")'
 
check: CoxHD
	R CMD check CoxHD

build: CoxHD
	cp -rf CoxHD CoxHD.build
	git describe --abbrev=0 | xargs -I%  git rev-list %..HEAD --count | xargs -I% sed '/^Version/ s/$$/\.%/g' CoxHD/DESCRIPTION > CoxHD.build/DESCRIPTION
	git log --pretty=format:'%ai' -n 1 HEAD | xargs -I% sed -i'~' 's/Date:.*/Date: %/g' CoxHD.build/DESCRIPTION
	R CMD build --no-build-vignettes CoxHD.build

install: CoxHD build
	R CMD INSTALL CoxHD.build
	
clean: 
	rm -rf CoxHD.build CoxHD.Rcheck
