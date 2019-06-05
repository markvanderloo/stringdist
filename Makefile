
doc: 
	R -s -e "pkgload::load_all('pkg');roxygen2::roxygenize('pkg')"

pkg: doc
	R CMD build pkg

check: doc
	R CMD build pkg
	R CMD check *.tar.gz

cran: doc
	R CMD build pkg
	R CMD check --as-cran *.tar.gz

test: doc
	R -s -e "tinytest::build_install_test('pkg')"

manual: doc
	R CMD Rd2pdf --force -o manual.pdf ./pkg

revdep: pkg
	rm -rf revcheck
	mkdir revcheck
	mv *.tar.gz revcheck
	R -s -e "out <- tools::check_packages_in_dir('revcheck',reverse=list(which='most'),Ncpus=3); print(summary(out)); saveRDS(out, file='revcheck/output.RDS')"


