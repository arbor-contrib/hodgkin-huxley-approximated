.DEFAULT_GOAL := hh_approximated.svg

BUILD_CATALOGUE_SCRIPT := build-catalogue

hh_approx-catalogue.so: $(wildcard mechanisms/*.mod)
	$(BUILD_CATALOGUE_SCRIPT) hh_approx mechanisms

hh_approximated.svg: hh_approx-catalogue.so
	./arbor_hh_approximated.py --save hh_approximated.svg --catalogue ./$<

.PHONY: clean
clean:
	rm -f *svg *so
