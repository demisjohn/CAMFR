camfr: FORCE
	cd camfr ; make

test:
	cd testsuite ; make

FORCE:

clean:
	rm -f *~ core MANIFEST
	python setup.py clean
	rm -R build
	rm -R dist
	cd camfr ; make clean
	cd examples ; make clean
	cd visualisation ; make clean
	cd testsuite ; make clean
	cd docs ; make clean
