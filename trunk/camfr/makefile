camfr: FORCE
	scons

test: FORCE
	cd testsuite
	make

FORCE:

clean:
	rm -f *~ *.pyc core MANIFEST
	python setup.py clean
	rm -f -R build
	rm -f -R dist
	scons -c
	cd examples ; make clean
	cd visualisation ; make clean
	cd testsuite ; make clean
	cd docs ; make clean
