camfr: FORCE
	python setup.py build

install:
	python setup.py install

test: FORCE
	cd testsuite
	make

distrib:
	rm -f *.tgz Exclude
	cd docs ; make pdf
	rm -f -r dist ../camfr_dist camfr_dist camfr_*
	mkdir ../camfr_dist
	cp -r * ../camfr_dist
	mv ../camfr_dist .
	cd camfr_dist ; make clean
	rm -r camfr_dist/docs/*
	cp docs/camfr.pdf camfr_dist/docs
	V=`python camfrversion.py` && mv camfr_dist camfr_$${V}
	V=`python camfrversion.py` && find camfr_$${V} -type d -print | egrep '/,|%$$|~$$|CVS|build' > Exclude
	V=`python camfrversion.py` && find camfr_$${V} ! -type d -print | egrep '/,|%$$|\#|~$$|dblite|pyc|\.old$$|/core$$|\.orig$$' >> Exclude
	V=`python camfrversion.py` && tar cvfzX camfr_$${V}.tgz Exclude camfr_$${V}
	V=`python camfrversion.py` && rm -r camfr_$${V}
	rm -f Exclude 

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
