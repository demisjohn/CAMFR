# Simple shell script to replace camfr_work by camfr in all
# the files in a given directory.
#
# To change the module name, the following needs to be changed:
#
#  imports in testsuite
#  imports in camfr/geometry.py
#  imports in camfr/visualisation
#  setup.py: name, extra_path

for name in *
do
  sed s/camfr_work/camfr/g <$name >tmpfile
  mv tmpfile $name 
done

rm -f tmpfile