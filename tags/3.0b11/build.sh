#!/bin/sh


echo 
echo "## Build WebLogo release" 
python -c "import weblogolib; print weblogolib.description"  

echo 
echo "## Check Subversion status (Working copy checked in and up to date?) " 
svn status 

echo 
echo "## Lines of code" 
wc -l *.py weblogolib/*.py | grep 'total' 

echo "# Cleaning previous "
rm -rd dist/_extract_/


echo 
echo "## Build API docs :" 
epydoc -q -o apidocs/ -n WebLogo -u http://code.google.com/p/weblogo/ --docformat plaintext --no-frames --no-private weblogolib    || exit



echo
echo "## PYTHON 2.4 ##"
echo "## Build source distribution :"
python2.4 setup.py -q sdist                                     || exit

echo
echo "## Extract source distribution"
mkdir -p dist/_extract_
tar -zxf dist/weblogo*.tar.gz -C dist/_extract_                 || exit

echo "## Change directory to source"
cd dist/_extract_/weblogo* 

echo "## Run unit tests  :"
python2.4 ./test_weblogo.py                                     || exit 

echo
echo "## Build "
python2.4 setup.py build                                        || exit

echo
echo "## Install "
python2.4 setup.py install --home ../install                    || exit
echo

echo "## Clean up..."
cd ../../..
rm -rd dist/_extract_/                                          || exit

echo "## Success"

