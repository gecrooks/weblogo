#!/bin/sh


echo 
echo "## Build a WebLogo release" 
echo "# Use setup.py to install from a built distribution."
echo "# Requires epydoc to build documentation. "
python -c "import weblogolib; print weblogolib.description"  

echo 
echo "## Check Subversion status (Working copy checked in and up to date?) " 
svn status 

echo 
echo "## Lines of code" 
wc -l *.py weblogolib/*.py corebio/*.py corebio/*/*.py test_corebio/*.py weblogo transformseq | grep 'total' 

echo 
echo "## Code Tags" 
grep 'FIXME\|TODO' *.py weblogolib/*.* weblogolib/htdocs/*.* 


echo
echo "# Cleaning previous "
rm -rd dist/_extract_/


echo 
echo "## Rebuild API docs :" 
epydoc -q -o apidocs/ -n WebLogo -u http://code.google.com/p/weblogo/ --docformat plaintext --no-frames --no-private weblogolib    || exit

echo 
echo "## Rebuild examples :" 
cd weblogolib/htdocs/examples
sh build_examples.sh
cd ../../..     ls


echo
echo "## CoreBio version"
python -c 'import corebio; print corebio.__version__'

echo
echo "## Build source distribution :"
python setup.py -q sdist                                     || exit

echo
echo "## Extract source distribution"
mkdir -p dist/_extract_
tar -zxf dist/weblogo*.tar.gz -C dist/_extract_                 || exit

echo "## Change directory to source"
cd dist/_extract_/weblogo* 

echo "## Run unit tests  :"
python ./test_weblogo.py                                     || exit 
python ./test_corebio.py                                     || exit 

echo
echo "## Build "
python setup.py build                                        || exit

echo
echo "## Install "
python setup.py install --home ../install                    || exit
echo

echo "## Clean up..."
cd ../../..
rm -rd dist/_extract_/                                          || exit

echo
echo "## Run build test  :"
sh build_test.sh                                                || exit 


echo "## Success"

