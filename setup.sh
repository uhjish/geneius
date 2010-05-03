echo "Removing artifacts..."
rm libgeneius/seq/*.so
echo "Compiling C-language dependencies"
python setup.py build_ext --inplace
echo "Done. Just enjoy..."
