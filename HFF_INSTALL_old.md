
pip install numpy --no-index
pip install scipy --no-index
pip install astropy --no-index
pip install emcee --no-index
pip install dynesty
pip install h5py --no-index

cd
git clone https://github.com/cconroy20/fsps
cd
gedit .bashrc
# ADD:'export SPS_HOME="/home/cam/hff/fsps/"' at end of file
cd
cd $SPS_HOME/src
make clean
make all

cd
git clone https://github.com/dfm/python-fsps
cd python-fsps
python setup.py install

cd
git clone https://github.com/bd-j/sedpy
cd sedpy
python setup.py install

cd
git clone https://github.com/bd-j/prospector
cd prospector
python setup.py install

## See here for more resources: ##
To install propsector on ComputeCanada resource:
[see https://docs.computecanada.ca/wiki/Python for initialization]
[also this for virtualenv info: https://dev.to/bowmanjd/python-tools-for-managing-virtual-environments-3bko#howto]
