
pip install numpy --no-index
pip install scipy --no-index
pip install astropy --no-index
pip install emcee --no-index
pip install dynesty
pip install h5py --no-index

cd
gedit .bashrc
# ADD:'export SPS_HOME="/home/clawlorf/fsps/"' at end of file

cd
git clone https://github.com/cconroy20/fsps.git $SPS_HOME # v3.2

cd
python -m pip install fsps # v0.4.1

cd
git clone https://github.com/bd-j/sedpy
# now copy all the 'hff_f*.par' files from /home/clawlorf/HFF/hff_filters/ into sedpy/sedpy/data/filters/
cd sedpy
python -m pip install . # v0.2.1

cd
git clone https://github.com/bd-j/prospector
# now edit the constants in prospector/prospect/sources/constants.py, particularly the cosmology
cd prospector
python -m pip install . # v1.0.0, but with updates for plotting SEDs with the faster sedpy FilterSets


## See here for more resources: ##
To install propsector on ComputeCanada resource:
[see https://docs.computecanada.ca/wiki/Python]
[also this for virtualenv info: https://dev.to/bowmanjd/python-tools-for-managing-virtual-environments-3bko#howto]
