========================
multicell
========================

.. {# pkglts, doc


.. image:: https://badge.fury.io/py/multicell.svg
    :alt: PyPI version
    :target: https://badge.fury.io/py/multicell

.. #}


Python library to run biological simulations of 3D, multicellular tissues.

----------------------------------------------
Installation instructions for Ubuntu 16.04.1
----------------------------------------------

Install Ubuntu and update it. ::

    sudo apt update
    sudo apt upgrade
    
Install required packages. ::
    
    sudo apt-get install git gitk python-pip scons zlib1g-dev libpng12-dev python-qt4 python-qt4-gl qt4-dev-tools libreadline-dev freeglut3 bison flex g++ libqt4-dev libqt4-opengl-dev libqhull-dev libreadline-dev pkg-config freeglut3-dev libann-dev liblapack-dev libmpfr-dev libblas-dev freeglut3-dev libboost-all-dev libeigen2-dev pyqt4-dev-tools python-sip-dev python-tk
    
Update pip (a Python module manager). ::

    pip install -U pip
    
Install Python modules with ``pip``. ::
    
    pip install ipython==5.5 numpy scipy matplotlib jupyter pandas soappy scikit-image scikit-learn sphinx --user
    
Install ``deploy`` and ``sconsx``. ::
    
    git clone https://github.com/openalea/deploy.git
    cd deploy
    git checkout 9ab468b5316ecdddcf38e17fce901cba09c018b5
    python setup.py develop --user
    cd ..
    git clone https://github.com/openalea/sconsx.git
    cd sconsx
    git checkout 1ebd690c65a87d41c5be4eb68b51e6451fb27993
    python setup.py develop --user
    cd ..
    
Edit the ``.pam-environment`` file, which contains a definition of environment variables for your session. ::
    
    gedit .pam-environment
    
Paste the following if the file is empty, otherwise, edit accordingly. ::
    
    PATH DEFAULT=${PATH}:${HOME}/.local/bin
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/.local/bin

Save, log out of your session, and log back in for the changes to take effect.

Install ``openalea``. ::

    git clone https://github.com/openalea/openalea.git
    cd openalea
    git checkout 7986945878f7fb78095265d5d9b97c32fbcc1f97
    python multisetup.py develop --user
    cd ..
    git clone https://github.com/openalea/openalea-components.git
    cd openalea-components
    git checkout 4f86e0a0a157ea9428a348f11ffb53e012122a69
    python multisetup.py develop --user
    cd ..
    
Download ``vplants``. ::
    
    wget https://github.com/jldinh/vplants/archive/refs/tags/v1.2.0.tar.gz
    tar xzvf vplants-1.2.0.tar.gz
    cd vplants-1.2.0
    
Edit ``multisetup.py`` (e.g. ``gedit multisetup.py``): in the `dirs` section, comment every line (add a leading #), except ‘tissue’ (remove the leading #), save and close. Install vplants. ::

    python multisetup.py develop --user
    cd ..
    
Install ``odesparse``, a numerical integrator. ::
    
    wget http://www.simuplant.org/files/simuplant-server.zip
    unzip simuplant-server.zip
    cd RootGUIServer/odesparse_1
    python setup.py install --user
    cd ../..
    
Finally, install ``multicell``. ::
    
    git clone https://github.com/jldinh/multicell.git
    cd multicell
    python setup.py develop --user

To start the examples, run::

    jupyter notebook
    
Browse to the ``examples`` folder, inside of the ``multicell`` folder, and open one of the files.
