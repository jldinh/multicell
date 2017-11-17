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

Install Ubuntu and update it.

::

    sudo apt update
    sudo apt upgrade
    
Install required packages.

::
    
    sudo apt-get install git gitk python-pip scons zlib1g-dev libpng12-dev python-qt4 python-qt4-gl qt4-dev-tools libreadline-dev freeglut3 bison flex g++ libqt4-dev libqt4-opengl-dev libqhull-dev libreadline-dev pkg-config freeglut3-dev libann-dev liblapack-dev libmpfr-dev libblas-dev freeglut3-dev libboost-all-dev libeigen2-dev pyqt4-dev-tools python-sip-dev python-tk
    
Install Python modules.

::
    
    pip install -U pip
    pip install ipython numpy scipy matplotlib jupyter pandas soappy scikit-image scikit-learn sphinx --user
    git clone https://github.com/openalea/deploy.git
    cd deploy
    python setup.py develop --user
    cd ..
    git clone https://github.com/openalea/sconsx.git
    cd sconsx
    python setup.py develop --user
    cd ..
    nano .pam-environment
    PATH DEFAULT=${PATH}:${HOME}/.local/bin
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/.local/bin
    Ctrl+X, Y, Enter
    Relog
    git clone https://github.com/openalea/openalea.git
    cd openalea
    python multisetup.py develop --user
    cd ..
    git clone https://github.com/openalea/openalea-components.git
    cd openalea-components
    python multisetup.py develop --user
    cd ..
    wget https://gforge.inria.fr/frs/download.php/34177/vplants-1.2.0.tar.gz
    tar xzvf vplants-1.2.0.tar.gz
    cd vplants
    
Edit ``multisetup.py`` (e.g. ``gedit multisetup.py``): comment ‘container’ (add a leading #), uncomment ‘tissue’ (remove it), save and close.

::

    python multisetup.py develop --user
    wget http://www.simuplant.org/files/simuplant-server.zip
    unzip simuplant-server.zip
    cd RootGUIServer/odesparse_1
    python setup.py install --user
    cd ..
    git clone https://github.com/jldinh/multicell.git
    cd multicell
    python setup.py develop --user
