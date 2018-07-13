#!/bin/sh

cp -r $SRC_DIR/src/*.py $PREFIX/bin
ln -s $PREFIX/bin/NGSweep.py $PREFIX/bin/ngsweep
chmod +x $PREFIX/bin/ngsweep
