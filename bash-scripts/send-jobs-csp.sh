#!/bin/bash

cd MgB2H8
echo -e "1\nyes" | python3 PyMCSP.py
echo -e "2\nyes" | python3 PyMCSP.py
cd ..

cd Zr2Co11
echo -e "1\nyes" | python3 PyMCSP.py
cd ..

cd UH7
echo -e "1\nyes" | python3 PyMCSP.py
echo -e "2\nyes" | python3 PyMCSP.py
cd ..

