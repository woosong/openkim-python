all: build

$(KIM_DIR)/KIM_API/python:
	mkdir $(KIM_DIR)/KIM_API/python

install_kim: $(KIM_DIR)/KIM_API/python
	cp -r tests $(KIM_DIR)/EXAMPLE_PYTHON
	cp -r virial $(KIM_DIR)/EXAMPLE_PYTHON/test_Cu_FCC_virial
	cp *.i *.c *.h *.py Makefile README $(KIM_DIR)/KIM_API/python

build:
	python setup.py build

install:
	python setup.py install
