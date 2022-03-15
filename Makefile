.PHONY: build dist redist install install-from-source clean uninstall

build:
	CYTHONIZE=1 ./setup.py build

install:
	CYTHONIZE=1 pip install .

clean:
	$(RM) -r build dist src/*.egg-info
	$(RM) -r src/xfeat/{basic.cpp,main.cpp}
	$(RM) -r .pytest_cache
	find . -name __pycache__ -exec rm -r {} +
	#git clean -fdX

uninstall:
	pip uninstall xfeat
