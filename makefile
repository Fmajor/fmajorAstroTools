all:
	python3 setup.py build
	python3 setup.py develop
local:
	python3 setup.py build
	python3 setup.py develop --user
linux:
	python3 setup.py clean
	python3 setup.py build
	sudo python3 setup.py develop
clean:
	python3 setup.py clean
	#python3 setup.py install --record files3.txt
	#cat files3.txt | xargs rm -rf
