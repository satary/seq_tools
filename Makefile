
clean:
	rm -rf build
	rm -rf dist
	rm -rf *.egg-info

sdist:
	rm -rf dist
	rm -rf *.egg-info
	python setup.py sdist

pypi-push:
	twine upload dist/*

gh-push:
	git add .
	git commit -m "Update"
	git push

upd: gh-push sdist pypi-push