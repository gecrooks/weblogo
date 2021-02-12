# top-level pyquil Makefile

NAME = weblogo
FILES = weblogo setup.py tests

USER = ec2-user
HOST = weblogo.threeplusone.com
DIR = /home/ec2-user/weblogo


# Kudos: Adapted from Auto-documenting default target 
# https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-12s\033[0m %s\n", $$1, $$2}'

test:		## Run unittests with current enviroment
	@pytest 

testall:	## Run full test suite
	@tox

coverage:	## Report test coverage
	@pytest --cov=weblogo --cov-report term-missing

lint:		## Lint check python source
	@isort --check $(FILES)  ||  echo "isort:   FAILED!"
	@black --check --quiet $(FILES)    || echo "black:   FAILED!"
	@flake8 $(FILES)

delint:   ## Run isort and black to delint project
	@echo
	isort $(FILES)
	@echo
	black $(FILES)
	@echo

typecheck:	## Static typechecking 
	@mypy $(FILES) --ignore-missing-imports --follow-imports=skip

untyped:	## Report type errors and untyped functions
	@mypy $(FILES) --ignore-missing-imports --follow-imports=skip --disallow-untyped-defs

docs:		## Build documentation
	$(MAKE) -C docs html && open docs/_build/html/index.html

examples: 	## Build examples
	cd weblogo/htdocs/examples && ./build_examples.sh

login:		## ssh into remote server
	ssh  ${USER}@${HOST}

sync:		## Sync remote server
	ssh ${USER}@${HOST} 'cd weblogo && git pull'

reinstall:	## Reinstall weblogo on remote server
	# pip install needed to update versions with latest git tag
	ssh ${USER}@${HOST} 'cd weblogo && sudo pip3 install -e .'

restart:  	## Restart remote server
	ssh ${USER}@${HOST} 'sudo systemctl restart httpd.service'


.PHONY: help
.PHONY: docs