# top-level pyquil Makefile

PACKAGENAME = weblogo
FILES = weblogo setup.py
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

lint:		## Delint python source
	@flake8 flake8 setup.py weblogo tests

typecheck:	## Static typechecking 
	@mypy $(FILES)  --ignore-missing-imports --follow-imports=skip

untyped:	## Report type errors and untyped functions
	@mypy $(FILES) --ignore-missing-imports --follow-imports=skip --disallow-untyped-defs

docs:		## Build documentation
	$(MAKE) -C docs html

.PHONY: help
.PHONY: docs