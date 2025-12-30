
.DEFAULT_GOAL := help

NAME = weblogo
FILES = $(NAME) tests docs/conf.py setup.py

USER = -i ~/.ssh/aws_kaiju.pem ec2-user
HOST = weblogo.threeplusone.com
DIR = /home/ec2-user/weblogo


# Kudos: Adapted from Auto-documenting default target
# https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-12s\033[0m %s\n", $$1, $$2}'


.PHONY: test
test:  ## Run unittests
	uv run python -m pytest --disable-pytest-warnings

.PHONY: coverage
coverage:  ## Report test coverage
	@echo
	uv run python -m pytest --disable-pytest-warnings --cov-report term-missing --cov $(NAME)
	@echo

.PHONY: lint
lint:  ## Lint check python source
	@uv run ruff check

.PHONY: delint
delint:
	@echo
	uv run ruff format

.PHONY: typecheck
typecheck:  ## Static typechecking
	uv run mypy $(NAME)

.PHONY: docs
docs:  ## Build documentation
	uv run sphinx-build -b html docs docs/_build/html

.PHONY: docs-open
docs-open:  ## Build documentation and open in webbrowser
	uv run sphinx-build -b html docs docs/_build/html
	open docs/_build/html/index.html

.PHONY: examples
examples:  ## Build examples
	cd weblogo/htdocs/examples && ./build_examples.sh

.PHONY: login
login:  ## ssh into remote server
	ssh  ${USER}@${HOST}

.PHONY: sync
sync:  ## Sync remote server
	ssh ${USER}@${HOST} 'cd weblogo && git pull'


.PHONY: restart
restart:  	## Restart remote server
	ssh ${USER}@${HOST} 'sudo systemctl restart httpd.service'
