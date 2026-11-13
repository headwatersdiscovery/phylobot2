.PHONY: website website-serve docker-cpus ml-image ml-test no-running-jobs golden-full grasp-image clean e2e submit env up down logs ps storage-init test lint fmt typecheck check golden-generate pipeline golden tools-image

# The runner mounts the project at its real path; see docker-compose.yml.
export PROJECT_ROOT := $(CURDIR)

# Settings require a secret key even for host-side tooling (mypy's Django plugin).
HOST_ENV = DJANGO_SECRET_KEY=host-tooling-not-secret PHYLOBOT_MODE=local

env:  ## Create .env from the example if missing; check it has every variable the example has
	@test -f .env || (cp .env.example .env && chmod 600 .env && echo "Created .env from .env.example")
	@./scripts/env-secrets.sh
	@./scripts/check-env.sh
	@./scripts/docker-cpus.sh --check

docker-cpus:  ## Show Docker's VM CPUs and the PHYLOBOT_DOCKER_CPUSET to set in .env
	@./scripts/docker-cpus.sh

up: no-running-jobs env tools-image ml-image  ## Build and start the full local stack (refuses while jobs run; FORCE=1)
	@# Storage first: web's health check needs the bucket, and worker and searches wait for
	@# a healthy web, so on a first install the bucket must exist before web starts.
	docker compose up -d --wait db garage
	./scripts/garage-init.sh
	docker compose up --build -d

no-running-jobs:  # make up recreates the worker, which stops its jobs
	@./scripts/check-no-jobs.sh
	@echo "PhyloBot: http://localhost:8000  (health: http://localhost:8000/healthz)"

down:  ## Stop the stack (data volumes are kept)
	docker compose down

logs:
	docker compose logs -f

submit:  ## Queue a job through the web stack: make submit IN=file.fasta [ARGS="--name x"]
	@test -n "$(IN)" || (echo "usage: make submit IN=file.fasta [ARGS=...]" && exit 1)
	docker compose cp $(IN) web:/tmp/submit.fasta
	docker compose exec web python manage.py submit_job /tmp/submit.fasta $(ARGS)

ps:
	docker compose ps

storage-init:  ## Garage setup: layout, bucket, access key (make up runs it; safe to repeat)
	./scripts/garage-init.sh

test:  ## All tests except golden, in a container against the compose Postgres
	docker compose run --rm --build test

lint:
	uv run ruff check .
	uv run ruff format --check .

fmt:
	uv run ruff check --fix .
	uv run ruff format .

typecheck:
	$(HOST_ENV) uv run mypy src tests manage.py

golden-generate:  ## Simulate golden datasets (skips existing; owner sign-off to regenerate)
	$(HOST_ENV) uv run python tests/golden/generate.py $(if $(FORCE),--force)

RUNNER = docker compose run --rm runner

tools-image: grasp-image  ## Build phylobot-tools (our pipeline steps); amd64 like the tool images
	docker build --platform linux/amd64 -f docker/tools.Dockerfile -t phylobot-tools:local .

ml-image:  ## Build phylobot-ml (tools + torch, metapredict, ESM-2 weights)
	docker build --platform linux/amd64 -f docker/ml.Dockerfile -t phylobot-ml:local .

ml-test:  ## The neural property predictors' tests, inside the ml image (needs torch and weights)
	docker build --platform linux/amd64 -f docker/ml.Dockerfile --target ml-test -t phylobot-ml:test .
	docker run --rm --platform linux/amd64 -v $(CURDIR)/tests:/opt/src/tests:ro -w /opt/src \
		phylobot-ml:test pytest -p no:django -p no:cacheprovider -q tests/test_properties_ml.py

grasp-image:  ## Build phylobot-grasp (downloads bnkit.jar, checks its SHA-256; never push it)
	docker build --platform linux/amd64 -f docker/grasp.Dockerfile -t phylobot-grasp:2.0.0 docker

# make pipeline IN=input.fasta OUT=results/name [PARAMS=params.json] [ARGS="--aligners mafft"]
#               [WORK=dir]  Nextflow scratch (default: work/, which grows with every run;
#                           `make clean` empties it). Golden tests use their own and
#                           delete it when they pass.
pipeline: tools-image ml-image  ## Run the pipeline (IN and/or PARAMS; see comment above)
	@./scripts/docker-cpus.sh --check
	@test -n "$(IN)$(PARAMS)" || (echo "usage: make pipeline IN=input.fasta OUT=results/name [PARAMS=params.json] [ARGS=...]" && exit 1)
	$(RUNNER) nextflow run pipeline/main.nf -profile docker \
		$(if $(IN),--input $(IN)) $(if $(PARAMS),-params-file $(PARAMS)) \
		--outdir $(or $(OUT),results/$(basename $(notdir $(or $(IN),$(PARAMS))))) \
		$(if $(WORK),-work-dir $(WORK)) $(ARGS)

golden:  ## Golden tests, quick tier: the first replicate of each scenario
	$(HOST_ENV) uv run pytest -m "golden and not golden_full" -s tests/golden

golden-full:  ## Golden tests, all replicates, plus per-scenario averages and method comparison
	$(HOST_ENV) uv run pytest -m golden -s tests/golden

e2e: up  ## Drive the running stack end to end over HTTP (submits, follows, deletes a job)
	$(HOST_ENV) uv run pytest -m e2e -s

clean:  ## Remove local scratch (work/, results/, Nextflow logs); never job data or volumes
	./scripts/clean.sh

# phylobot.com (website/README.md): one static page built from the app's own About text,
# skins and backdrop. Runs in the web image, as you, so website/dist/ is yours.
WEBSITE = docker compose run --rm --no-deps --user "$$(id -u):$$(id -g)" -e PYTHONDONTWRITEBYTECODE=1 \
	-v $(CURDIR)/website:/app/website

website:  ## Build phylobot.com into website/dist/ [ARGS="--youtube-id <id>" to preview the video]
	$(WEBSITE) web python website/build.py $(ARGS)

website-serve: website  ## Build it, then serve it at http://localhost:8001 until Ctrl-C
	@echo "phylobot.com preview: http://localhost:8001"
	$(WEBSITE) -p 127.0.0.1:8001:8001 web python -m http.server 8001 --bind 0.0.0.0 --directory website/dist

check: lint typecheck test  ## Everything CI would run
