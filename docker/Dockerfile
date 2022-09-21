FROM python:3.7-buster AS dev

RUN pip install poetry

COPY pyproject.toml /app/build/ukb-pret/
COPY poetry.toml /app/build/ukb-pret/
COPY poetry.lock /app/build/ukb-pret/
COPY src /app/build/ukb-pret/src

WORKDIR /app/build/ukb-pret
RUN poetry build

FROM python:3.7-buster AS local

WORKDIR /app
COPY --from=dev /app/build/ukb-pret/dist/ukb_pret*.whl /app

RUN pip install ukb_pret*.whl
