FROM python:3.7

LABEL name "foundation-xml-cnv"
LABEL version "1.0.0"
LABEL maintainer "LifeOmic <development@lifeomic.com>"

RUN mkdir -p /opt/app
WORKDIR /opt/app
COPY . /opt/app
RUN pip install -r requirements.txt

ENTRYPOINT ["python", "/opt/app/src/convert.py"]
