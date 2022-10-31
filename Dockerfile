FROM continuumio/miniconda3 AS build
#Copy environment specification to current directory
COPY conda.yml requirements.txt pyproject.toml setup.py README.md ./
COPY couplet ./couplet/
COPY bin ./bin/

#create conda environment
RUN conda env create  --file conda.yml --name couplet
# push conda activate command in bashrc to use the new environment:
RUN echo "conda activate couplet" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]
RUN echo "couplet activated"
#setup utilizing setup.py in couplet 
RUN python3 ./setup.py install
RUN echo "conda deactivate " >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]
RUN echo "couplet deactivated"
# Install conda-pack:
RUN conda install -c conda-forge conda-pack
# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda-pack -n couplet -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

# We've put venv in same path it'll be in final image,
# so now fix up paths:
RUN /venv/bin/conda-unpack

# The runtime-stage image; we can use ubuntu as the
# base image since the Conda env also includes Python
# for us.
FROM ubuntu:latest AS runtime
RUN apt-get update -y && apt-get install -y procps
# Copy /venv from the previous stage:
COPY --from=build /venv /venv
COPY bin/*.py /venv/bin/

#When image is run, run the code with the environment
# activated:
SHELL ["/bin/bash", "-c"]

 #          python3 -c "import couplet; print( couplet.__file__);  print('success!')" 

#export path to /venv/bin/
ENV PATH="/venv/bin:${PATH}"
