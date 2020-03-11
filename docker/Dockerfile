# Base Image
FROM ubuntu:18.04

# Metadata
LABEL base.image="ubuntu:18.04"
LABEL software="bayesian_hdx"
LABEL software.version="v2.1"
LABEL description="Bayesian Analysis of HDX-MS Data."
LABEL website="https://github.com/salilab/bayesian_hdx"
LABEL documentation="https://github.com/salilab/bayesian_hdx"
LABEL license="GNU Lesser General Public License (LGPL)"
LABEL tags="HDX"

# Maintainer
MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

ENV URL=https://github.com/salilab/bayesian_hdx
ENV DST=/usr/local/

USER root

RUN apt-get clean all && \
    apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y apt-utils && \
    apt-get install -y tzdata && \
    apt-get install -y software-properties-common && \
    apt-get install -y git python3 python3-pip && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN ln -sfn /usr/bin/python3.6 /usr/bin/python
RUN mkdir $DST

RUN cd /$DST/ && \
    git clone $URL && \
    chmod a+x $DST/bin/* && \
    pip3 install -r $DST/requirements.txt

RUN adduser --disabled-password --gecos '' ubuntu
RUN chmod a+rwx /home/ubuntu/
RUN mkdir /home/ubuntu/bin
RUN chown -R ubuntu /home/ubuntu
USER ubuntu
ENV PYTHONPATH=${DST}/pyext/src
ENV PATH=${DST}/bin:$PATH

WORKDIR /data/
