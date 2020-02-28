# Specify OS (Ubuntu 16.04 'Xenial Xerus')
FROM ubuntu:16.04

RUN apt-get update

# Install latest version of R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get install -y software-properties-common
RUN apt-get install apt-transport-https
RUN add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
RUN apt-get update
RUN apt-get install -y r-base
RUN apt-get install -y r-recommended
RUN apt-get install -y r-base-html
RUN apt-get install -y r-base

RUN apt-get install -y git
RUN apt-get install -y curl
RUN apt-get install -y apt-utils
RUN apt-get install -y libcurl4-gnutls-dev
RUN apt-get install -y libssl-dev
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y libxml2-dev

ADD install-dada.R .
RUN ./install-dada.R

RUN mkdir -p /opt/dada2
WORKDIR /opt/dada2
COPY run_dada2.R /opt/dada2/

ENTRYPOINT ["Rscript", "/opt/dada2/run_dada2.R"]
CMD ["--help"]
