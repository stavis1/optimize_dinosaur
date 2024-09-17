FROM r-base:4.4.1

RUN mkdir /data
RUN mkdir /xcms

COPY xcms_quantify_features.R /xcms/xcms_quantify_features.R

RUN apt-get update -y && \
    apt-get install -y libcurl4-openssl-dev && \
    apt-get install -y libnetcdf-dev

RUN Rscript -e 'install.packages(c("optparse", "configr"))'
RUN Rscript /xcms/xcms_quantify_features.R --install

