FROM ubuntu:16.10

LABEL maintainer "alanh@unc.edu"


# ...put your own build instructions here...

RUN apt-get update && apt-get install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev r-base git

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://archive.linux.duke.edu/cran'; options(repos = r);" > ~/.Rprofile 
RUN Rscript -e 'install.packages (c("mclust","flsa","foreach","ggplot2","gridExtra","inline","Rcpp","devtools","optparse")); ' \
    && Rscript -e 'source("https://bioconductor.org/biocLite.R");biocLite(); biocLite("DNAcopy");' \
    && Rscript -e 'library(devtools) ; install_github("ChenMengjie/SynthEx");'


ADD https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz bedtools-2.26.0.tar.gz
RUN tar -zxf bedtools-2.26.0.tar.gz && cd bedtools2 && make && make install


# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /bedtools*

COPY Rprofile.site /etc/R/
COPY SynthEx.R /SynthEx.R

ENTRYPOINT ["Rscript", "/SynthEx.R"]
CMD ["--help"]
