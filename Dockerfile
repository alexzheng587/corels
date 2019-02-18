FROM centos:7.5.1804

RUN yum -y install gmp-devel gcc make gcc-c++

COPY src/ /src
COPY data /data

WORKDIR /src

ENTRYPOINT ["/bin/sh"]
