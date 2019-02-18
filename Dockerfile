FROM centos:7.5.1804

RUN yum -y install gmp-devel make centos-release-scl gdb vim
RUN yum -y install devtoolset-7-gcc-c++

COPY src/ /src
COPY data /data

WORKDIR /src

CMD scl enable devtoolset-7 /bin/bash
