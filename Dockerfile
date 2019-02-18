FROM centos:7.5.1804

RUN yum -y install gmp-devel make centos-release-scl devtoolset-7-gcc-c++ gdb vim

COPY src/ /src
COPY data /data

WORKDIR /src

ENTRYPOINT ["scl enable devtoolset-7 bash"]
