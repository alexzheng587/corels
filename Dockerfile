# commands to run:

# sudo docker run -it -v src:/src -v data:/data corels (ON LOCAL MACHINE)
# scl enable devtoolset-7 bash (IN DOCKER IMAGE) 

FROM centos:7.5.1804

RUN yum -y install centos-release-scl

RUN yum-config-manager --enable rhel-server-rhscl-7-rpms

RUN yum -y install devtoolset-7 gmp-devel

RUN scl enable devtoolset-7 bash

#COPY src/ /src
#COPY data /data

WORKDIR /src

ENTRYPOINT ["/bin/bash"]
