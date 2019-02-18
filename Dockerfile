FROM centos7.5.1804

RUN yum install gmp-devel gcc

COPY src/ /src
COPY data /data

WORKDIR /src

RUN make clean && make

ENTRYPOINT ["/bin/sh"]
