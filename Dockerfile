FROM ubuntu:18.04 as intermediate
LABEL stage=intermediate

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git
#RUN mkdir /lpsubpclones && cd /lpsubpclones && git clone https://github.com/supark87/lpsubp.git 
RUN git clone https://github.com/supark87/lpsubp.git 

FROM ubuntu:18.04
COPY --from=intermediate /lpsubp/ /data/
WORKDIR /data/
#RUN git clone https://github.com/supark87/lpsubp.git
#COPY . /data/
RUN apt-get update && apt-get install -y ncbi-blast+=2.6.0-1
RUN apt-get update && apt-get install -y muscle
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.7 \
    python3-pip \
    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install matplotlib 
RUN pip3 install biopython
RUN pip3 install wheel
RUN pip3 install ipython
RUN pip3 install hdbscan==0.8.13
RUN pip3 install -r /data/requirements1.txt

#CMD ["python3", "/data/pipe.py","/data/test2/"]

ENTRYPOINT ["python3","/data/typer_test_ver1.py"]
