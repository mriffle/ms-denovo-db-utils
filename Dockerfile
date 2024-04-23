FROM python:3.10

ADD python_scripts/process_casanovo_results.py /usr/local/bin/process_casanovo_results.py
ADD python_scripts/process_comet_results.py /usr/local/bin/process_comet_results.py
ADD entrypoint.sh /usr/local/bin/entrypoint.sh

RUN apt-get update && \
    apt-get install -y procps && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD []
