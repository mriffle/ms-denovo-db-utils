FROM python:3.10-slim

ADD python_scripts/process_casanovo_results.py /usr/local/bin/process_casanovo_results.py
ADD python_scripts/process_comet_results.py /usr/local/bin/process_comet_results.py
ADD python_scripts/collate_into_fasta.py /usr/local/bin/collate_into_fasta.py
ADD python_scripts/build_reset_input.py /usr/local/bin/build_reset_input.py
ADD generate_reverse_decoys.py /usr/local/bin/generate_reverse_decoys.py
ADD bash_scripts/split_fasta.sh /usr/local/bin/split_fasta.sh

ADD entrypoint.sh /usr/local/bin/entrypoint.sh

RUN apt-get update && \
    apt-get install -y procps && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean && \
    chmod 755 /usr/local/bin/entrypoint.sh && \
    chmod 755 /usr/local/bin/split_fasta.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD []
