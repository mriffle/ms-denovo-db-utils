FROM python:3.10

ADD python_scripts/process_casanovo_results.py /usr/local/bin/process_casanovo_results.py
ADD python_scripts/process_comet_results.py /usr/local/bin/process_comet_results.py

ADD entrypoint.sh /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD []
