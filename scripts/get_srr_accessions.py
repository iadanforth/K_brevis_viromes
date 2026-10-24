#!/usr/bin/env python3

import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

accessions = """
SRX23656680
SRX23656682
SRX23656683
SRX23656684
SRX23656667
SRX23656668
SRX23656669
SRX23656670
SRX23656671
SRX23656672
SRX23656673
SRX23656674
SRX23656675
SRX23656676
SAMN50254270
SAMN50254271
""".split()

base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def fetch_xml(url):
    max_attempts = 6

    for attempt in range(max_attempts):
        # Keep requests spaced out, including retries.
        time.sleep(1.1)

        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "SRA-accession-lookup/1.0"},
            )
            with urllib.request.urlopen(request, timeout=30) as response:
                return ET.fromstring(response.read())

        except urllib.error.HTTPError as error:
            retryable = error.code == 429 or 500 <= error.code < 600
            if not retryable or attempt == max_attempts - 1:
                raise

            retry_after = error.headers.get("Retry-After")
            try:
                delay = float(retry_after)
            except (TypeError, ValueError):
                delay = min(60, 5 * (2 ** attempt))

            print(
                f"NCBI returned HTTP {error.code}; retrying in {delay:g}s",
                file=sys.stderr,
                flush=True,
            )
            time.sleep(delay)


for accession in accessions:
    search_params = urllib.parse.urlencode({
        "db": "sra",
        "term": accession,
        "retmode": "xml",
    })
    search_url = f"{base_url}esearch.fcgi?{search_params}"

    search_xml = fetch_xml(search_url)
    ids = [
        node.text
        for node in search_xml.findall(".//Id")
        if node.text
    ]

    if not ids:
        print(f"{accession}\tNO SRA RECORD FOUND", flush=True)
        continue

    fetch_params = urllib.parse.urlencode({
        "db": "sra",
        "id": ",".join(ids),
        "retmode": "xml",
    })
    fetch_url = f"{base_url}efetch.fcgi?{fetch_params}"

    fetch_result = fetch_xml(fetch_url)
    runs = sorted({
        node.attrib["accession"]
        for node in fetch_result.findall(".//RUN")
        if "accession" in node.attrib
    })

    print(
        f"{accession}\t{','.join(runs) if runs else 'NO RUN FOUND'}",
        flush=True,
    )
