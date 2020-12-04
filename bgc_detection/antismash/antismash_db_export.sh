#!/usr/bin/env bash
# Download antiSMASH clusters based on accession ids
# Input: one accession ID (or part of ID) per line
# Will include partial accession matches, so it's possible that additional results will be obtained

echoerr() { cat <<< "$@" 1>&2; }

set -e

while read line
do
    echoerr 'Fetching' $line;
    curl 'https://antismash-db.secondarymetabolites.org/api/v1.0/export' \
            --fail \
            --retry 10 \
            --retry-delay 5 \
            -H 'Accept-Encoding: gzip, deflate, br' \
            -H 'accept: text/csv' \
            -H 'Content-Type: application/json;charset=UTF-8' \
            -H 'Connection: keep-alive' \
            --data-binary "{\"query\":{\"search\":\"cluster\",\"return_type\":\"csv\",\"terms\":{\"term_type\":\"expr\",\"category\":\"acc\",\"term\":\"$line\"}}}" \
            --compressed;
done < "${1:-/dev/stdin}"
