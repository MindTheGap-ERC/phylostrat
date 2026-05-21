#!/bin/bash

# Usage:
# ./download_yoda.sh https://doi.org/10.xxxx/xxxxx

DOI_URL="$1"

if [ -z "$DOI_URL" ]; then
    echo "Usage: $0 <DOI_URL>"
    exit 1
fi

OUTDIR="yoda_download"

mkdir -p "$OUTDIR"

echo "Resolving DOI..."

# 1. Resolve DOI and save landing page
curl -L "$DOI_URL" -o page.html

# 2. Extract the View Contents URL
VIEW_URL=$(grep -oP 'id="viewContents"[^>]*href="\K[^"]+' page.html)

if [ -z "$VIEW_URL" ]; then
    echo "Could not find View Contents URL."
    exit 1
fi

echo "View contents URL:"
echo "$VIEW_URL"

# 3. Build original/ URL
ORIGINAL_URL="${VIEW_URL%/}/original/"

echo "Downloading from:"
echo "$ORIGINAL_URL"

# 4. Download everything inside original/
wget -r -np -nH --cut-dirs=3 \
    --reject "index.html*" \
    -P "$OUTDIR" \
    "$ORIGINAL_URL"

echo "Download completed."
echo "Files saved in: $OUTDIR"
