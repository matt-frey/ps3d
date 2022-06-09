#!/bin/bash

dir="$(dirname $0)"

echo "Run all unit tests ..."

for file in "$dir"/*; do
    if [[ "$file" != "$0" ]]; then
        $file
    fi
done

echo "Done."
